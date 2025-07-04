# ================================= GENERAL MODEL ANALYTICS ============================================
# Analyses various model-related quantities towards discerning some insight
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R
#   - 5b(i).CoxDiscreteTime_Advanced.R
#   - 5b(ii).CoxDiscreteTime_Basic.R

# -- Inputs:
#   - datCredit_real | Prepared from script 2f.
#
# -- Outputs:
#   - Cumulative baseline hazard rates by performance spells (2 different groupings)
#   - Event probabilities
#   - Hazard rates
#   - Kaplan-Meyer analysis by performance spells
#   - datSurv objects | Respective to each setting, containiing survival, cumulative hezard, & event probabilities
# ------------------------------------------------------------------------------------------------------




# ----------------- 1. Load & prepare data for analysis

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# - Weigh default cases heavier. as determined interactively based on calibration success (script 6e)
datCredit_train[, Weight := ifelse(DefaultStatus1==1,10,1)]
datCredit_valid[, Weight := ifelse(DefaultStatus1==1,10,1)] # for merging purposes




# ----------------- 2. Fit a discrete-time hazard model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic discrete-time hazard model
# - Initialize variables
vars_basic <- c("-1", "Time_Binned", "log(TimeInPerfSpell):PerfSpell_Num_binned",
                "Arrears", "InterestRate_Nom", "M_Inflation_Growth_6")

# - Fit discrete-time hazard model with selected variables
modLR_basic <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars_basic, collapse = " + "))),
                    data=datCredit_train, family="binomial", weights = Weight)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced discrete-time hazard model
# - Initialize variables
vars <- c("-1", "Time_Binned*PerfSpell_Num_binned", #"log(TimeInPerfSpell):PerfSpell_Num_binned",
          "g0_Delinq_SD_4", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom",
          "M_Inflation_Growth_6","M_DTI_Growth")

# - Fit discrete-time hazard model with selected variables
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial", weights = Weight)




# ----------------- 3. Comparative study of various plots of interesting survival quantities

# --- Preliminaries
# - Create pointer to the appropriate data object 
#datCredit <- rbind(datCredit_train_PWPST, datCredit_valid_PWPST)
datCredit <- datCredit_valid_PWPST

# --- Handle left-truncated spells by adding a starting record 
# This is necessary for calculating certain survival quantities later
datAdd <- subset(datCredit, Counter == 1 & TimeInPerfSpell > 1)
datAdd[, TimeInPerfSpell := TimeInPerfSpell-1]
datAdd[, Counter := 0]
datCredit <- rbind(datCredit, datAdd)
setorder(datCredit, PerfSpell_Key, TimeInPerfSpell) # Preparing for survival model scoring


# --- Calculate survival quantities of interest
# Predict hazard h(t) = P(T=t | T>= t) in discrete-time
datCredit[!is.na(PerfSpell_Num), Hazard_adv := predict(modLR, newdata=.SD[], type = "response")]
datCredit[!is.na(PerfSpell_Num), Hazard_bas := predict(modLR_basic, newdata=.SD[], type = "response")]
# Derive survival probability S(t) = \prod ( 1- hazard)
datCredit[!is.na(PerfSpell_Num), Survival_adv := cumprod(1-Hazard_adv), by=list(PerfSpell_Key)]
datCredit[!is.na(PerfSpell_Num), Survival_bas := cumprod(1-Hazard_bas), by=list(PerfSpell_Key)]
# Derive discrete density, or event probability f(t) = S(t-1) . h(t)
datCredit[!is.na(PerfSpell_Num), EventRate_adv := shift(Survival_adv, type="lag", n=1, fill=1) - Survival_adv, by=list(PerfSpell_Key)]
datCredit[!is.na(PerfSpell_Num), EventRate_bas := shift(Survival_bas, type="lag", n=1, fill=1) - Survival_bas, by=list(PerfSpell_Key)]

# - Remove added rows
datCredit <- subset(datCredit, Counter > 0) %>% setorder(LoanID, Date)


# --- Aggregate to time in performing spell and calculate various survival quantities
datAggr <- datCredit[, list(Defaults = sum(DefaultStatus1), AtRisk = .N, Hazard = sum(DefaultStatus1) / .N,
                            EventRate_bas = mean(EventRate_adv, na.rm=T)), 
                     by=list(TimeInPerfSpell)][order(TimeInPerfSpell),]
# - Calculate more survival quantities
datAggr[, Survival := 1]
datAggr[, Survival := cumprod(shift(Survival, type="lag",n=1,fill=1)*(1-Hazard))] # Kaplan-Meier
datAggr[, CumulLifeDist := 1-Survival]
datAggr[, EventRate := shift(Survival,type="lag", n=1,fill=1) - Survival]
datAggr[, CumulHazard := cumsum(Hazard)]
datAggr[, CumulEventRate_bas := cumsum(EventRate_bas)]
# - Cumulative default-related aggregates
datAggr[, CumulDefaults := cumsum(Defaults)]
datAggr[, Pop_Init := AtRisk[1]]
datAggr[, CumulDefRate := CumulDefaults / Pop_Init]
datAggr[, DefRate_Marg := CumulDefRate - shift(CumulDefRate,type="lag",n=1, fill=NA)]


# - Plots: cumulative default rate
plot(datAggr[TimeInPerfSpell <= 300, CumulDefRate], type="b")
plot(datAggr[TimeInPerfSpell <= 300, DefRate_Marg], type="b")
# - Plots: survival-related
plot(datAggr[TimeInPerfSpell <= 300, CumulEventRate_bas], type="b")
plot(datAggr[, CumulEventRate_bas], type="b")
plot(datAggr[TimeInPerfSpell <= 300, Hazard], type="b")
plot(datAggr[TimeInPerfSpell <= 300, EventRate], type="b")
# - Plots: survival distributions
plot(datAggr[TimeInPerfSpell <= 300, Survival], type="b")
plot(datAggr[TimeInPerfSpell <= 300, CumulLifeDist], type="b")


### RESULTS: It seems there are fundamental differences in how one might derive the term-structure of default risk,
# at least ito of its derivation from the cumulative default rate [DefRate_Marg], vs die survival-related [EventRate].
# While both quantitites seem to have the same range of values, and even the same locale at earlier ages,
# they differ markedly for later ages. If we take [EventRate] as the true quantity, then [DefRate_Marg] vastly
# underestimates default risk at later spell ages.