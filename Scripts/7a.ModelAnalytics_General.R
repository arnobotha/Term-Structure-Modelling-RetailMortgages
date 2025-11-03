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

# - Create subset of performing spells towards modelling each spell as a single record | Classical PD-modelling
datTrain_classic <- subset(datCredit_train, PerfSpell_Counter==1)
datValid_classic <- subset(datCredit_valid, PerfSpell_Counter==1)

# - Assign target/response field for the classical PD-model
datTrain_classic[, Defaulted := ifelse(PerfSpellResol_Type_Hist=="Defaulted", 1, 0)]
datValid_classic[, Defaulted := ifelse(PerfSpellResol_Type_Hist=="Defaulted", 1, 0)]

# - Create a copy of PerfSpell_Age to serve as an input into the classical PD-model
datTrain_classic[, PerfSpell_Age2 := PerfSpell_Age]
datValid_classic[, PerfSpell_Age2 := PerfSpell_Age]


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


# ------ Classical Logistic Regression PD-model

# - Initialize variables
vars_classic <- c("PerfSpell_Age2*PerfSpell_Num_binned",
          "g0_Delinq_SD_4", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom",
          "M_Inflation_Growth_6","M_DTI_Growth")

# - Fit logistic regression model onto entire spell, similar to application credit scoring with an open-ended outcome window
modLR_classic <- glm( as.formula(paste("Defaulted ~", paste(vars_classic, collapse = " + "))),
              data=datTrain_classic, family="binomial")
summary(modLR_classic)
roc(response=datTrain_classic$Defaulted, 
    predictor=predict(modLR_classic, datTrain_classic, type="response")) # 80% AUC




# ----------------- 3. Comparative study of various plots of interesting survival quantities

# --- Preliminaries
# - Create pointer to the appropriate data object 
#datCredit <- rbind(datCredit_train, datCredit_valid)
datCredit <- datCredit_valid

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

# - Score using classical PD-model each instance of TimeInPerfSpell as PerfSpell_Age
datCredit[, PerfSpell_Age2 := TimeInPerfSpell]
datCredit[!is.na(PerfSpell_Num), Hazard_PD := predict(modLR_classic, newdata=.SD[], type = "response")]
datCredit[!is.na(PerfSpell_Num), Survival_PD := cumprod(1-Hazard_PD), by=list(PerfSpell_Key)]
datCredit[!is.na(PerfSpell_Num), EventRate_PD := shift(Survival_PD, type="lag", n=1, fill=1) - Survival_PD, by=list(PerfSpell_Key)]



# --- 1a. Aggregate to time in performing spell and calculate various survival quantities: TimeInPerfSpell
datAggr <- datCredit[, list(Defaults = sum(DefaultStatus1), AtRisk = .N, Hazard = mean(DefaultStatus1),
                            EventRate_PD = mean(EventRate_PD, na.rm=T),
                            EventRate_adv = mean(EventRate_adv, na.rm=T),
                            Hazard_adv = mean(Hazard_adv, na.rm=T)), 
                     by=list(TimeInPerfSpell)][order(TimeInPerfSpell),]
# - Calculate more survival quantities
datAggr[, Survival := 1]
datAggr[, Survival := cumprod(shift(Survival, type="lag",n=1,fill=1)*(1-Hazard))] # Kaplan-Meier
datAggr[, CumulLifeDist := 1-Survival]
datAggr[, CumulHazard := cumsum(Hazard)]
datAggr[, CumulHazard_adv := cumsum(Hazard_adv)]
datAggr[, EventRate := shift(Survival,type="lag", n=1,fill=1) - Survival]
# - Cumulative default-related aggregates
datAggr[, CumulDefaults := cumsum(Defaults)]
datAggr[, Pop_Init := AtRisk[1]]
datAggr[, CumulDefRate := CumulDefaults / Pop_Init]
datAggr[, DefRate_Marg := CumulDefRate - shift(CumulDefRate,type="lag",n=1, fill=NA)]


# -- Plots
# - Plot: Base
plot(datAggr[TimeInPerfSpell <= 300, AtRisk], type="b")
# - Plots: survival distributions
plot(datAggr[TimeInPerfSpell <= 300, Survival], type="b")
plot(datAggr[TimeInPerfSpell <= 300, CumulLifeDist], type="b")
plot(datAggr[, CumulLifeDist], type="b")
# - Plots: survival-related
plot(datAggr[TimeInPerfSpell <= 300, CumulHazard], type="b")
lines(datAggr[TimeInPerfSpell <= 300, CumulHazard_adv], type="b", col="red")
plot(datAggr[TimeInPerfSpell <= 300, Hazard], type="b")
# Event rates
plot(datAggr[TimeInPerfSpell <= 300, EventRate], type="b")
lines(datAggr[TimeInPerfSpell <= 300, EventRate_adv], type="b", col="red")
lines(datAggr[TimeInPerfSpell <= 300, EventRate_PD], type="b", col="cyan")
# Cumulative lifetime distribution F(t)
plot(datAggr[TimeInPerfSpell <= 300, 1- exp(-CumulHazard_adv)], type="b", col="red")
lines(datAggr[TimeInPerfSpell <= 300, 1- exp(-CumulHazard)], type="b") 
# - Plots: cumulative default rate
plot(datAggr[TimeInPerfSpell <= 300, CumulDefRate], type="b")
plot(datAggr[TimeInPerfSpell <= 300, DefRate_Marg], type="b")

### RESULTS: It seems there are fundamental differences in how one might derive the term-structure of default risk,
# at least ito of its derivation from the cumulative default rate [DefRate_Marg], vs die survival-related [EventRate].
# While both quantities seem to have the same range of values, and even the same locale at earlier times,
# they differ markedly for later times If we take [EventRate] as the true quantity, then [DefRate_Marg] vastly
# underestimates default risk at later spell times



# --- 1b. Aggregate to time in performing spell and calculate various survival quantities: PerfSpell_Age
datAggr2 <- datCredit[PerfSpell_Counter==1, list(Defaults = sum(PerfSpellResol_Type_Hist=="Defaulted"), N = .N,
                                                 EventRate_bas = mean(EventRate_bas, na.rm=T), 
                                                 EventRate_adv = mean(EventRate_adv, na.rm=T), 
                                                 EventRate_PD = mean(EventRate_PD, na.rm=T),
                                                 Hazard_bas = mean(Hazard_bas, na.rm=T),
                                                 Hazard_adv = mean(Hazard_adv, na.rm=T),
                                                 EventRate_adv_med = median(EventRate_adv, na.rm=T), 
                                                 Hazard_adv_med = median(Hazard_adv, na.rm=T)),                      
                     by=list(PerfSpell_Age)][order(PerfSpell_Age),]
# - Calculate base survival quantities
datAggr2[, AtRisk := sum(N) - cumsum(Defaults)]
datAggr2[, Hazard := Defaults / AtRisk]
# - Calculate more survival quantities
datAggr2[, Survival := 1]
datAggr2[, Survival := cumprod(shift(Survival, type="lag",n=1,fill=1)*(1-Hazard))] # Kaplan-Meier
datAggr2[, CumulLifeDist := 1-Survival]
datAggr2[, CumulHazard := cumsum(Hazard)]
datAggr2[, CumulHazard_bas := cumsum(Hazard_bas)]
datAggr2[, CumulHazard_adv := cumsum(Hazard_adv)]
datAggr2[, CumulHazard_adv_med := cumsum(Hazard_adv_med)]
datAggr2[, EventRate := shift(Survival,type="lag", n=1,fill=1) - Survival]
# - Cumulative default-related aggregates
datAggr2[, CumulDefaults := cumsum(Defaults)]
datAggr2[, Pop_Init := AtRisk[1]]
datAggr2[, CumulDefRate := CumulDefaults / Pop_Init]
datAggr2[, DefRate_Marg := CumulDefRate - shift(CumulDefRate,type="lag",n=1, fill=NA)]


# -- Plots
# - Plot: Base
plot(datAggr2[PerfSpell_Age <= 300, N], type="b")
plot(datAggr2[PerfSpell_Age <= 300, AtRisk], type="b")
# - Plots: survival distributions
plot(datAggr2[PerfSpell_Age <= 300, Survival], type="b")
plot(datAggr2[PerfSpell_Age <= 300, CumulLifeDist], type="b")
plot(datAggr2[, CumulLifeDist], type="b")
# - Plots: survival-related
plot(datAggr2[PerfSpell_Age <= 300, CumulHazard], type="b")
lines(datAggr2[PerfSpell_Age <= 300, CumulHazard_adv_med], type="b", col="green")
plot(datAggr2[PerfSpell_Age <= 300, CumulHazard_adv], type="b", col="red")
plot(datAggr2[PerfSpell_Age <= 300, CumulHazard_bas], type="b", col="blue")
plot(datAggr2[PerfSpell_Age <= 300, Hazard], type="b")
# Event rate
plot(datAggr2[PerfSpell_Age <= 300, EventRate], type="b")
lines(datAggr2[PerfSpell_Age <= 300, EventRate_adv], type="b", col="red")
lines(datAggr2[PerfSpell_Age <= 300, EventRate_adv_med], type="b", col="green")
plot(datAggr2[PerfSpell_Age <= 300, EventRate_bas], type="b", col="blue") # very far away from empirical
plot(datAggr2[PerfSpell_Age <= 300, EventRate_PD], type="b", col="cyan") # even further away from empirical
# Cumulative lifetime distribution F(t)
plot(datAggr2[PerfSpell_Age <= 300, 1- exp(-CumulHazard_adv)], type="b", col="red")
lines(datAggr2[PerfSpell_Age <= 300, 1- exp(-CumulHazard)], type="b") 
lines(datAggr2[PerfSpell_Age <= 300, 1- exp(-CumulHazard_adv_med)], type="b", col="green")
# - Plots: cumulative default rate
plot(datAggr2[PerfSpell_Age <= 300, CumulDefRate], type="b") # same as CumulLifeDist
plot(datAggr2[PerfSpell_Age <= 300, DefRate_Marg], type="b") # same as EventRate

# --- Distribution analysis at specific points in spell duration
ptInspect <- 12
datInspect <- datCredit[PerfSpell_Counter==1 & PerfSpell_Age==ptInspect,]
describe(datInspect$EventRate_adv); hist(datInspect$EventRate_adv, breaks="FD")
hist(datInspect[EventRate_adv <= 0.003, EventRate_adv], breaks="FD")
### RESULTS: Highly right-skewed distribution with extreme outliers. 
# For t=12, majority of distribution (5%-95%) lies within [0.0000009612 ; 0.0028248805 ]


### RESULTS: When aggregated by spell age, the difference between [DefRate_Marg] and [EventRate] resolves to 0.
# The same happens with [CumulLifeDist] and [CumulDefRate]
# When swapping the mean with the median, the aggregates switches from over-prediction to under-prediction.
# Even though a strong case can be made for using the median, it is no longer the expected value E(x), and
# therefore disqualifies itself.
