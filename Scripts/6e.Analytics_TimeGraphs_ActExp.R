# ======================= TIME GRAPHS OF EVENT RATES: ACTUAL VS EXPECTED ===============================
# Calculate and compare time graphs of the actual vs expected resolution & 12-month default rates
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
#   - 2f.Data_Fusion1.R
#   - 3b.Data_Fusion2_PWP_ST.R
#   - 5a(i).CoxPropHaz_PWP_Advanced.R
#   - 5a(i).CoxPropHaz_PWP_Basic.R
#   - 5b(i).CoxDiscreteTime_Advanced.R
#   - 5b(ii).CoxDiscreteTime_Basic.R

# -- Inputs:
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b
#
# -- Outputs:
#   - <Analytics> | Graphs
# ------------------------------------------------------------------------------------------------------






# ----------------- 1. Load & prepare data for analysis

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# - Weigh cases heavier during known downturns
datCredit_train[, Weight := ifelse( #(Date >= as.Date("2007-06-30") & Date <= as.Date("2009-03-31")) | 
                                     #(Date >= as.Date("2019-06-30") & Date <= as.Date("2020-05-31")) | 
                                     DefaultStatus1==1,
                                   10,1)]


# ----------------- 2a. Fit a Cox regression model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic Cox-regression model
# - Initialize variables
vecVars_PWPST_bas <- c("Arrears", "InterestRate_Nom")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_basic <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                           paste(vecVars_PWPST_bas,collapse=" + "), 
                                           " + strata(PerfSpell_Num_binned)")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST_basic); AIC(cox_PWPST_basic); concordance(cox_PWPST_basic)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced Cox-regression model
# - Initialize variables
vecVars_PWPST_adv <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_adv <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                         paste(vecVars_PWPST_adv,collapse=" + "), 
                                         " + strata(PerfSpell_Num_binned)")),
                       id=PerfSpell_Key, datCredit_train, ties="efron")
summary(cox_PWPST_adv); AIC(cox_PWPST_adv); concordance(cox_PWPST_adv)




# ----------------- 2b. Fit a discrete-time hazard model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic discrete-time hazard model
# - Initialize variables
vars_basic <- c("-1", "Time_Binned", "log(TimeInPerfSpell):PerfSpell_Num_binned",
                "Arrears", "InterestRate_Nom", "M_Inflation_Growth_6")

# - Fit discrete-time hazard model with selected variables
modLR_basic <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars_basic, collapse = " + "))),
                    data=datCredit_train, family="binomial", weights = Weight)
summary(modLR_basic)

# Robust (sandwich) standard errors
robust_se <- vcovHC(modLR_basic, type="HC0")

# - Summary with robust SEs
coeftest(modLR_basic, vcov.=robust_se)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced discrete-time hazard model
# - Initialize variables
vars <- c("-1", "Time_Binned*PerfSpell_Num_binned", #"log(TimeInPerfSpell):PerfSpell_Num_binned",
          "g0_Delinq_SD_4", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom", "BalanceToPrincipal",
          "M_Inflation_Growth_6","M_DTI_Growth")

# - Fit discrete-time hazard model with selected variables
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial", weights = Weight)
summary(modLR)

# Robust (sandwich) standard errors
robust_se_adv <- vcovHC(modLR, type="HC0")

# - Summary with robust SEs
coeftest(modLR, vcov.=robust_se_adv)



# ----------------- 3. Constructing expected term-structure of default risk | Discrete-time hazard models

# --- Preliminaries
# - Create pointer to the appropriate data object 
datCredit <- rbind(datCredit_train_PWPST, datCredit_valid_PWPST)


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


# --- Create 12-month PD by summing event probabilities across a rolling window
# NOTE: This step deliberately spans both performing and default spells
# NOTE: Need to specify a (k+1)-window for the "frollapply()" function, e.g., a 12-month outcome implies 13 elements
# Uses the custom function "imputLastKnown" defined in script 0 to fill in NAs caused towards the end of the window
maxDate <- max(datCredit[,Date], na.rm=T) - years(1) # Dates larger than maxDate do not have 12-month default because of the end of the sampling window
datCredit[, PD_12_bas := ifelse(Date<=maxDate,imputeLastKnown(frollapply(
  x=EventRate_bas, n=13, align="left", FUN=function(x){
    sum(x,na.rm=T)
  })),NA), by=list(LoanID)]
datCredit[, PD_12_adv := ifelse(Date<=maxDate,imputeLastKnown(frollapply(
  x=EventRate_adv, n=13, align="left", FUN=function(x){
    sum(x,na.rm=T)
  })),NA), by=list(LoanID)]

# Assign event probability 1 to defaulted cases
# NOTE: This does NOT influence the eventual default rate time graphs
datCredit[is.na(EventRate_bas), PD_12_bas:=1]
datCredit[is.na(EventRate_adv), PD_12_adv:=1]


# --- Graph 12-month expected vs actual default rate

# - Aggregate to reporting time level
datAggr_act <- datCredit[DefaultStatus1==0, list(EventRate = mean(DefaultStatus1_lead_12_max,na.rm=T), Type="a_Act"),, by=list(Date)]
datAggr_bas <- datCredit[DefaultStatus1==0, list(EventRate = mean(PD_12_bas,na.rm=T), Type="b_ExpBasic"),, by=list(Date)]
datAggr_adv <- datCredit[DefaultStatus1==0, list(EventRate = mean(PD_12_adv,na.rm=T), Type="c_ExpAdvanced"),, by=list(Date)]
datAggr <- rbind(datAggr_act, datAggr_bas, datAggr_adv)[Date<=maxDate,]

# - Annotations
(MAE_bas <- mean(abs(datAggr_act$EventRate - datAggr_bas$EventRate), na.rm=T))
(MAE_adv <- mean(abs(datAggr_act$EventRate - datAggr_adv$EventRate), na.rm=T))

# - Graphing parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(9, "Set1")[c(3,1,2)]
vLabel <- c("a_Act" = "Actual/Empirical rate", 
            "b_ExpBasic" = "Expected rate: Basic model",
            "c_ExpAdvanced" = "Expected rate: Advanced model")

(gPlot <- ggplot(datAggr, aes(x=Date, y=EventRate, group=Type)) + theme_minimal() + 
    labs(x="Reporting time (ccyymm)", y="12-month default rate (%)") +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Type, linetype=Type), size=0.5) + 
    #geom_point(aes(colour=Type, shape=Type), size=1) + 
    #annotations
    annotate(geom="text", y=0.05, x=as.Date("2015-04-30"), family=chosenFont, size=3,
             label=paste0("MAE (basic): ", round(MAE_bas*100,digits=3), "%")) + 
    annotate(geom="text", y=0.04, x=as.Date("2014-12-31"), family=chosenFont, size=3,
             label=paste0("MAE (advanced): ", round(MAE_adv*100,digits=3), "%")) + 
    # scale options
    scale_colour_manual(name="", values=vCol, labels=vLabel) + 
    scale_shape_discrete(name="", labels=vLabel) + 
    scale_linetype_discrete(name="", labels=vLabel) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# Save graph
dpi <- 200
ggsave(gPlot, file=paste0(genFigPath, "DefaultRate_ActVsExp_CoxDisc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")






########################## SCRATCH ####################################
# --- cumulative default probability
datAggr <- datCredit[, list(Defaults = sum(DefaultStatus1), AtRisk = .N, Hazard = sum(DefaultStatus1) / .N,
                            EventRate_bas = mean(EventRate_adv, na.rm=T)), 
                     by=list(TimeInPerfSpell)][order(TimeInPerfSpell),]
# - Survival quantities
datAggr[, Survival := 1]
datAggr[, Survival := cumprod(shift(Survival, type="lag",n=1,fill=1)*(1-Hazard))] # Kaplan-Meier
datAggr[, CumulLifeDist := 1-Survival]
datAggr[, EventRate := shift(Survival,type="lag", n=1,fill=1) - Survival]
datAggr[, CumulHazard := cumsum(Hazard)]
datAggr[, CumulEventRate_bas := cumsum(EventRate_bas)]
# - Cumulative default-related aggregates
datAggr[, CumulDefaults := cumsum(Defaults)]
#datAggr[, Pop_Init := length(unique(datCredit$LoanID))]
datAggr[, Pop_Init := AtRisk[1]]
datAggr[, CumulDefRate := CumulDefaults / Pop_Init]
datAggr[, DefRate_Marg := CumulDefRate - shift(CumulDefRate,type="lag",n=1, fill=NA)]
# Plots
plot(datAggr[TimeInPerfSpell <= 300, CumulDefRate], type="b")
plot(datAggr[TimeInPerfSpell <= 300, CumulEventRate_bas], type="b")

plot(datAggr[TimeInPerfSpell <= 300, DefRate_Marg], type="b")
plot(datAggr[TimeInPerfSpell <= 300, Hazard], type="b")
plot(datAggr[TimeInPerfSpell <= 300, EventRate], type="b")
plot(datAggr[TimeInPerfSpell <= 300, Survival], type="b")
plot(datAggr[TimeInPerfSpell <= 300, CumulLifeDist], type="b")

