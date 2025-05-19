# ======================================= INPUT SPACE: Advanced  =======================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the TPWPST model
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

# -- Inputs:
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b
#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------





# ------ 0. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# - remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Empty model | Stepwise forward selection procedure
cox_PWPST_empty <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ 1")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST_empty)



# ------ 1. Delinquency-themed variables

# ------ 1.1 Which time window length is the best in calculating delinquency volatility?

# - Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: g0_Delinq_SD_5, g0_Delinq_SD_6, g0_Delinq_SD_9, g0_Delinq_SD_4 
# Best Harrell's C-statistics: g0_Delinq_SD_4, g0_Delinq_SD_5, g0_Delinq_SD_6
# All variables are statistically significant on their own

### CONCLUSION: Decreasing trend in Harrell's c across window length (earlier is better). 
# Best variable: g0_Delinq_SD_4. Variable has >90% value in Harrell's c



# ------ 1.2 Which lag order is the best in calculating the portfolio-level fraction of the non-defaulted proportion with any delinquency?

# - Initialize variables to be tested
vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Any_Aggr_Prop_Lag_1", "g0_Delinq_Any_Aggr_Prop_Lag_2",
          "g0_Delinq_Any_Aggr_Prop_Lag_3", "g0_Delinq_Any_Aggr_Prop_Lag_4", "g0_Delinq_Any_Aggr_Prop_Lag_5",
          "g0_Delinq_Any_Aggr_Prop_Lag_6", "g0_Delinq_Any_Aggr_Prop_Lag_9", "g0_Delinq_Any_Aggr_Prop_Lag_12" )

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath) 
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: g0_Delinq_Any_Aggr_Prop_Lag_2, g0_Delinq_Any_Aggr_Prop_Lag_1 , g0_Delinq_Any_Aggr_Prop_Lag_3
# Best Harrell's C-statistics: g0_Delinq_Any_Aggr_Prop_Lag_3, g0_Delinq_Any_Aggr_Prop_Lag_2, g0_Delinq_Any_Aggr_Prop_Lag_1
# All variables are statistically significant on their own

### CONCLUSION: Decreasing trend in Harrell's c across lag order (earlier is better), though differences are slight at best
# Best variable: g0_Delinq_Any_Aggr_Prop_Lag_3. Variable has >60% value in Harrell's c



# ------ 1.3 Which lag order is the best in calculating the portfolio-level fraction of defaulted accounts?

# - Initialize variables to be tested
vars <- c("DefaultStatus1_Aggr_Prop", "DefaultStatus1_Aggr_Prop_Lag_1", "DefaultStatus1_Aggr_Prop_Lag_2",
          "DefaultStatus1_Aggr_Prop_Lag_3", "DefaultStatus1_Aggr_Prop_Lag_4", "DefaultStatus1_Aggr_Prop_Lag_5",
          "DefaultStatus1_Aggr_Prop_Lag_6", "DefaultStatus1_Aggr_Prop_Lag_9", "DefaultStatus1_Aggr_Prop_Lag_12" )

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath) 
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: DefaultStatus1_Aggr_Prop, DefaultStatus1_Aggr_Prop_Lag_1, DefaultStatus1_Aggr_Prop_Lag_2, DefaultStatus1_Aggr_Prop_Lag_3 
# Best Harrell's C-statistics: DefaultStatus1_Aggr_Prop, DefaultStatus1_Aggr_Prop_Lag_1, DefaultStatus1_Aggr_Prop_Lag_2, DefaultStatus1_Aggr_Prop_Lag_3
# All variables are statistically significant on their own

### CONCLUSION: Decreasing trend in Harrell's c across lag order (earlier is better)
# Best variable: g0_Delinq_Any_Aggr_Prop_Lag_3. Variable has >50% value in Harrell's c



# ------ 1.4 How do other portfolio-level delinquency-themed variables fare as single-factor models?

# - Initialize variables to be tested
vars <- c("g0_Delinq_Ave", "ArrearsToBalance_Aggr_Prop", "CuringEvents_Aggr_Prop" )

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath) 
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: g0_Delinq_Ave, ArrearsToBalance_Aggr_Prop, CuringEvents_Aggr_Prop
# Best Harrell's C-statistics: g0_Delinq_Ave, ArrearsToBalance_Aggr_Prop, CuringEvents_Aggr_Prop
# All variables are statistically significant on their own

### CONCLUSION: 
# Best variable: g0_Delinq_Ave, ArrearsToBalance_Aggr_Prop. Variable has >60% value in Harrell's c



# ------ 1.5 How do account-level delinquency-themed variables fare as single-factor models?

# - Feature engineering
datCredit_train[, g0_Delinq_Lag_1 := shift(g0_Delinq,fill=0),by=LoanID]
datCredit_valid[, g0_Delinq_Lag_1 := shift(g0_Delinq,fill=0),by=LoanID]

# - Initialize variables to be tested
vars <- c("TimeInPerfSpell", "g0_Delinq_fac", "g0_Delinq", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3_Change_Ind",
          "g0_Delinq_Num", "TimeInDelinqState", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "PerfSpell_g0_Delinq_Num", "Arrears")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath) 
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: TimeInDelinqState, slc_acct_roll_ever_24_imputed_mean, Arrears, PerfSpell_g0_Delinq_Num, g0_Delinq_Num
# Best Harrell's C-statistics: TimeInDelinqState, Arrears, PerfSpell_g0_Delinq_Num, g0_Delinq_Num, slc_acct_roll_ever_24_imputed_mean
# Following were not statistically significant or led to unstable models: g0_Delinq, g0_Delinq_fac

### CONCLUSION: 
# Best variable: TimeInDelinqState, PerfSpell_g0_Delinq_Num, g0_Delinq_Num, slc_acct_roll_ever_24_imputed_mean
# Variable has >90% value in Harrell's c



# ------ 1.6 Combining insights: delinquency-themed variables

# - Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_Any_Aggr_Prop_Lag_3", "DefaultStatus1_Aggr_Prop", "g0_Delinq_Ave", "ArrearsToBalance_Aggr_Prop",
          "TimeInDelinqState", "PerfSpell_g0_Delinq_Num", "g0_Delinq_Num", "slc_acct_roll_ever_24_imputed_mean", "Arrears")

# - Correlation analysis towards obtaining clusters of correlated variables
corrAnalysis(datCredit_train, vars, corrThresh = 0.6, method = 'spearman')
### RESULTS: Some variables are expectedly correlated, though the following correlations indicate a choice:
# PerfSpell_g0_Delinq_Num  and  g0_Delinq_Num
# g0_Delinq_Ave  and  ArrearsToBalance_Aggr_Prop and g0_Delinq_Any_Aggr_Prop_Lag_3  

# - Full model | Stepwise forward selection procedure
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST)
# Statistical insignificance detected for: g0_Delinq_Any_Aggr_Prop_Lag_3

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
cox_PWPST_step <- stepAIC(cox_PWPST_empty, scope = list(lower = ~ 1, 
                                                        upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                          direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(cox_PWPST_step); AIC(cox_PWPST_step); concordance(cox_PWPST_step);
proc.time() - ptm # IGNORE: elapsed runtime; 32m
### RESULTS: strata(PerfSpell_Num_binned) could not be included into the stepAIC() function, so the stepwise results
# are only preliminary and not definitive.
# Selected variables include: g0_Delinq_SD_4 + TimeInDelinqState + slc_acct_roll_ever_24_imputed_mean + 
# g0_Delinq_Ave + ArrearsToBalance_Aggr_Prop + DefaultStatus1_Aggr_Prop + g0_Delinq_Num
# AIC: 133259.8; Harrell's c: 0.9975 

# - Domain expertise
# Model is unstable due to some covariates having zero-values for exp(coef)
# Add Arrears for account-level delinquency
# Remove ArrearsToBalance_Aggr_Prop due to its high standard error
# Remove TimeInDelinqState since it renders model as unstable

# - Final variables
vars <- c("g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave",
          "DefaultStatus1_Aggr_Prop", "g0_Delinq_Num", "Arrears")

cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
### RESULTS: All variables are statistically significant
# AIC: 124354.4; Harrell's c: 0.9958    




# ------ 2. Other portfolio-level variables

# ------ 2.1 Which lag order is the best in calculating the median interest rate of the portfolio?

# - Initialize variables to be tested
vars <- c("InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_1", "InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_3", "InterestRate_Margin_Aggr_Med_9")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: InterestRate_Margin_Aggr_Med, InterestRate_Margin_Aggr_Med_1, InterestRate_Margin_Aggr_Med_2
# Best Harrell's C-statistics: InterestRate_Margin_Aggr_Med, InterestRate_Margin_Aggr_Med_1, InterestRate_Margin_Aggr_Med_2
# All variables are statistically significant on their own

### CONCLUSION: Decreasing trend in Harrell's c across lag order (earlier is better). 
# Best variable: InterestRate_Margin_Aggr_Med. Variable has >60% value in Harrell's c



# ------ 2.2 How do other portfolio-level (non-delinquency) variables fare as single-factor models?

# - Initialize variables to be tested
vars <- c("InstalmentToBalance_Aggr_Prop", "AgeToTerm_Aggr_Mean", "PerfSpell_Maturity_Aggr_Mean",
          "CreditLeverage_Aggr", "Ave_Margin_Aggr", "NewLoans_Aggr_Prop")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: AgeToTerm_Aggr_Mean, InstalmentToBalance_Aggr_Prop, PerfSpell_Maturity_Aggr_Mean, NewLoans_Aggr_Prop
# Best Harrell's C-statistics: AgeToTerm_Aggr_Mean, InstalmentToBalance_Aggr_Prop, PerfSpell_Maturity_Aggr_Mean, NewLoans_Aggr_Prop
# All variables are statistically significant on their own

### CONCLUSION: Only certain variables resulted in stable models
# Best variable: AgeToTerm_Aggr_Mean, InstalmentToBalance_Aggr_Prop. Variable has >60% value in Harrell's c



# ------ 2.3 Combining insights: Delinquency-themed and portfolio-level variables

# - Initialize variables to be tested
vars.min <- c("g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave",
              "DefaultStatus1_Aggr_Prop", "g0_Delinq_Num", "Arrears")
vars <- c("AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "PerfSpell_Maturity_Aggr_Mean", "NewLoans_Aggr_Prop")

# - Correlation analysis towards obtaining clusters of correlated variables
corrAnalysis(datCredit_train, vars, corrThresh = 0.6, method = 'spearman')
### RESULTS: Some variables are expectedly correlated, though the following correlations indicate a choice:
# PerfSpell_Maturity_Aggr_Mean  and  NewLoans_Aggr_Prop 

# - Starting model | Stepwise forward selection procedure
cox_PWPST_empty <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars.min,collapse=" + "))),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST_empty); AIC(cox_PWPST_empty); concordance(cox_PWPST_empty)

# - Full/nested model | Stepwise forward selection procedure
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(c(vars,vars.min), collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
# Very large coefficient (and standard errors) for: InstalmentToBalance_Aggr_Prop 
# Statistical insignificance detected for: InstalmentToBalance_Aggr_Prop

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
cox_PWPST_step <- stepAIC(cox_PWPST_empty, scope = list(lower = as.formula(paste("~", paste(vars.min, collapse = " + "))), 
                                                        upper = as.formula(paste("~", paste(c(vars.min, vars), collapse = " + ")))), 
                          direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(cox_PWPST_step); AIC(cox_PWPST_step); concordance(cox_PWPST_step);
proc.time() - ptm # IGNORE: elapsed runtime; 32m
### RESULTS: strata(PerfSpell_Num_binned) could not be included into the stepAIC() function, so the stepwise results
# are only preliminary and not definitive.
# Selected variables include: <minimum variables> + AgeToTerm_Aggr_Mean + 
# PerfSpell_Maturity_Aggr_Mean + NewLoans_Aggr_Prop

# - Final variables
vars <- c("g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave",
          "DefaultStatus1_Aggr_Prop", "g0_Delinq_Num", "Arrears", "AgeToTerm_Aggr_Mean",
          "PerfSpell_Maturity_Aggr_Mean", "NewLoans_Aggr_Prop")

# - Final model
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
### RESULTS: All variables are statistically significant
# AIC: 123595.8; Harrell's c: 0.9956   





# ------ 3. Account-level variables

# ------ 3.1 How do various non-delinquency account-level variables fare as single-factor models?

# - Initialize variables to be tested
vars <- c("Principal_Real", "Principal", "InterestRate_Margin", "LN_TPE", "pmnt_method_grp",
          "Balance_Real", "Balance", "Instalment_Real", "InterestRate_Nom", "AgeToTerm",
          "BalanceToPrincipal", "slc_acct_pre_lim_perc_imputed_med")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: slc_acct_pre_lim_perc_imputed_med + InterestRate_Nom  + Principal + BalanceToPrincipal + InterestRate_Margin + LN_TPE + Balance_Real + Principal_Real + Balance + AgeToTerm  
# Best Harrell's C-statistics: BalanceToPrincipal + pmnt_method_grp + slc_acct_pre_lim_perc_imputed_med + InterestRate_Nom
# Statistical insignificance detected for: pmnt_method_grp, Instalment_Real

### CONCLUSION: Only certain variables resulted in stable models
# Best variable: BalanceToPrincipal + pmnt_method_grp + slc_acct_pre_lim_perc_imputed_med + InterestRate_Nom
# Variables have 64%-77% values in Harrell's c



# ------ 3.2 Combining insights: Delinquency-themed, portfolio-level, and account-level variables

# - Initialize variables to be tested
vars.min <- c("g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave",
              "DefaultStatus1_Aggr_Prop", "g0_Delinq_Num", "Arrears", "AgeToTerm_Aggr_Mean",
              "PerfSpell_Maturity_Aggr_Mean", "NewLoans_Aggr_Prop")
vars <- c("BalanceToPrincipal", "pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med", "InterestRate_Nom")

# - Starting model | Stepwise forward selection procedure
cox_PWPST_empty <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars.min,collapse=" + "))),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST_empty); AIC(cox_PWPST_empty); concordance(cox_PWPST_empty)

# - Full/nested model | Stepwise forward selection procedure
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(c(vars,vars.min), collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
# Very large coefficient (and standard errors) for: NewLoans_Aggr_Prop 
# Statistical insignificance detected for: slc_acct_pre_lim_perc_imputed_med
# AIC: 113470.1. Harrell's C: 0.9985 

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
cox_PWPST_step <- stepAIC(cox_PWPST_empty, scope = list(lower = as.formula(paste("~", paste(vars.min, collapse = " + "))), 
                                                        upper = as.formula(paste("~", paste(c(vars.min, vars), collapse = " + ")))), 
                          direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(cox_PWPST_step); AIC(cox_PWPST_step); concordance(cox_PWPST_step);
proc.time() - ptm # IGNORE: elapsed runtime; 5m
### RESULTS: strata(PerfSpell_Num_binned) could not be included into the stepAIC() function, so the stepwise results
# are only preliminary and not definitive.
# Procedure could not obtain convergence in models and subsequently failed.
# Selected variables include: <minimum variables> + BalanceToPrincipal + pmnt_method_grp + slc_acct_pre_lim_perc_imputed_med + InterestRate_Nom

# - Empty model | Stepwise forward selection procedure
cox_PWPST_empty <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ 1")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
cox_PWPST_step <- stepAIC(cox_PWPST_empty, scope = list(lower = ~ 1, 
                                                        upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                          direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(cox_PWPST_step); AIC(cox_PWPST_step); concordance(cox_PWPST_step);
proc.time() - ptm # IGNORE: elapsed runtime; 7m
### RESULTS: strata(PerfSpell_Num_binned) could not be included into the stepAIC() function, so the stepwise results
# are only preliminary and not definitive.
# All variables were selected and are statistically significant





# ------ 4. Macroeconomic variables

# ------ 4.1 Which lag order is the best for: M_Repo_Rate

# - Initialize variables to be tested
vars <- c("M_Repo_Rate", "M_Repo_Rate_1 ", "M_Repo_Rate_2", "M_Repo_Rate_3", "M_Repo_Rate_6", "M_Repo_Rate_9", "M_Repo_Rate_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: M_Repo_Rate_6, M_Repo_Rate_3, M_Repo_Rate_2, M_Repo_Rate_9, M_Repo_Rate_1
# Best Harrell's C-statistics: M_Repo_Rate_9, M_Repo_Rate_12, M_Repo_Rate_6, M_Repo_Rate_3
# All variables are statistically significant

### CONCLUSION: It seems that slightly longer lags fare better than shorter ones, but only up to point (6-months)
# Best variables: M_Repo_Rate_6, M_Repo_Rate_9
# Variables have ~63% values in Harrell's c



# ------ 4.2 Which lag order is the best for: M_Inflation_Growth

# - Initialize variables to be tested
vars <- c("M_Inflation_Growth", "M_Inflation_Growth_1 ", "M_Inflation_Growth_2", "M_Inflation_Growth_3", 
          "M_Inflation_Growth_6", "M_Inflation_Growth_9", "M_Inflation_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: M_Inflation_Growth_6, M_Inflation_Growth_3, M_Inflation_Growth_2, M_Inflation_Growth_1,
# Best Harrell's C-statistics: M_Inflation_Growth_6, M_Inflation_Growth_3, M_Inflation_Growth_9, M_Inflation_Growth_2
# All variables are statistically significant

### CONCLUSION: Aside from 6-month lag, shorter lag orders are better fitting than longer lag orders
# Best variable: M_Inflation_Growth_6, M_Inflation_Growth_3
# Variables have ~63% values in Harrell's c



# ------ 4.3 Which lag order is the best for: M_RealGDP_Growth

# - Initialize variables to be tested
vars <- c("M_RealGDP_Growth", "M_RealGDP_Growth_1 ", "M_RealGDP_Growth_2", "M_RealGDP_Growth_3", 
          "M_RealGDP_Growth_6", "M_RealGDP_Growth_9", "M_RealGDP_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: M_RealGDP_Growth_12, M_RealGDP_Growth_9, M_RealGDP_Growth_6, M_RealGDP_Growth_3
# Best Harrell's C-statistics: M_RealGDP_Growth_12, M_RealGDP_Growth_9, M_RealGDP_Growth_6, M_RealGDP_Growth_3
# All variables are statistically significant

### CONCLUSION: An increasing trend appears across lag orders (longer is better)
# Best variable: M_RealGDP_Growth_12, M_RealGDP_Growth_9
# Variables have ~62% values in Harrell's c



# ------ 4.4 Which lag order is the best for: M_RealIncome_Growth

# - Initialize variables to be tested
vars <- c("M_RealIncome_Growth", "M_RealIncome_Growth_1 ", "M_RealIncome_Growth_2", "M_RealIncome_Growth_3", 
          "M_RealIncome_Growth_6", "M_RealIncome_Growth_9", "M_RealIncome_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: M_RealIncome_Growth_12, M_RealIncome_Growth_9, M_RealIncome_Growth_6, M_RealIncome_Growth_3
# Best Harrell's C-statistics: M_RealIncome_Growth_12, M_RealIncome_Growth_9, M_RealIncome_Growth_6, M_RealIncome_Growth_3
# All variables are statistically significant

### CONCLUSION: An increasing trend appears across lag orders (longer is better)
# Best variable: M_RealIncome_Growth_12, M_RealIncome_Growth_9
# Variables have ~60% values in Harrell's c



# ------ 4.5 Which lag order is the best for: M_DTI_Growth

# - Initialize variables to be tested
vars <- c("M_DTI_Growth", "M_DTI_Growth_1 ", "M_DTI_Growth_2", "M_DTI_Growth_3", 
          "M_DTI_Growth_6", "M_DTI_Growth_9", "M_DTI_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: M_DTI_Growth, M_DTI_Growth_1, M_DTI_Growth_2, M_DTI_Growth_3, M_DTI_Growth_9
# Best Harrell's C-statistics: M_DTI_Growth_9, M_DTI_Growth_6, M_DTI_Growth_1, M_DTI_Growth, 
# All variables are statistically significant

### CONCLUSION: A decreasing trend appears across lag orders (shorter is better), but the reverse is true for Harrell's c
# Best variable: M_DTI_Growth_9, M_DTI_Growth_6, M_DTI_Growth_1
# Variables have ~64% values in Harrell's c



# ------ 4.6 Which lag order is the best for: M_Emp_Growth

# - Initialize variables to be tested
vars <- c("M_Emp_Growth", "M_Emp_Growth_1 ", "M_Emp_Growth_2", "M_Emp_Growth_3", 
          "M_Emp_Growth_6", "M_Emp_Growth_9", "M_Emp_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: M_Emp_Growth_12, M_Emp_Growth_9, M_Emp_Growth_6, M_Emp_Growth_3
# Best Harrell's C-statistics: M_Emp_Growth_12, M_Emp_Growth_9, M_Emp_Growth_6, M_Emp_Growth_3
# Statistical insignificance detected: M_Emp_Growth 

### CONCLUSION: An increasing trend appears across lag orders (longer is better)
# Best variable: M_Emp_Growth_12, M_Emp_Growth_9
# Variables have ~60% values in Harrell's c



# ------ 4.7 Combining insights: Macroeconomic variables

# - Initialize variables to be tested
vars <- c("M_Repo_Rate_6", "M_Repo_Rate_9", "M_Inflation_Growth_6", "M_Inflation_Growth_3",
          "M_RealGDP_Growth_12", "M_RealGDP_Growth_9", "M_DTI_Growth_9", "M_DTI_Growth_6", "M_DTI_Growth_1",
          "M_Emp_Growth_12", "M_Emp_Growth_9", "M_RealIncome_Growth_12", "M_RealIncome_Growth_9" )

# - Empty model | Stepwise forward selection procedure
cox_PWPST_empty <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ 1")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
cox_PWPST_step <- stepAIC(cox_PWPST_empty, scope = list(lower = ~ 1, 
                                                        upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                          direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(cox_PWPST_step); AIC(cox_PWPST_step); concordance(cox_PWPST_step);
proc.time() - ptm # IGNORE: elapsed runtime; 12m
### RESULTS: strata(PerfSpell_Num_binned) could not be included into the stepAIC() function, so the stepwise results
# are only preliminary and not definitive.
# Selected variables: M_DTI_Growth_9 + M_Inflation_Growth_6 + M_Repo_Rate_6



# ------ 4.8 Combining insights: Delinquency-themed, portfolio-level, account-level, and macroeconomic variables

# - Initialize variables to be tested
# No minimum specified deliberately to test domain expertise
vars <- c("g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave",
          "DefaultStatus1_Aggr_Prop", "g0_Delinq_Num", "Arrears", "AgeToTerm_Aggr_Mean",
          "PerfSpell_Maturity_Aggr_Mean", "NewLoans_Aggr_Prop",
          "BalanceToPrincipal", "pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med", "InterestRate_Nom",
          "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
cox_PWPST_step <- stepAIC(cox_PWPST_empty, scope = list(lower = ~ 1, 
                                                        upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                          direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(cox_PWPST_step); AIC(cox_PWPST_step); concordance(cox_PWPST_step);
proc.time() - ptm # IGNORE: elapsed runtime; 167m
### RESULTS: strata(PerfSpell_Num_binned) could not be included into the stepAIC() function, so the stepwise results
# are only preliminary and not definitive.
# Selected variables: g0_Delinq_SD_4 + TimeInDelinqState + slc_acct_roll_ever_24_imputed_mean + g0_Delinq_Ave
#   g0_Delinq_Num + Arrears + AgeToTerm_Aggr_Mean
#   PerfSpell_Maturity_Aggr_Mean + BalanceToPrincipal + pmnt_method_grp + InterestRate_Nom
#   M_DTI_Growth_9 + M_Inflation_Growth_6
# Deselected ones: DefaultStatus1_Aggr_Prop, slc_acct_pre_lim_perc_imputed_med, NewLoans_Aggr_Prop, M_Repo_Rate_6

# - Using domain expertise to refit model
# Add repo rate due to structural correlation with inflation rate, which may break down during crises, e.g., covid-19
# Add slc_acct_arr_dir_3_Change_Ind to align with TFD-model's input space
# Remove (PerfSpell_Maturity_Aggr_Mean, g0_Delinq_Num) due to statistical insignificance
# Add PerfSpell_Num for strategic purposes
vars <- c("g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", #"g0_Delinq_Num", 
          "Arrears", "PerfSpell_Num", "AgeToTerm_Aggr_Mean", 
          "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
          "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Final model
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(c(vars), collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
# AIC: 121610.7 Harrell's C: 0.9966     





# ------ 5 Stratum-based interaction effects

# ------ 5.1 Which specific interaction effects, chosen by greatest Harrell's c, are best?
# An interaction effect implies coefficients specific to a given stratum for a given variable

# - Initialize variables to be tested
vars <- c("factor(PerfSpell_Num_binned) * g0_Delinq_SD_4", "factor(PerfSpell_Num_binned) * Arrears", 
          "factor(PerfSpell_Num_binned) * BalanceToPrincipal","factor(PerfSpell_Num_binned) * pmnt_method_grp")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genObjPath)
### RESULTS: Best AIC-results: g0_Delinq_SD_4, pmnt_method_grp, Arrears, BalanceToPrincipal
# Best Harrell's C-statistics: g0_Delinq_SD_4, Arrears, BalanceToPrincipal, pmnt_method_grp
# Some p-values could not be reliably obtained, indicating model instability

# - Using domain expertise to refit model interactively
vars <- c("g0_Delinq_SD_4*factor(PerfSpell_Num_binned)", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave",
          "Arrears", "PerfSpell_Num", "AgeToTerm_Aggr_Mean", 
          "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
          "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Final model
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(c(vars), collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
# AIC: 113359.9 Harrell's C: 0.9985  
### RESULTS: Interacting the stratum with (g0_Delinq_SD_4, Arrears, pmnt_method_grp) 
# led to an unstable model with statistical insignificance in the affected coefficients
# However, some factor-levels had some significance for g0_Delinq_SD_4, which suggests that
# delinquency volatility may have more pronounced effects in subsequent spells; an idea for feature engineering

### CONCLUSION: The model could not be improved with interaction effects





# ------ 6. Final model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# - remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Initialize variables
vars2 <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Build model based on variables
cox_PWPST <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars2,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)
### RESULTS: AIC: 121610.7 Harrel's C: 0.9966  

c <- coefficients(cox_PWPST)
(c <- data.table(Variable=names(c),Coefficient=c))
#1:                     g0_Delinq_SD_4   6.69540818813
#2: slc_acct_roll_ever_24_imputed_mean   0.77131142478
#3:                      g0_Delinq_Ave  -4.21817975335
#4:                            Arrears   0.00002067274
#5:                      PerfSpell_Num   0.11336658180
#6:                AgeToTerm_Aggr_Mean  -7.20254802087
#7:                 BalanceToPrincipal   0.17036000871
#8:        pmnt_method_grpMISSING_DATA   1.21272161739
#9:     pmnt_method_grpSalary/Suspense   0.73089509413
#10:           pmnt_method_grpStatement   0.34801337044
#11:                   InterestRate_Nom   5.27075922146
#12:      slc_acct_arr_dir_3_Change_Ind   0.89997373162
#13:                     M_DTI_Growth_9  -2.07779970769
#14:               M_Inflation_Growth_6   4.48567772506
#15:                      M_Repo_Rate_6 -11.66007375132

# -Test Goodness of fit using bootstrapped B-statistics (1-KS statistic) over single-factor models
csTable_PWPST <- csTable(datCredit_train, vars2, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", seedVal=1, numIt=10,
                         fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="DefaultStatus1", genPath=genPath)
### RESULTS: Top single-factor models: PerfSpell_Num + M_Inflation_Growth_6  + Arrears + AgeToTerm_Aggr_Mean  
# Results do not vary much from each other, at least not meaningfully

aicTable_PWPST <- aicTable(datCredit_train, variables=vars2, fldSpellID="PerfSpell_Key",
                           TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genPath)
### RESULTS: Top single-factor models: g0_Delinq_SD_4 + slc_acct_roll_ever_24_imputed_mean + slc_acct_arr_dir_3_Change_Ind  + Arrears
# Where the first 3 results have AIC values significantly different from the rest.

# Test accuracy using Harrell's c-statistic over single-factor models
concTable_PWPST <- concTable(datCredit_train, datCredit_valid, variables=vars2, 
                             fldSpellID="PerfSpell_Key", TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genPath)
### RESULTS: Top single-factor models (>85%): g0_Delinq_SD_4 + Arrears + slc_acct_roll_ever_24_imputed_mean + slc_acct_arr_dir_3_Change_Ind

# - Combine results into a single object
Table_PWPST <- concTable_PWPST[,1:2] %>% left_join(aicTable_PWPST, by ="Variable") %>% 
  left_join(data.table(csTable_PWPST$Results), by="Variable")

# - Test Goodnes-of-fit using Cox-Snell, having measured distance between residual distribution and unit exponential using KS-statistic
GoF_CoxSnell_KS(cox_PWPST, datCredit_train, GraphInd=TRUE, legPos=c(0.6,0.4), panelTitle="Prentice-Williams-Peterson (PWP) model: Advanced",
                fileName = paste0(genFigPath, "KS_Test_CoxSnellResiduals_Exp_PWPST", ".png"), dpi=280) # 0.6167
### RESULTS: Goodness of fit for the model seems to be a bit low.

# Save objects
pack.ffdf(paste0(genObjPath,"PWPST_Univariate_Models"), Table_PWPST)
pack.ffdf(paste0(genPath,"PWPST_Cox_Model"), cox_PWPST)

# - Cleanup
rm(datCredit_train, datCredit_valid, cox_PWPST, cox_PWPST_empty, cox_PWPST_step,
   aicTable_PWPST, concTable_PWPST, csTable_PWPST, c); gc()
