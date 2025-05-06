
### AB [2025-05-05]: Deprecated since the Cox-model from LifetimePD-TermStructure-RecurrentEvents 
# codebase is taken given its much greater maturity.


# ====================================== Cox PH Model ======================================
# The training dataset is used to fit various Cox PH models with variables split according
# grouped into themes. The models are evaluated at each fit and the significant set of 
# variables are subsequently chosen. The last modelling theme uses all the significant 
# variables from the previous themes to obtain the final set of significant variables.
# The modelling exercise is done first for default and then the competing risks, i.e.
# early settlement, write-off, and paid-up.
# -----------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller
# -----------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 2d.Data_Fusion.R

# -- Inputs:
#   - datCredit_real | Prepared from scipt 2d.
#
# -- Outputs:
#   - dat_train | Training dataset which is sampled (per performance spell) from datCredit_real
#   - dat_valid | Validation dataset which is sampled (per performance spell) from datCredit_real
#   - Some graphs to explore the validility of the training- and validation datset
# -----------------------------------------------------------------------------------------


# ------ 1. Preliminaries
# --- Loaing in the required dataset
if (!exists('dat_train')) unpack.ffdf(paste0(genPath,"dat_train"), tempPath)
if (!exists('dat_valid')) unpack.ffdf(paste0(genPath,"dat_valid"), tempPath)

# --- Remove some unessary variables to improve performance
dat_train <- dat_train[,-(63:200)]
dat_valid <- dat_valid[,-(63:200)] # use dplyr to select columns

# --- Transforming the Default status into a numeric variable to ensure compatibility with time-dependent AUC/ ROC functions
dat_train[,DefaultStatus1 := as.numeric(DefaultStatus1)]; dat_train[,DefaultStatus1 := ifelse(DefaultStatus1==1,0,1)]
dat_valid[,DefaultStatus1 := as.numeric(DefaultStatus1)]; dat_valid[,DefaultStatus1 := ifelse(DefaultStatus1==1,0,1)]

# --- Grouping the data according to PerfSpell_Key and TimeInPerfSpell
dat_train <- dat_train %>% group_by(PerfSpell_Key, TimeInPerfSpell)
dat_valid <- dat_valid %>% group_by(PerfSpell_Key, TimeInPerfSpell)




# ------------ Default Models
# ------ 1. Theme: Account-level information
# --- Raw Variables
# - Fit a Cox model with the raw variables
cph_def_ALI_raw1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          Term + Receipt_Inf_wins + Balance_wins + Instalment_wins + Principal_wins +
                          InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind + FurtherLoan_Amt_wins +
                          Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins + BalanceToTerm_wins,
                        data=dat_train, id=PerfSpell_Key)
### RESULTS:~  Insignificant variables:
###            [FurtherLoan_Amt_wins]
###            Remove insignificant variables and refit model

cph_def_ALI_raw2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          Term + Receipt_Inf_wins + Balance_wins + Instalment_wins + Principal_wins +
                          InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                          Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins + BalanceToTerm_wins,
                        data=dat_train, id=PerfSpell_Key)
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_def_ALI_raw2 <- swiss_model(dat_valid, cph_def_ALI_raw2, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                   AUC_explain = FALSE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36))

### RESULTS:~ [Balance_wins] + [Principal_wins] + [BalanceToTerm_wins] have very high (>10) VIFs
###           Remove [Balance_wins] and [BalanceToTerm_wins] and refit model
start_time_def_ALI_raw <- Sys.time()
cph_def_ALI_raw3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          Term + Receipt_Inf_wins + Instalment_wins + Principal_wins +
                          InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                          Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins,
                        data=dat_train, id=PerfSpell_Key)
end_time_def_ALI_raw <- Sys.time()
(RunTime_def_ALI_Raw <- end_time_def_ALI_raw - start_time_def_ALI_raw)
### RESULTS:~  No insignificant variables.
###            Runtime = 8.0202742
###            Proceed to assessment.

assess_def_ALI_raw3 <- swiss_model(dat_valid, cph_def_ALI_raw3, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                   AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60))
### RESULTS:~ No unreasonably high VIFs
###           AUC-values: 12-months = 0.9222941
###                       24-months = 0.9140966
###                       36-months = 0.9043979
###                       60-months = 0.8796845

# --- Scaled Variables
# - Fit a Cox model with the scaled variables
cph_def_ALI_scaled1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             Term + Receipt_Inf_wins_mm_scaled + Balance_wins_norm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled +
                             InterestRate_Nom +InterestRate_Margin + FurtherLoan_Ind + FurtherLoan_Amt_wins_mm_scaled +
                             Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled + BalanceToTerm_wins_norm_scaled,
                            data=dat_train, id=PerfSpell_Key)
### RESULTS:~  Insignificant variables:
###            [FurtherLoan_Amt_wins_mm_scaled]
###            Remove [FurtherLoan_Amt_wins_mm_scaled] and refit model
cph_def_ALI_scaled2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             Term + Receipt_Inf_wins_mm_scaled + Balance_wins_norm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled +
                             InterestRate_Nom +InterestRate_Margin + FurtherLoan_Ind +
                             Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled + BalanceToTerm_wins_norm_scaled,
                          data=dat_train, id=PerfSpell_Key)
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_def_ALI_scaled1 <- swiss_model(dat_valid, cph_def_ALI_scaled2, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                      AUC_explain = FALSE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60))
### RESULTS:~ [Balance_wins_norm_scaled] + [Principal_wins_mm_scaled] + [BalanceToTerm_wins_norm_scaled] have very high (>10) VIFs
###           Remove [Balance_wins_norm_scaled] and [BalanceToTerm_wins_norm_scaled] and refit model

start_time_def_ALI_scaled <- Sys.time()
cph_def_ALI_scaled3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled +
                             InterestRate_Nom +InterestRate_Margin + FurtherLoan_Ind +
                             Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled,
                           data=dat_train, id=PerfSpell_Key)
end_time_def_ALI_scaled <- Sys.time()
(RunTime_def_ALI_Scaled <- end_time_def_ALI_scaled - start_time_def_ALI_scaled)
### RESULTS:~ No insignificant variables.
###           Runtime = 8.561359 mins.
###           Proceed to assessment.

assess_def_ALI_scaled2 <- swiss_model(dat_valid, cph_def_ALI_scaled3, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60))
### RESULTS:~ No unreasonably high VIFs
###           AUC-values: 12-months = 0.9222941
###                       24-months = 0.9140966
###                       36-months = 0.9043979
###                       60-months = 0.8796845

# --- Raw vs Scaled Variables

assess_def_ALI_raw2$AUC - assess_def_ALI_scaled2$AUC # AUC differences

RunTime_def_ALI_Raw - RunTime_def_ALI_Scaled # Runtime difference

summary(cph_def_ALI_raw3)$coefficients[,3] - summary(cph_def_ALI_scaled3)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference in AUCs, doesn't provide clear cut guidance
###           Raw variables result in superior convergence time
###           Raw variables result in smaller standard errors of the coefficients
###           Use raw variables



# ------ 2. Theme: Macroeconomic variables
# - Fit a Cox model with the variables
cph_def_mac1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                      M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + M_Emp_Growth +
                      M_RealGDP_Growth + M_RealIncome_Growth,
                     data=dat_train, id=PerfSpell_Key)
### RESULTS:~  Insignificant variables:
###            [M_RealGDP_Growth]
###            Remove insignificant variables and refit model

cph_def_mac2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                      M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + M_Emp_Growth +
                      M_RealIncome_Growth,
                    data=dat_train, id=PerfSpell_Key)
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_def_mac1 <- swiss_model(dat_valid, cph_def_mac1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                               AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60))
### RESULTS:~ [M_Emp_Growth] + [M_Real_GDP_Growth] + [M_RealIncome_Growth] have high VIFs; ignore since we know that they are highly correlated, but necessary to add in model.
###           AUC-values: 12-months = 0.5854601
###                       24-months = 0.5911637
###                       36-months = 0.5900309
###                       60-months = 0.5732540



# ------ 3. Theme: Delinquency/ Performance Data
# - Fit a Cox model with the appropriate variables
cph_def_del1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                      PerfSpell_Num + slc_acct_roll_ever_24_imputed + #g0_Delinq + slc_acct_arr_dir_3
                    ,data=dat_train, id=PerfSpell_Key)
### NOTE:~ Tried to fit g0_Delinq, but model convergence is impossible.
###        Binned all [g0_Delinq]=3 into [g0_Delinq]=2, from which model convergence was  achieved. However, the proportional hazard was very large and the variable was very insignificant...?
###        [slc_acct_dir_3] also results in no model convergence and was subsequnetly removed.
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_def_del1 <- swiss_model(dat_valid, cph_def_del1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                               AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60))
### RESULTS:~ No unreasonably high VIFs
###           AUC-values: 12-months = 0.9235616
###                       24-months = 0.8448925
###                       36-months = 0.8492486
###                       60-months = 0.831844


# ------ 4. Behavioral Variables
# --- Raw variables (without missing value indicators)
# - Fit a Cox model with the raw variables
start_time_def_beh_raw <- Sys.time()
cph_def_beh_raw1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          slc_pmnt_method + slc_acct_pre_lim_perc_imputed +
                          slc_past_due_amt_imputed_wins + slc_acct_prepaid_perc_dir_12_imputed_wins,
                        data=dat_train, id=PerfSpell_Key)
end_time_def_beh_raw <- Sys.time()
(RunTime_def_beh_Raw <- end_time_def_beh_raw - start_time_def_beh_raw)
### RESULTS:~  No insignificant variables.
###            Runtime = 7.169934 mins.
###            Proceed to assessment.

(assess_def_beh_raw <- swiss_model(dat_valid, cph_def_beh_raw1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                   AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# --- Scaled Variables (without missing value indicators)
# - Fit a Cox model with the scaled variables
start_time_def_beh_scaled <- Sys.time()
cph_def_beh_scaled1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             slc_pmnt_method + slc_acct_pre_lim_perc_imputed +
                             slc_past_due_amt_imputed_wins_mm_scaled + slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled,
                           data=dat_train, id=PerfSpell_Key)
end_time_def_beh_scaled <- Sys.time()
(RunTime_def_beh_Scaled <- end_time_def_beh_scaled - start_time_def_beh_scaled)
### RESULTS:~  No insignificant variables.
###            Runtime = 6.610058 mins.
###            Proceed to assessment.

(assess_def_beh_scaled <- swiss_model(dat_valid, cph_def_beh_scaled1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# --- Raw Variables (with missing value indicators)
# - Fit a Cox model with the raw variables (with value indicators)
start_time_def_beh_raw_val <- Sys.time()
cph_def_beh_raw_val1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc +
                              slc_past_due_amt_imputed_wins:value_ind_slc_past_due_amt + slc_acct_prepaid_perc_dir_12_imputed_wins:value_ind_slc_acct_prepaid_perc_dir_12
                            , data=dat_train, id=PerfSpell_Key)
end_time_def_beh_raw_val <- Sys.time()
(RunTime_def_beh_Raw_val <- end_time_def_beh_raw_val - start_time_def_beh_raw_val)
### RESULTS:~  No insignificant variables.
###            Runtime = 6.222053 mins.
###            Proceed to assessment.

(assess_def_beh_raw_val <- swiss_model(dat_valid, cph_def_beh_raw_val1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                       AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# - Fit a Cox model with the raw variables (with value indicators)
start_time_def_beh_scaled_val <- Sys.time()
cph_def_beh_scaled_val1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                 slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc +
                                 slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                              , data=dat_train, id=PerfSpell_Key)
end_time_def_beh_scaled_val <- Sys.time()
(RunTime_def_beh_Scaled_val <- end_time_def_beh_scaled_val - start_time_def_beh_scaled_val)
### RESULTS:~  No insignificant variables.
###            Runtime = 6.222053 mins.
###            Remove insignificant variables and refit model

(assess_def_beh_scaled_val <- swiss_model(dat_valid, cph_def_beh_scaled_val1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                          AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# --- Raw vs scaled (without missing value indicators)
assess_def_beh_raw$AUC - assess_def_beh_scaled$AUC # AUC differences

RunTime_def_beh_Raw - RunTime_def_beh_Scaled # Runtime difference

summary(cph_def_beh_raw1)$coefficients[,3] - summary(cph_def_beh_scaled1)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference in AUCs, doesn't provide clear cut guidance
###           Scaled variables result in superior convergence time (~30 seconds)
###           Scaled variables result in smaller standard errors of the coefficients
###           Compare against missing value indicators

# --- Raw (with missing value indicators) vs scaled (with missing value indicators)
assess_def_beh_raw_val$AUC - assess_def_beh_scaled_val$AUC # AUC differences

RunTime_def_beh_Raw_val - RunTime_def_beh_Scaled_val # Runtime difference

summary(cph_def_beh_raw_val1)$coefficients[,3] - summary(cph_def_beh_scaled_val1)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference, doesn't provide clear cut guidance
###           Scaled variables result in superior convergence time (~15 seconds)
###           Scaled variables result in smaller standard errors of the coefficients
###           Compare against missing value indicators

# --- Raw (without missing value indicators) vs scaled (with missing value indicators)
assess_def_beh_raw$AUC - assess_def_beh_scaled_val$AUC # AUC differences

RunTime_def_beh_Raw - RunTime_def_beh_Scaled_val # Runtime difference

summary(cph_def_beh_raw_val1)$coefficients[,3] - summary(cph_def_beh_scaled_val1)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference, doesn't provide clear cut guidance
###           Scaled variables (with missing value indicators) result in superior convergence time (~60 seconds)
###           Scaled variables (with missing value indicators) result in smaller standard errors of the coefficients
###           Use scaled variables (with missing value indicators)




# ------ 5. Combined Variables
# --- Raw variables
# - Fit a Cox model with all chosen variables
start_time_def_com_raw1 <- Sys.time()
cph_def_com_raw1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               Term + Receipt_Inf_wins + Instalment_wins + Principal_wins + # ALI
                               InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                               Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins +
                               M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                               M_Emp_Growth + M_RealIncome_Growth +
                               PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                               slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                               slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                               slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                          , data=dat_train, id=PerfSpell_Key)
end_time_def_com_raw1 <- Sys.time()
(RunTime_def_com_raw1 <- end_time_def_com_raw1 - start_time_def_com_raw1)
### RESULTS:~ Convergence achieved, but coefficients/ hazard ratios may be inaccurate
###           Insignificant variables:
###           [Term] + [Receipt_Inf_wins]+ [InterestRate_Nom] + [FurtherLoan_Ind] + [Redrawn_Amt_wins] + [LN_TPE]
###           [M_Inflation_Growth] + [M_DTI_Growth] + [M_Emp_Growth] + [M_RealIncome_Growth] +
###           [slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove [InterestRate_Margin] and refit

start_time_def_com_raw2 <- Sys.time()
cph_def_com_raw2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          Term + Receipt_Inf_wins + Instalment_wins + Principal_wins + # ALI
                          InterestRate_Nom + FurtherLoan_Ind +
                          Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins +
                          M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                          M_Emp_Growth + M_RealIncome_Growth +
                          PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                          slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                          slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                          slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                      , data=dat_train, id=PerfSpell_Key)
end_time_def_com_raw2 <- Sys.time()
(RunTime_def_com_raw2 <- end_time_def_com_raw2 - start_time_def_com_raw2)
### RESULTS:~ Convergence successful.
###           Insignificant variables:
###           [InterestRate_Nom] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove insignificant variables and refit model

start_time_def_com_raw3 <- Sys.time()
cph_def_com_raw3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          Term + Receipt_Inf_wins + Instalment_wins + Principal_wins + # ALI
                          FurtherLoan_Ind +
                          Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins +
                          M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                          M_Emp_Growth + M_RealIncome_Growth +
                          PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                          slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                          slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt
                      , data=dat_train, id=PerfSpell_Key)
end_time_def_com_raw3 <- Sys.time()
(RunTime_def_com_raw3 <- end_time_def_com_raw3 - start_time_def_com_raw3)
### RESULTS:~ Convergence successful.
###           No insignificant variables.
###           Proceed to assessment.

(assess_def_com_raw <- swiss_model(dat_valid, cph_def_com_raw3, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                   AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9859523
# AUC-24 months = 0.9876597
# AUC-36 months = 0.9883425
# AUC-60 months = 0.9845827



# --- Scaled variables
start_time_def_com_scaled1 <- Sys.time()
cph_def_com_scaled1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled + # ALI
                             InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                             Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled +
                             M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                             M_Emp_Growth + M_RealIncome_Growth +
                             PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                             slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                             slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                             slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                          , data=dat_train, id=PerfSpell_Key)
end_time_def_com_scaled1 <- Sys.time()
(RunTime_def_com_scaled1 <- end_time_def_com_scaled1 - start_time_def_com_scaled1)
### RESULTS:~ Convergence achieved, but coefficients/ hazard ratios may be inaccurate
###           Insignificant variables:
###           [Term] + [Receipt_Inf_wins]+ [InterestRate_Nom] + [FurtherLoan_Ind] + [Redrawn_Amt_wins] + [LN_TPE]
###           [M_Inflation_Growth] + [M_DTI_Growth] + [M_Emp_Growth] + [M_RealIncome_Growth] +
###           [slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove [InterestRate_Margin] and refit

start_time_def_com_scaled2 <- Sys.time()
cph_def_com_scaled2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled + # ALI
                             InterestRate_Nom + FurtherLoan_Ind +
                             Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled +
                             M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                             M_Emp_Growth + M_RealIncome_Growth +
                             PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                             slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                             slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                             slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                          , data=dat_train, id=PerfSpell_Key)
end_time_def_com_scaled2 <- Sys.time()
(RunTime_def_com_scaled2 <- end_time_def_com_scaled2 - start_time_def_com_scaled2)
### RESULTS:~ Convergence achieved, but coefficients/ hazard ratios may be inaccurate
###           Insignificant variables:
###           [InterestRate_Nom] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove insignificant variables and refit model

start_time_def_com_scaled3 <- Sys.time()
cph_def_com_scaled3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled + # ALI
                             FurtherLoan_Ind +
                             Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled +
                             M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                             M_Emp_Growth + M_RealIncome_Growth +
                             PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                             slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                             slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt
                          , data=dat_train, id=PerfSpell_Key)
end_time_def_com_scaled3 <- Sys.time()
(RunTime_def_com_scaled3 <- end_time_def_com_scaled3 - start_time_def_com_scaled3)
### RESULTS:~ Convergence successful.
###           No insignificant variables.
###           Proceed to assessment.

(assess_def_com_scaled <- swiss_model(dat_valid, cph_def_com_raw3, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9859523
# AUC-24 months = 0.9876597
# AUC-36 months = 0.9883425
# AUC-60 months = 0.9845827


# --- Raw -vs scaled variables
assess_def_com_raw$AUC - assess_def_com_scaled$AUC # AUC differences

RunTime_def_com_raw3 - RunTime_def_com_scaled3 # Runtime difference

summary(cph_def_com_raw3)$coefficients[,3] - summary(cph_def_com_scaled3)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Nol difference, doesn't provide clear cut guidance.
###           Scaled variables result in superior convergence time (~5 seconds)
###           Scaled variables result in smaller standard errors of the coefficients
###           Use scaled variables


(assess_def_com_scaled <- swiss_model(dat_valid, cph_def_com_scaled3, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "DefaultStatus1"), predict.time = c(12,24,36,60)))
### RESULTS:~ PH assumption may be violated for some variables.
###           Macroeconomic variables have a bit high VIFs (expected); [Principal_wins_mm_scaled] also has a bit high VIF.
###           

# PH_Test <-cox.zph(cph_def_com_scaled3, transform='log')
# ggcoxzph(PH_Test)

# Model_ANOVA <- anova(cph_def_com_scaled3, test = 'chisq')

# --- Save final model to disk
pack.ffdf(paste0(genPath,"cph_def_fin"), cph_def_com_scaled3)






# ------------ Early Settlement Models
# ------ 1. Theme: Account-level information
# --- Raw Variables
# - Fit a Cox model with the raw variables
cph_set_ALI_raw1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            Term + Receipt_Inf_wins + Balance_wins + Instalment_wins + Principal_wins +
                            InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind + FurtherLoan_Amt_wins +
                            Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins + BalanceToTerm_wins,
                          data=dat_train, id=PerfSpell_Key)
### RESULTS:~  Insignificant variables:
###            [FurtherLoan_Amt_wins]
###            Remove insignificant variables and refit model

cph_set_ALI_raw2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            Term + Receipt_Inf_wins + Balance_wins + Instalment_wins + Principal_wins +
                            InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                            Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins + BalanceToTerm_wins,
                          data=dat_train, id=PerfSpell_Key)
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_set_ALI_raw1 <- swiss_model(dat_valid, cph_set_ALI_raw2, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                   AUC_explain = FALSE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36))

### RESULTS:~ [Balance_wins] + [Principal_wins] + [BalanceToTerm_wins] have very high (>10) VIFs
###           Remove [Balance_wins] and [BalanceToTerm_wins] and refit model
start_time_set_ALI_raw <- Sys.time()
cph_set_ALI_raw3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            Term + Receipt_Inf_wins + Instalment_wins + Principal_wins +
                            InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                            Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins,
                          data=dat_train, id=PerfSpell_Key)
end_time_set_ALI_raw <- Sys.time()
(RunTime_set_ALI_Raw <- end_time_set_ALI_raw - start_time_set_ALI_raw)
### RESULTS:~  No insignificant variables.
###            Runtime = 8.0202742
###            Proceed to assessment.

assess_set_ALI_raw2 <- swiss_model(dat_valid, cph_set_ALI_raw3, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                   AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60))
### RESULTS:~ No unreasonably high VIFs
###           AUC-values: 12-months = 0.9222941
###                       24-months = 0.9140966
###                       36-months = 0.9043979
###                       60-months = 0.8796845

# --- Scaled Variables
# - Fit a Cox model with the scaled variables
cph_set_ALI_scaled1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               Term + Receipt_Inf_wins_mm_scaled + Balance_wins_norm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled +
                               InterestRate_Nom +InterestRate_Margin + FurtherLoan_Ind + FurtherLoan_Amt_wins_mm_scaled +
                               Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled + BalanceToTerm_wins_norm_scaled,
                             data=dat_train, id=PerfSpell_Key)
### RESULTS:~  Insignificant variables:
###            [FurtherLoan_Amt_wins_mm_scaled]
###            Remove [FurtherLoan_Amt_wins_mm_scaled] and refit model
cph_set_ALI_scaled2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               Term + Receipt_Inf_wins_mm_scaled + Balance_wins_norm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled +
                               InterestRate_Nom +InterestRate_Margin + FurtherLoan_Ind +
                               Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled + BalanceToTerm_wins_norm_scaled,
                             data=dat_train, id=PerfSpell_Key)
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_set_ALI_scaled1 <- swiss_model(dat_valid, cph_set_ALI_scaled2, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                      AUC_explain = FALSE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60))
### RESULTS:~ [Balance_wins_norm_scaled] + [Principal_wins_mm_scaled] + [BalanceToTerm_wins_norm_scaled] have very high (>10) VIFs
###           Remove [Balance_wins_norm_scaled] and [BalanceToTerm_wins_norm_scaled] and refit model

start_time_set_ALI_scaled <- Sys.time()
cph_set_ALI_scaled3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled +
                               InterestRate_Nom +InterestRate_Margin + FurtherLoan_Ind +
                               Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled,
                             data=dat_train, id=PerfSpell_Key)
end_time_set_ALI_scaled <- Sys.time()
(RunTime_set_ALI_Scaled <- end_time_set_ALI_scaled - start_time_set_ALI_scaled)
### RESULTS:~ No insignificant variables.
###           Runtime = 8.561359 mins.
###           Proceed to assessment.

assess_set_ALI_scaled2 <- swiss_model(dat_valid, cph_set_ALI_scaled3, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60))
### RESULTS:~ No unreasonably high VIFs
###           AUC-values: 12-months = 0.9222941
###                       24-months = 0.9140966
###                       36-months = 0.9043979
###                       60-months = 0.8796845

# --- Raw vs Scaled Variables

assess_set_ALI_raw2$AUC - assess_set_ALI_scaled2$AUC # AUC differences

RunTime_set_ALI_Raw - RunTime_set_ALI_Scaled # Runtime difference

summary(cph_set_ALI_raw3)$coefficients[,3] - summary(cph_set_ALI_scaled3)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference in AUCs, doesn't provide clear cut guidance
###           Raw variables result in superior convergence time
###           Raw variables result in smaller standard errors of the coefficients
###           Use raw variables



# ------ 2. Theme: Macroeconomic variables
# - Fit a Cox model with the variables
cph_set_mac1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                        M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + M_Emp_Growth +
                        M_RealGDP_Growth + M_RealIncome_Growth,
                      data=dat_train, id=PerfSpell_Key)
### RESULTS:~  Insignificant variables:
###            [M_RealGDP_Growth]
###            Remove insignificant variables and refit model

cph_set_mac2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                        M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + M_Emp_Growth +
                        M_RealIncome_Growth,
                      data=dat_train, id=PerfSpell_Key)
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_set_mac1 <- swiss_model(dat_valid, cph_set_mac1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                               AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60))
### RESULTS:~ [M_Emp_Growth] + [M_Real_GDP_Growth] + [M_RealIncome_Growth] have high VIFs; ignore since we know that they are highly correlated, but necessary to add in model.
###           AUC-values: 12-months = 0.5854601
###                       24-months = 0.5911637
###                       36-months = 0.5900309
###                       60-months = 0.5732540



# ------ 3. Theme: Delinquency/ Performance Data
# - Fit a Cox model with the variables
cph_set_del1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                        PerfSpell_Num + slc_acct_roll_ever_24_imputed # + g0_Delinq + slc_acct_arr_dir_3
                      ,data=dat_train, id=PerfSpell_Key)
### NOTE:~ Tried to fit g0_Delinq, but model convergence is impossible.
###        Binned all [g0_Delinq]=3 into [g0_Delinq]=2, from which model convergence was  achieved. However, the proportional hazard was very large and the variable was very insignificant...?
###        [slc_acct_dir_3] also results in no model convergence and was subsequnetly removed.
### RESULTS:~  No insignificant variables.
###            Proceed to assessment.

assess_set_del1 <- swiss_model(dat_valid, cph_set_del1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                               AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60))
### RESULTS:~ No unreasonably high VIFs
###           AUC-values: 12-months = 0.9235616
###                       24-months = 0.8448925
###                       36-months = 0.8492486
###                       60-months = 0.831844


# ------ 4. Behavioral Variables
# --- Raw variables (without missing value indicators)
# - Fit a Cox model with the raw variables
start_time_set_beh_raw <- Sys.time()
cph_set_beh_raw1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            slc_pmnt_method + slc_acct_pre_lim_perc_imputed +
                            slc_past_due_amt_imputed_wins + slc_acct_prepaid_perc_dir_12_imputed_wins,
                          data=dat_train, id=PerfSpell_Key)
end_time_set_beh_raw <- Sys.time()
(RunTime_set_beh_Raw <- end_time_set_beh_raw - start_time_set_beh_raw)
### RESULTS:~  No insignificant variables.
###            Runtime = 7.169934 mins.
###            Proceed to assessment.

(assess_set_beh_raw <- swiss_model(dat_valid, cph_set_beh_raw1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                   AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# --- Scaled Variables (without missing value indicators)
# - Fit a Cox model with the scaled variables
start_time_set_beh_scaled <- Sys.time()
cph_set_beh_scaled1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               slc_pmnt_method + slc_acct_pre_lim_perc_imputed +
                               slc_past_due_amt_imputed_wins_mm_scaled + slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled,
                             data=dat_train, id=PerfSpell_Key)
end_time_set_beh_scaled <- Sys.time()
(RunTime_set_beh_Scaled <- end_time_set_beh_scaled - start_time_set_beh_scaled)
### RESULTS:~  No insignificant variables.
###            Runtime = 6.610058 mins.
###            Proceed to assessment.

(assess_set_beh_scaled <- swiss_model(dat_valid, cph_set_beh_scaled1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# --- Raw Variables (with missing value indicators)
# - Fit a Cox model with the raw variables (with value indicators)
start_time_set_beh_raw_val <- Sys.time()
cph_set_beh_raw_val1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                                slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc +
                                slc_past_due_amt_imputed_wins:value_ind_slc_past_due_amt + slc_acct_prepaid_perc_dir_12_imputed_wins:value_ind_slc_acct_prepaid_perc_dir_12
                              , data=dat_train, id=PerfSpell_Key)
end_time_set_beh_raw_val <- Sys.time()
(RunTime_set_beh_Raw_val <- end_time_set_beh_raw_val - start_time_set_beh_raw_val)
### RESULTS:~  No insignificant variables.
###            Runtime = 6.222053 mins.
###            Proceed to assessment.

(assess_set_beh_raw_val <- swiss_model(dat_valid, cph_set_beh_raw_val1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                       AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# - Fit a Cox model with the raw variables (with value indicators)
start_time_set_beh_scaled_val <- Sys.time()
cph_set_beh_scaled_val1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                                   slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc +
                                   slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                                 , data=dat_train, id=PerfSpell_Key)
end_time_set_beh_scaled_val <- Sys.time()
(RunTime_set_beh_Scaled_val <- end_time_set_beh_scaled_val - start_time_set_beh_scaled_val)
### RESULTS:~  No insignificant variables.
###            Runtime = 6.222053 mins.
###            Remove insignificant variables and refit model

(assess_set_beh_scaled_val <- swiss_model(dat_valid, cph_set_beh_scaled_val1, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                          AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9466732
# AUC-24 months = 0.9574479
# AUC-36 months = 0.9624574
# AUC-60 months = 0.9607957

# --- Raw vs scaled (without missing value indicators)
assess_set_beh_raw$AUC - assess_set_beh_scaled$AUC # AUC differences

RunTime_set_beh_Raw - RunTime_set_beh_Scaled # Runtime difference

summary(cph_set_beh_raw1)$coefficients[,3] - summary(cph_set_beh_scaled1)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference in AUCs, doesn't provide clear cut guidance
###           Scaled variables result in superior convergence time (~30 seconds)
###           Scaled variables result in smaller standard errors of the coefficients
###           Compare against missing value indicators

# --- Raw (with missing value indicators) vs scaled (with missing value indicators)
assess_set_beh_raw_val$AUC - assess_set_beh_scaled_val$AUC # AUC differences

RunTime_set_beh_Raw_val - RunTime_set_beh_Scaled_val # Runtime difference

summary(cph_set_beh_raw_val1)$coefficients[,3] - summary(cph_set_beh_scaled_val1)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference, doesn't provide clear cut guidance
###           Scaled variables result in superior convergence time (~15 seconds)
###           Scaled variables result in smaller standard errors of the coefficients
###           Compare against missing value indicators

# --- Raw (without missing value indicators) vs scaled (with missing value indicators)
assess_set_beh_raw$AUC - assess_set_beh_scaled_val$AUC # AUC differences

RunTime_set_beh_Raw - RunTime_set_beh_Scaled_val # Runtime difference

summary(cph_set_beh_raw_val1)$coefficients[,3] - summary(cph_set_beh_scaled_val1)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Marginal difference, doesn't provide clear cut guidance
###           Scaled variables (with missing value indicators) result in superior convergence time (~60 seconds)
###           Scaled variables (with missing value indicators) result in smaller standard errors of the coefficients
###           Use scaled variables (with missing value indicators)




# ------ 5. Combined Variables
# --- Raw variables
# - Fit a Cox model with all chosen variables
start_time_set_com_raw1 <- Sys.time()
cph_set_com_raw1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            Term + Receipt_Inf_wins + Instalment_wins + Principal_wins + # ALI
                            InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                            Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins +
                            M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                            M_Emp_Growth + M_RealIncome_Growth +
                            PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                            slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                            slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                            slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                          , data=dat_train, id=PerfSpell_Key)
end_time_set_com_raw1 <- Sys.time()
(RunTime_set_com_raw1 <- end_time_set_com_raw1 - start_time_set_com_raw1)
### RESULTS:~ Convergence achieved, but coefficients/ hazard ratios may be inaccurate
###           Insignificant variables:
###           [Term] + [Receipt_Inf_wins]+ [InterestRate_Nom] + [FurtherLoan_Ind] + [Redrawn_Amt_wins] + [LN_TPE]
###           [M_Inflation_Growth] + [M_DTI_Growth] + [M_Emp_Growth] + [M_RealIncome_Growth] +
###           [slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove [InterestRate_Margin] and refit

start_time_set_com_raw2 <- Sys.time()
cph_set_com_raw2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            Term + Receipt_Inf_wins + Instalment_wins + Principal_wins + # ALI
                            InterestRate_Nom + FurtherLoan_Ind +
                            Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins +
                            M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                            M_Emp_Growth + M_RealIncome_Growth +
                            PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                            slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                            slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                            slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                          , data=dat_train, id=PerfSpell_Key)
end_time_set_com_raw2 <- Sys.time()
(RunTime_set_com_raw2 <- end_time_set_com_raw2 - start_time_set_com_raw2)
### RESULTS:~ Convergence successful.
###           Insignificant variables:
###           [InterestRate_Nom] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove insignificant variables and refit model

start_time_set_com_raw3 <- Sys.time()
cph_set_com_raw3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                            Term + Receipt_Inf_wins + Instalment_wins + Principal_wins + # ALI
                            FurtherLoan_Ind +
                            Redrawn_Amt_wins + Redraw_Ind + LN_TPE + AgeToTerm_wins +
                            M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                            M_Emp_Growth + M_RealIncome_Growth +
                            PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                            slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                            slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt
                          , data=dat_train, id=PerfSpell_Key)
end_time_set_com_raw3 <- Sys.time()
(RunTime_set_com_raw3 <- end_time_set_com_raw3 - start_time_set_com_raw3)
### RESULTS:~ Convergence successful.
###           No insignificant variables.
###           Proceed to assessment.

(assess_set_com_raw <- swiss_model(dat_valid, cph_set_com_raw3, PH_assumption = FALSE, VIF = FALSE, anova_explain = FALSE,
                                   AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9859523
# AUC-24 months = 0.9876597
# AUC-36 months = 0.9883425
# AUC-60 months = 0.9845827



# --- Scaled variables
start_time_set_com_scaled1 <- Sys.time()
cph_set_com_scaled1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled + # ALI
                               InterestRate_Nom + InterestRate_Margin + FurtherLoan_Ind +
                               Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled +
                               M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                               M_Emp_Growth + M_RealIncome_Growth +
                               PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                               slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                               slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                               slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                             , data=dat_train, id=PerfSpell_Key)
end_time_set_com_scaled1 <- Sys.time()
(RunTime_set_com_scaled1 <- end_time_set_com_scaled1 - start_time_set_com_scaled1)
### RESULTS:~ Convergence achieved, but coefficients/ hazard ratios may be inaccurate
###           Insignificant variables:
###           [Term] + [Receipt_Inf_wins]+ [InterestRate_Nom] + [FurtherLoan_Ind] + [Redrawn_Amt_wins] + [LN_TPE]
###           [M_Inflation_Growth] + [M_DTI_Growth] + [M_Emp_Growth] + [M_RealIncome_Growth] +
###           [slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove [InterestRate_Margin] and refit

start_time_set_com_scaled2 <- Sys.time()
cph_set_com_scaled2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled + # ALI
                               InterestRate_Nom + FurtherLoan_Ind +
                               Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled +
                               M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                               M_Emp_Growth + M_RealIncome_Growth +
                               PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                               slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                               slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt + 
                               slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12
                             , data=dat_train, id=PerfSpell_Key)
end_time_set_com_scaled2 <- Sys.time()
(RunTime_set_com_scaled2 <- end_time_set_com_scaled2 - start_time_set_com_scaled2)
### RESULTS:~ Convergence achieved, but coefficients/ hazard ratios may be inaccurate
###           Insignificant variables:
###           [InterestRate_Nom] + [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled:value_ind_slc_acct_prepaid_perc_dir_12]
###           Remove insignificant variables and refit model

start_time_set_com_scaled3 <- Sys.time()
cph_set_com_scaled3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=EarlySettle_Ind==1)~
                               Term + Receipt_Inf_wins_mm_scaled + Instalment_wins_mm_scaled + Principal_wins_mm_scaled + # ALI
                               FurtherLoan_Ind +
                               Redrawn_Amt_wins_mm_scaled + Redraw_Ind + LN_TPE + AgeToTerm_wins_mm_scaled +
                               M_Repo_Rate + M_Inflation_Growth + M_Inflation_Growth + M_DTI_Growth + # Macroeconomic
                               M_Emp_Growth + M_RealIncome_Growth +
                               PerfSpell_Num + slc_acct_roll_ever_24_imputed + # Delinquency/ performance
                               slc_pmnt_method + slc_acct_pre_lim_perc_imputed:value_ind_slc_acct_pre_lim_perc + # Behavioral
                               slc_past_due_amt_imputed_wins_mm_scaled:value_ind_slc_past_due_amt
                             , data=dat_train, id=PerfSpell_Key)
end_time_set_com_scaled3 <- Sys.time()
(RunTime_set_com_scaled3 <- end_time_set_com_scaled3 - start_time_set_com_scaled3)
### RESULTS:~ Convergence successful.
###           No insignificant variables.
###           Proceed to assessment.

(assess_set_com_scaled <- swiss_model(dat_valid, cph_set_com_raw3, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
# AUC-12 months = 0.9859523
# AUC-24 months = 0.9876597
# AUC-36 months = 0.9883425
# AUC-60 months = 0.9845827


# --- Raw -vs scaled variables
assess_set_com_raw$AUC - assess_set_com_scaled$AUC # AUC differences

RunTime_set_com_raw3 - RunTime_set_com_scaled3 # Runtime difference

summary(cph_set_com_raw3)$coefficients[,3] - summary(cph_set_com_scaled3)$coefficients[,3] # Difference in coefficients' standard errors (since both models have the same covariates)
### RESULTS:~ Nol difference, doesn't provide clear cut guidance.
###           Scaled variables result in superior convergence time (~5 seconds)
###           Scaled variables result in smaller standard errors of the coefficients
###           Use scaled variables


(assess_set_com_scaled <- swiss_model(dat_valid, cph_set_com_scaled3, PH_assumption = FALSE, VIF = TRUE, anova_explain = FALSE,
                                      AUC_explain = TRUE, AUC_input.v = c("TimeInPerfSpell", "EarlySettle_Ind"), predict.time = c(12,24,36,60)))
### RESULTS:~ PH assumption may be violated for some variables.
###           Macroeconomic variables have a bit high VIFs (expected); [Principal_wins_mm_scaled] also has a bit high VIF.
###           

# PH_Test <-cox.zph(cph_set_com_scaled3, transform='log')
# ggcoxzph(PH_Test)

# Model_ANOVA <- anova(cph_set_com_scaled3, test = 'chisq')

# --- Save final model to disk
pack.ffdf(paste0(genPath,"cph_set_fin"), cph_set_com_scaled3)













# ------ Testing why g0_Delinq results in no model convergence and yields a very high coefficient and proportional hazard value
unique(dat_valid$g0_Delinq)
# 0, 1, 2, 3
# 
dat_valid[DefaultStatus1==1, .N]
# 37 462
dat_valid[g0_Delinq==3, .N]
# 37 462
dat_valid[g0_Delinq==3 & DefaultStatus1==1, .N]
# 37 462
all.equal(dat_valid[g0_Delinq==3,PerfSpell_Key], dat_valid[DefaultStatus1==1,PerfSpell_Key])
# 37 462
### RESUTLS:~ All records where g0_Delinq = 3 indicates DefaultStatus1 = 1
