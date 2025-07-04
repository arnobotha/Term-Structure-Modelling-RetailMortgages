# ======================================= INPUT SPACE: Cox PO (PWP) ADVANCED============================
# Divide data into thematic groups and perform data analysis on them to compile an input space for 
# and advanced discrete-time hazard model, using the PWPST-definition for recurrency.
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





# ------ 1. Preliminaries

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]
# remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Fit an "empty" model as a performance gain, used within some diagnostic functions
modLR_base <- glm(PerfSpell_Event ~ 1, data=datCredit_train, family="binomial")

# - Fit a baseline model for stepwise forward selection
modLR_base2 <- glm(PerfSpell_Event ~ log(TimeInPerfSpell)*PerfSpell_Num_binned, 
                   data=datCredit_train, family="binomial")
# Insight interactively mined from modelling theme 2a-b was used in fitting this model.



# ------ 2a. Modelling theme: Handling time towards embedding baseline hazard h_0(t)
# Which way in embedding time is best?

# --- Single-factor models: Single time variable
modLR <- glm(PerfSpell_Event ~  TimeInPerfSpell,
                   data = datCredit_train, family="binomial")
summary(modLR)
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  177,477; McFadden R^2:  0.26%; AUC:  59.12%

# --- Single-factor models: Single time variable (transform)
modLR <- glm(PerfSpell_Event ~  log(TimeInPerfSpell),
             data = datCredit_train, family="binomial")
summary(modLR)
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  176,610; McFadden R^2:  0.75%; AUC:  59.12%
# Log-transformation affords a better goodness-of-fit, though AUC understandably remains unchanged

# --- Single-factor models: Single time variable (binned)
modLR <- glm(PerfSpell_Event ~  Time_Binned,
             data = datCredit_train, family="binomial")
summary(modLR)
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  175,173; McFadden R^2:  1.58%; AUC:  62.83%.
# Binned version of time substantially outperforms previous iteration.
# Binning was interactively tweaked such that most resulting bins are statistically significant

# --- Single-factor models: Single time variable (splines)
# NOTE: Splines set to the same number of bins in [Time_Binned] for approximate comparison
modLR <- glm(PerfSpell_Event ~  ns(TimeInPerfSpell, df=20),
             data = datCredit_train, family="binomial")
summary(modLR)
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   175,033; McFadden R^2:  1.66%; AUC:  63.25%.
# A marginal improvement upon the previous iteration, though still significant.

### CONCLUSION: Binning and or natural splines allows for incorporating nonlinear effects,
# which is to be expected for the baseline hazard function



# ------ 2b. Modelling theme: Handling time towards embedding baseline hazard h_0(t) across spell number
# Does a stratified model work better, i.e., a separate baseline per spell number

# --- Single-factor models: Single time variable (transform)
modLR <- glm(PerfSpell_Event ~  log(TimeInPerfSpell)*PerfSpell_Num_binned,
             data = datCredit_train, family="binomial")
summary(modLR)
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   171,796; McFadden R^2:  3.46%; AUC:  62.50%.
# Yes, an interaction effect between the time in spell and spell number is not only statistically significant,
# but also improves the model across all metrics.

# --- Single-factor models: Single time variable (binned) with interaction effect with spell number
modLR <- glm(PerfSpell_Event ~  Time_Binned + log(TimeInPerfSpell):PerfSpell_Num_binned,
             data = datCredit_train, family="binomial")
summary(modLR)
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   171,715; McFadden R^2:  3.52%; AUC:  66.94%
# Even better performance using the binned version of time, as expected




# ------ 3. Delinquency-themed variables

# ------ 3.1 Which time window length is the best in calculating delinquency volatility?

# - Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: g0_Delinq_SD_5, g0_Delinq_SD_4, g0_Delinq_SD_6, g0_Delinq_SD_9, g0_Delinq_SD_12
# Best C-statistics: g0_Delinq_SD_4, g0_Delinq_SD_5, g0_Delinq_SD_6, g0_Delinq_SD_9, g0_Delinq_SD_12

### CONCLUSION: No statistical significance, though it seems shorter periods are better than longer ones, given AIC.
# Choose first three for stepwise selection, based on AIC-differences



# ------ 3.2 Which lag order is the best in calculating the portfolio-level fraction of the non-defaulted proportion with any delinquency?

# - Initialize variables to be tested
vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Any_Aggr_Prop_Lag_1", "g0_Delinq_Any_Aggr_Prop_Lag_2",
          "g0_Delinq_Any_Aggr_Prop_Lag_3", "g0_Delinq_Any_Aggr_Prop_Lag_4", "g0_Delinq_Any_Aggr_Prop_Lag_5",
          "g0_Delinq_Any_Aggr_Prop_Lag_6", "g0_Delinq_Any_Aggr_Prop_Lag_9", "g0_Delinq_Any_Aggr_Prop_Lag_12" )

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: g0_Delinq_Any_Aggr_Prop_Lag_2, g0_Delinq_Any_Aggr_Prop_Lag_1, g0_Delinq_Any_Aggr_Prop_Lag_3, g0_Delinq_Any_Aggr_Prop
# Best C-statistics: g0_Delinq_Any_Aggr_Prop_Lag_1, g0_Delinq_Any_Aggr_Prop_Lag_2, g0_Delinq_Any_Aggr_Prop, g0_Delinq_Any_Aggr_Prop_Lag_3

### CONCLUSION: Earlier lags are better, though the AIC-differences were minuscule, as are concordance-differences
# Choose top 3 for stepwise selection



# ------ 3.3 Which lag order is the best in calculating the portfolio-level fraction of defaulted accounts?

# - Initialize variables to be tested
vars <- c("DefaultStatus1_Aggr_Prop", "DefaultStatus1_Aggr_Prop_Lag_1", "DefaultStatus1_Aggr_Prop_Lag_2",
          "DefaultStatus1_Aggr_Prop_Lag_3", "DefaultStatus1_Aggr_Prop_Lag_4", "DefaultStatus1_Aggr_Prop_Lag_5",
          "DefaultStatus1_Aggr_Prop_Lag_6", "DefaultStatus1_Aggr_Prop_Lag_9", "DefaultStatus1_Aggr_Prop_Lag_12" )

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: DefaultStatus1_Aggr_Prop, DefaultStatus1_Aggr_Prop_Lag_1, DefaultStatus1_Aggr_Prop_Lag_2, etc.
# Best C-statistics: DefaultStatus1_Aggr_Prop, DefaultStatus1_Aggr_Prop_Lag_1, DefaultStatus1_Aggr_Prop_Lag_2, etc.

### CONCLUSION: Earlier is better, though the AIC-differences were minuscule. Concordance-differences were greater
# Choose top 3 for stepwise selection



# ------ 3.4 How do other portfolio-level delinquency-themed variables fare as single-factor models?

# - Initialize variables to be tested
vars <- c("g0_Delinq_Ave", "ArrearsToBalance_Aggr_Prop", "CuringEvents_Aggr_Prop" )

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: ArrearsToBalance_Aggr_Prop, g0_Delinq_Ave, CuringEvents_Aggr_Prop
# Best C-statistics: g0_Delinq_Ave, ArrearsToBalance_Aggr_Prop, CuringEvents_Aggr_Prop

### CONCLUSION: Use Top 2, leave CuringEvents_Aggr_Prop



# ------ 3.5 How do account-level delinquency-themed variables fare as single-factor models?

# - Initialize variables to be tested
vars <- c("TimeInPerfSpell", "g0_Delinq_fac", "g0_Delinq", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3_Change_Ind",
          "g0_Delinq_Num", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "PerfSpell_g0_Delinq_Num", "Arrears")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: g0_Delinq, g0_Delinq_fac, g0_Delinq_Lag_1, TimeInDelinqState, slc_acct_arr_dir_3, Arrears, slc_acct_roll_ever_24_imputed_mean, slc_acct_arr_dir_3_Change_Ind, g0_Delinq_Num 
# Best C-statistics: g0_Delinq_fac, TimeInDelinqState, g0_Delinq_Lag_1, Arrears, slc_acct_arr_dir_3, slc_acct_roll_ever_24_imputed_mean, g0_Delinq_Num, slc_acct_arr_dir_3_Change_Ind

### CONCLUSION: Discard g0_Delinq given its deep relation to the outcome; rather use variants.
# Choose: TimeInDelinqState, g0_Delinq_Lag_1, slc_acct_arr_dir_3, slc_acct_roll_ever_24_imputed_mean, g0_Delinq_Num




# ------ 3.6 Combining insights: delinquency-themed variables

# - Initialize variables to be tested
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned",
          # Portfolio-level (Delinquency-themed)
          "g0_Delinq_SD_4", "g0_Delinq_SD_5","g0_Delinq_SD_6",
          "g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Any_Aggr_Prop_Lag_1","g0_Delinq_Any_Aggr_Prop_Lag_2",
          "DefaultStatus1_Aggr_Prop", "DefaultStatus1_Aggr_Prop_Lag_1", "DefaultStatus1_Aggr_Prop_Lag_2",
          "g0_Delinq_Ave", "ArrearsToBalance_Aggr_Prop",
          # Account-level (Delinquency-themed)
          "TimeInDelinqState", "g0_Delinq_Lag_1","slc_acct_arr_dir_3",
          "slc_acct_roll_ever_24_imputed_mean","g0_Delinq_Num")

# - Full model | Stepwise forward selection procedure
modLR_full <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
                  data=datCredit_train, family="binomial")
summary(modLR_full);
evalLR(modLR_full, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   19,736; McFadden R^2:  88.93%; AUC:  99.96%.
# Warning: model did not converge, so performance measures are not reliable

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
modLR_step <- stepAIC(modLR_base2, scope = list(lower = ~ log(TimeInPerfSpell)*PerfSpell_Num_binned, 
                                                    upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                      direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(modLR_step)
evalLR(modLR_step, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
proc.time() - ptm # IGNORE: elapsed runtime; 147m
### RESULTS: AIC:   20,174; McFadden R^2:  88.68%; AUC:  99.96%.

# - Domain expertise
# Remove g0_Delinq_SD_5 and g0_Delinq_SD_6 since g0_Delinq_SD_4 is already present
# Remove g0_Delinq_Any_Aggr_Prop_Lag_1 and g0_Delinq_Any_Aggr_Prop_Lag_2 since g0_Delinq_Any_Aggr_Prop is already present
# Remove g0_Delinq_Any_Aggr_Prop due to statistical insignificance

# - Final variables
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned", "g0_Delinq_SD_4", "g0_Delinq_Ave", 
          "TimeInDelinqState", "g0_Delinq_Lag_1","slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean")
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
                   data=datCredit_train, family="binomial")
summary(modLR);
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   21,407; McFadden R^2:  87.98%; AUC:  99.96%.





# ------ 4. Other portfolio-level variables

# ------ 4.1 Which lag order is the best in calculating the median interest rate of the portfolio?

# - Initialize variables to be tested
vars <- c("InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_1", "InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_3", "InterestRate_Margin_Aggr_Med_9")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: InterestRate_Margin_Aggr_Med , InterestRate_Margin_Aggr_Med_1, InterestRate_Margin_Aggr_Med_2, etc.
# Best C-statistics: InterestRate_Margin_Aggr_Med , InterestRate_Margin_Aggr_Med_1, InterestRate_Margin_Aggr_Med_2, etc.

### CONCLUSION: Earlier lags are better. Differences in AIC and concordance are minuscule.



# ------ 4.2 How do other portfolio-level (non-delinquency) variables fare as single-factor models?

# - Initialize variables to be tested
vars <- c("InstalmentToBalance_Aggr_Prop", "AgeToTerm_Aggr_Mean", "PerfSpell_Maturity_Aggr_Mean",
          "CreditLeverage_Aggr", "Ave_Margin_Aggr", "NewLoans_Aggr_Prop")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: AgeToTerm_Aggr_Mean, InstalmentToBalance_Aggr_Prop, NewLoans_Aggr_Prop, PerfSpell_Maturity_Aggr_Mean
# Best C-statistics: AgeToTerm_Aggr_Mean, InstalmentToBalance_Aggr_Prop, NewLoans_Aggr_Prop, PerfSpell_Maturity_Aggr_Mean

### CONCLUSION: Choose AgeToTerm_Aggr_Mean + InstalmentToBalance_Aggr_Prop + NewLoans_Aggr_Prop



# ------ 4.3 Combining insights: Delinquency-themed and portfolio-level variables

# - Initialize variables to be tested
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned", "g0_Delinq_SD_4", "g0_Delinq_Ave", 
          "TimeInDelinqState", "g0_Delinq_Lag_1","slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_1",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop")

# - Full model | Stepwise forward selection procedure
modLR_full <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
                   data=datCredit_train, family="binomial")
summary(modLR_full);
evalLR(modLR_full, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   21,299; McFadden R^2:  88.05%; AUC:  99.96%.
# Warning: model did not converge

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
modLR_step <- stepAIC(modLR_base2, scope = list(lower = ~ log(TimeInPerfSpell)*PerfSpell_Num_binned, 
                                                upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                      direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(modLR_step)
evalLR(modLR_step, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
proc.time() - ptm # IGNORE: elapsed runtime; 73m
### RESULTS: AIC:   21,311; McFadden R^2:  88.04%; AUC:  99.96%.

# Final variables
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned", "g0_Delinq_SD_4", #"g0_Delinq_Ave", 
          "TimeInDelinqState", "g0_Delinq_Lag_1","slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          #"InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_1",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop")
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial")
summary(modLR);
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   21,311; McFadden R^2:  88.04%; AUC:  99.96%.





# ------ 5. Account-level variables

# ------ 5.1 How do various non-delinquency account-level variables fare as single-factor models?

# - Initialize numeric variables to be tested for correlation
vars <- c("Principal_Real", "Principal", "InterestRate_Margin", # "pmnt_method_grp",
          "Balance_Real", "Balance", "Instalment_Real", "InterestRate_Nom", "AgeToTerm",
          "BalanceToPrincipal", "slc_acct_pre_lim_perc_imputed_med")

# - Correlation analysis towards obtaining clusters of correlated variables
corrAnalysis(datCredit_train, vars, corrThresh = 0.6, method = 'spearman')
### RSULTS: Expected correlation found between Principal_real and Principal
# 3-way correlation between Principal_Real, Balance/Balance_Real and Instalment
# Expected correlation between Balance/Balance_Real and BalanceToPrincipal
# Correlation found between AgeToTerm and BalanceToPrincipal
# Correlation found between Balance and AgeToTerm

# - Initialize variables to be tested
vars <- c("Principal_Real", "Principal", "InterestRate_Margin", "pmnt_method_grp",
          "Balance_Real", "Balance", "Instalment_Real", "InterestRate_Nom", "AgeToTerm",
          "BalanceToPrincipal", "slc_acct_pre_lim_perc_imputed_med")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: pmnt_method_grp, slc_acct_pre_lim_perc_imputed_med, InterestRate_Nom, Principal, InterestRate_Margin , etc
# Best C-statistics: pmnt_method_grp, slc_acct_pre_lim_perc_imputed_med, BalanceToPrincipal, InterestRate_Nom , Principal, etc.

### CONCLUSION: AIC-differences suggest top 4, while C-differences suggest top 2 only (66-74%)
# Choose: pmnt_method_grp, slc_acct_pre_lim_perc_imputed_med, InterestRate_Nom, Principal, BalanceToPrincipal


# ------ 5.2 Combining insights: Delinquency-themed, portfolio-level, and account-level variables

# - Initialize variables to be tested
# Added g0_Delinq_Ave again given its strategic importance
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned", "g0_Delinq_SD_4", "g0_Delinq_Ave",
          "TimeInDelinqState", "g0_Delinq_Lag_1","slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med", "InterestRate_Nom", "Principal", "BalanceToPrincipal")

# - Full model | Stepwise forward selection procedure
modLR_full <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
                   data=datCredit_train, family="binomial")
summary(modLR_full);
evalLR(modLR_full, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   20,690   88.40%; McFadden R^2:  %; AUC:  99.96%.
# Warning: model did not converge

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
modLR_step <- stepAIC(modLR_base2, scope = list(lower = ~ log(TimeInPerfSpell)*PerfSpell_Num_binned, 
                                                upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                      direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(modLR_step)
evalLR(modLR_step, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
proc.time() - ptm # IGNORE: elapsed runtime; 175m
### RESULTS: AIC:   20,688;   88.40%; McFadden R^2:  %; AUC:  99.96%.

# - Domain expertise
# Removed TimeInDelinqState since it led to model instability
# Added g0_Delinq_Lag_1 given its strategic importance

# - Final variables
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned", "g0_Delinq_SD_4", "g0_Delinq_Lag_1",
          #"g0_Delinq_Ave", ,"slc_acct_pre_lim_perc_imputed_med", "TimeInDelinqState", 
          "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom", "Principal", "BalanceToPrincipal")
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial")
summary(modLR);
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   26,135;   85.33%; McFadden R^2:  %; AUC:  99.91%.




# ------ 6. Macroeconomic variables

# ------ 6.1 Which lag order is the best for: M_Repo_Rate

# - Initialize variables to be tested
vars <- c("M_Repo_Rate", "M_Repo_Rate_1 ", "M_Repo_Rate_2", "M_Repo_Rate_3", "M_Repo_Rate_6", "M_Repo_Rate_9", "M_Repo_Rate_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: M_Repo_Rate_6, M_Repo_Rate_3, M_Repo_Rate_2, M_Repo_Rate_9, etc.
# Best C-statistics: M_Repo_Rate_12, M_Repo_Rate_9, M_Repo_Rate_6, M_Repo_Rate_3, etc.

### CONCLUSION: The minuscle AIC-differences suggest mid-term lags are better. C-differneces are small toe (57%).
# Choose: M_Repo_Rate_6, M_Repo_Rate_3, M_Repo_Rate_9



# ------ 6.2 Which lag order is the best for: M_Inflation_Growth

# - Initialize variables to be tested
vars <- c("M_Inflation_Growth", "M_Inflation_Growth_1 ", "M_Inflation_Growth_2", "M_Inflation_Growth_3", 
          "M_Inflation_Growth_6", "M_Inflation_Growth_9", "M_Inflation_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: M_Inflation_Growth_6, M_Inflation_Growth_3, M_Inflation_Growth_2, M_Inflation_Growth_1, etc.
# Best C-statistics: M_Inflation_Growth_6, M_Inflation_Growth_3, M_Inflation_Growth_9, M_Inflation_Growth_2, M_Inflation_Growth_1, etc. 

### CONCLUSION: Earlier lags seem better, based on very small AIC-differences, whereas mid-term lags seem better
# according to C-differences (57%). Choose: M_Inflation_Growth_6, M_Inflation_Growth_3, M_Inflation_Growth_9



# ------ 6.3 Which lag order is the best for: M_RealGDP_Growth

# - Initialize variables to be tested
vars <- c("M_RealGDP_Growth", "M_RealGDP_Growth_1 ", "M_RealGDP_Growth_2", "M_RealGDP_Growth_3", 
          "M_RealGDP_Growth_6", "M_RealGDP_Growth_9", "M_RealGDP_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: M_RealGDP_Growth_12, M_RealGDP_Growth_9, M_RealGDP_Growth_6, M_RealGDP_Growth_3, etc.
# Best C-statistics: M_RealGDP_Growth_12, M_RealGDP_Growth_9, M_RealGDP_Growth_6, M_RealGDP_Growth_3, etc.

### CONCLUSION: Later lags seem better, based on very small AIC-differences and c-differences (57-54%)
# Choose: M_RealGDP_Growth_12, M_RealGDP_Growth_9, M_RealGDP_Growth_6



# ------ 6.4 Which lag order is the best for: M_RealIncome_Growth

# - Initialize variables to be tested
vars <- c("M_RealIncome_Growth", "M_RealIncome_Growth_1 ", "M_RealIncome_Growth_2", "M_RealIncome_Growth_3", 
          "M_RealIncome_Growth_6", "M_RealIncome_Growth_9", "M_RealIncome_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: M_RealIncome_Growth_12, M_RealIncome_Growth_9, M_RealIncome_Growth_6, M_RealIncome_Growth_3, etc.
# Best C-statistics: M_RealIncome_Growth_12, M_RealIncome_Growth_9, M_RealIncome_Growth_6, M_RealIncome_Growth_3, etc.

### CONCLUSION: Later lags seem better, based on very small AIC-differences and c-differences (55%-53%)



# ------ 6.5 Which lag order is the best for: M_DTI_Growth

# - Initialize variables to be tested
vars <- c("M_DTI_Growth", "M_DTI_Growth_1 ", "M_DTI_Growth_2", "M_DTI_Growth_3", 
          "M_DTI_Growth_6", "M_DTI_Growth_9", "M_DTI_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: M_DTI_Growth, M_DTI_Growth_1, M_DTI_Growth_2, M_DTI_Growth_9, M_DTI_Growth_6, etc.
# Best C-statistics: M_DTI_Growth, M_DTI_Growth_1, M_DTI_Growth_9, M_DTI_Growth_2, M_DTI_Growth_12

### CONCLUSION: Shorter lags seem better, based on very small AIC-differences, whereas no discnerable
# trend exist across c-differences (58%)
# Choose: M_DTI_Growth, M_DTI_Growth_1, M_DTI_Growth_2



# ------ 6.6 Which lag order is the best for: M_Emp_Growth

# - Initialize variables to be tested
vars <- c("M_Emp_Growth", "M_Emp_Growth_1 ", "M_Emp_Growth_2", "M_Emp_Growth_3", 
          "M_Emp_Growth_6", "M_Emp_Growth_9", "M_Emp_Growth_12")

# - Single-factor modelling results
# Goodness-of-fit
aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Discriminatory power (in-sample)
concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
### RESULTS: Best AIC-results: M_Emp_Growth_12, M_Emp_Growth_9, M_Emp_Growth_6, M_Emp_Growth_3
# Best C-statistics: M_Emp_Growth_12, M_Emp_Growth_9, M_Emp_Growth_6, M_Emp_Growth_3

### CONCLUSION: Longer lags seem better, based on very small AIC-differnces, affirmed by the c-differences (56%-55%)



# ------ 6.7 Combining insights: Macroeconomic variables

# - Initialize variables to be tested
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned",
          "M_Repo_Rate_6", "M_Repo_Rate_3", "M_Repo_Rate_9", 
          "M_Inflation_Growth_6", "M_Inflation_Growth_3", "M_Inflation_Growth_9",
          "M_RealGDP_Growth_12", "M_RealGDP_Growth_9", "M_RealGDP_Growth_6",
          "M_RealIncome_Growth_12", "M_RealIncome_Growth_9", "M_RealIncome_Growth_6",
          "M_DTI_Growth", "M_DTI_Growth_1", "M_DTI_Growth_2", 
          "M_Emp_Growth_12", "M_Emp_Growth_9", "M_Emp_Growth_6")

# - Full model | Stepwise forward selection procedure
modLR_full <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
                   data=datCredit_train, family="binomial")
summary(modLR_full);
evalLR(modLR_full, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  176,311; McFadden R^2:  0.94%; AUC:  58.88%.

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
modLR_step <- stepAIC(modLR_base2, scope = list(lower = ~ log(TimeInPerfSpell)*PerfSpell_Num_binned, 
                                                upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                      direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(modLR_step)
evalLR(modLR_step, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
proc.time() - ptm # IGNORE: elapsed runtime; 34m
### RESULTS: AIC:   169,531;   McFadden R^2:  4.73%; AUC:  68.60%.

# - Domain expertise
# Removed M_DTI_Growth_2 since its unlagged variant performed better in single-factor models

# - Final variables
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned",
          "M_Repo_Rate_3", "M_Inflation_Growth_6","M_DTI_Growth", "M_Emp_Growth_6")
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial")
summary(modLR);
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   169,562;  McFadden R^2:  4.72%; AUC:  68.56%.



# ------ 6.8 Combining insights: Delinquency-themed, portfolio-level, account-level, and macroeconomic variables

# - Initialize variables to be tested
vars <- c("log(TimeInPerfSpell)*PerfSpell_Num_binned", "g0_Delinq_SD_4", "g0_Delinq_Lag_1",
          #"g0_Delinq_Ave","slc_acct_pre_lim_perc_imputed_med", "TimeInDelinqState", 
          "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom", "Principal", "BalanceToPrincipal",
          "M_Repo_Rate_3", "M_Inflation_Growth_6","M_DTI_Growth", "M_Emp_Growth_6")

# - Full model | Stepwise forward selection procedure
modLR_full <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
                   data=datCredit_train, family="binomial")
summary(modLR_full);
evalLR(modLR_full, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:   26,003;   McFadden R^2:  85.41%; AUC:  99.91%.

# - Stepwise forward selection using BIC
ptm <- proc.time() # for runtime calculations (ignore)
modLR_step <- stepAIC(modLR_base2, scope = list(lower = ~ log(TimeInPerfSpell)*PerfSpell_Num_binned, 
                                                upper = as.formula(paste("~", paste(vars, collapse = " + ")))), 
                      direction = "both", k=log(datCredit_train[,.N]), maxit=50)
summary(modLR_step)
evalLR(modLR_step, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
proc.time() - ptm # IGNORE: elapsed runtime; 117m
### RESULTS: AIC:  26,008;  McFadden R^2:  85.41%; AUC:  99.91%.

# - Domain expertise
# Exchanged log(TimeInPerfSpell) with Time_Binned for final modelling iteration
# Removed intercept term
# Reweighted default cases
# Removed BalanceToPrincipal based on robust sandwich estimator's summary (given weighting of cases)

# - Weigh default cases heavier. as determined interactively based on calibration success (script 6e)
datCredit_train[, Weight := ifelse(DefaultStatus1==1,10,1)]

# - Final variables
vars <- c("-1", "Time_Binned*PerfSpell_Num_binned", #"log(TimeInPerfSpell):PerfSpell_Num_binned",
          "g0_Delinq_SD_4", "g0_Delinq_Lag_1", 
          #"g0_Delinq_Ave","slc_acct_pre_lim_perc_imputed_med", "TimeInDelinqState", #"Principal", "M_Emp_Growth_6"
          "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom",
          "M_Inflation_Growth_6","M_DTI_Growth")
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial", weights = Weight)
#summary(modLR);
# Robust (sandwich) standard errors
robust_se <- vcovHC(modLR, type="HC0")
# Summary with robust SEs
coeftest(modLR, vcov.=robust_se)

# Other diagnostics
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  86,046;  McFadden R^2:  51.70%; AUC:  99.94%.





# ------ 6. Final model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]
# remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Weigh default cases heavier. as determined interactively based on calibration success (script 6e)
datCredit_train[, Weight := ifelse(DefaultStatus1==1,10,1)]

# - Fit an "empty" model as a performance gain, used within some diagnostic functions
modLR_base <- glm(PerfSpell_Event ~ 1, data=datCredit_train, family="binomial")

# - Final variables
vars <- c("-1", "Time_Binned*PerfSpell_Num_binned", #"log(TimeInPerfSpell):PerfSpell_Num_binned",
  "g0_Delinq_SD_4", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
  "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
  "pmnt_method_grp", "InterestRate_Nom",
  "M_Inflation_Growth_6","M_DTI_Growth")
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial", weights = Weight)
#summary(modLR);
# Robust (sandwich) standard errors
robust_se <- vcovHC(modLR, type="HC0")
# Summary with robust SEs
coeftest(modLR, vcov.=robust_se)

# - Other diagnostics
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  86,046;  McFadden R^2:  51.70%; AUC:  99.94%.

# - Test goodness-of-fit using AIC-measure over single-factor models
aicTable_CoxDisc <- aicTable(datCredit_train, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Top variables: g0_Delinq_SD_4, g0_Delinq_Lag_1, slc_acct_arr_dir_3, slc_acct_roll_ever_24_imputed_mean, pmnt_method_grp, Time_Binned*PerfSpell_Num_binned

# Test accuracy using c-statistic over single-factor models
concTable_CoxDisc <- concTable(datCredit_train, datCredit_valid, vars, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Top variables: g0_Delinq_SD_4, g0_Delinq_Lag_1, slc_acct_arr_dir_3, slc_acct_roll_ever_24_imputed_mean, pmnt_method_grp, Time_Binned*PerfSpell_Num_binned

# - Combine results into a single object
Table_CoxDisc <- concTable_CoxDisc[,1:2] %>% left_join(aicTable_CoxDisc, by ="Variable")

# Save objects
pack.ffdf(paste0(genObjPath,"CoxDisc_advanced_fits"), Table_CoxDisc)

