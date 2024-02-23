# ====================================== Cox PH Model (Quasi-Complete Separation) ======================================
# Fitting Cox Proportional Hazard (PH) models using g0_Delinq, which is a variable known to cause complete/
# quasi-complete separation issues.
# -----------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller
# -----------------------------------------------------------------------------------------
# -- Script dependencies:
#   -  0a.Setup
#   - 2d.Data_Fusion.R

# -- Inputs:
#   - datCredit_real | Prepared from script 2d.
#
# -- Outputs:
#   - Analytics on Cox PH models containing various instances of the g0_Delinq variable. 
# -----------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# --- Loaing in the required dataset
if (!exists('dat_train')) unpack.ffdf(paste0(genPath,"dat_train"), tempPath)

# --- Remove some unnecessary variables to improve performance
colnames_dat_train <- colnames(dat_train)
dat_train <- dat_train %>% subset(select = -c(which(colnames_dat_train == "M_Inflation_Growth_1"):which(colnames_dat_train == "M_Repo_Rate_Vol_12")))

# --- Subsampling
# - Setting the seed
set.seed(123)
full_size <- 1000000
full_prop <- full_size/nrow(dat_train)

# First obtain the first observation of each performance spell in the dataset (in doing so we obtain all the unique PerfSpell_Keys)
dat_temp <- dat_train[PerfSpell_Counter==1, c("Date", "PerfSpell_Key", "PerfSpellResol_Type_Hist")]
# Now sample from the unique IDs
dat_full_keys <- dat_temp %>% slice_sample(prop=full_prop)
# Resample using full dataset with the PerfSpell_Keys of the temp dataset
dat_train2 <- dat_train[PerfSpell_Key %in% dat_full_keys$PerfSpell_Key]; rm(dat_train, dat_temp, dat_full_keys)

# --- Small sanity check
dat_train$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
### RESULTS:~ Censored = 45.29%; Defaulted = 12.13%; Paid-up = 0.53%; Settled = 41.76%; Written-off = 0.29%
dat_train2$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
### RESULTS:~ Censored = 45.44%; Defaulted = 12.10%; Paid-up = 0.47%; Settled = 41.68%; Written-off = 0.30%

# --- Setting the seed to ensure reporducability
set.seed(123)

# --- Experiments 1
# - Random assignment of levels where g0_Delinq=3
#  Applying random sampling to instances where g0_Delinq==3 with values c(0,1,2) with uniform probabilities (1/3)
dat_train2[, g0_Delinq2a := ifelse(g0_Delinq!=3, g0_Delinq, sample(x = 0:2, size=1, prob=rep(1/3, times = 3)))]
#  Applying random sampling to instances where g0_Delinq==3 with values c(0,1,2) with prior probabilities of g0_Delinq
g0_Outcomes <- table(dat_train2$g0_Delinq) %>% prop.table() # vector containnig the prior probabilities
dat_train2[, g0_Delinq2b := ifelse(g0_Delinq!=3, g0_Delinq, sample(0:2, size=1, prob=g0_Outcomes[1:3]))]


# - Assignment of g0_Delinq=3 with mean and median imputation
# Assigning the median value of g0_Delinq to instances where g0_Delinq==3
g0_Median <- median(dat_train2[g0_Delinq != 3, g0_Delinq])
dat_train2[ , g0_Delinq3 := ifelse(g0_Delinq!=3, g0_Delinq, g0_Median)]
# Assiging the mean value of g0_Delinq to instances where g0_Delinq==3
g0_Mean <- mean(dat_train2[g0_Delinq != 3, g0_Delinq])
dat_train2[ , g0_Delinq4 := ifelse(g0_Delinq!=3, g0_Delinq, g0_Mean)]


# - Binning Schemes
# Binning values where g0_Delinq==3 with g0_Delinq==2
dat_train2[, g0_Delinq5 := ifelse(g0_Delinq!=3, g0_Delinq, 2)]
# Binning values where g0_Delinq==2 or where g0_Delinq==3 with g0_Delinq==1
dat_train2[, g0_Delinq6 := ifelse(g0_Delinq %in% c(1,2,3), 1, 0)]


# - Indicator variables
# Indicator variable for each level
dat_train2[, g0_Delinq7_1 := ifelse(g0_Delinq == 1, 1, 0)]
dat_train2[, g0_Delinq7_2 := ifelse(g0_Delinq == 2, 1, 0)]
dat_train2[, g0_Delinq7_3 := ifelse(g0_Delinq == 3, 1, 0)]
# Indicator variable for each level, where levels 2 and 3 are grouped together
dat_train2[, g0_Delinq7_23 := ifelse(g0_Delinq %in% c(2,3), 1, 0)]
# Indicator variables for levels 1, 2, and 3 grouped and 0 alone
dat_train2[, g0_Delinq7_123 := ifelse(g0_Delinq %in% c(1,2,3), 1, 0)]


# - Facotorise g0_Delinq
# Factorise the raw variable
dat_train2[, g0_Delinq8a := factor(g0_Delinq)]
# Factorise the binned variable
dat_train2[, g0_Delinq8b := factor(g0_Delinq5)]


# --- Experiments 4
# - State specific variables (raw variables)
# Indicator for when a shift in the state of g0_Delinq occurs (target event) - special case for the first record of an account
dat_train2[, StateShift := ifelse(lag(g0_Delinq5, n=1)==g0_Delinq5,0,1), by=PerfSpell_Key]
dat_train2[, StateShift := ifelse(is.na(StateShift)==T, 1, StateShift)]
# State number
dat_train2[, State_Num := cumsum(StateShift), by=PerfSpell_Key]
# State key (unique ID per state)
dat_train2[, State_Key := paste0(PerfSpell_Key, "_", State_Num)]
# State counter (Time in state)
dat_train2[ , TimeInState := 1:.N, by=State_Key]
# Performance spell example: 3000003205066_2
# lookup <- dat_train2[PerfSpell_Key=="3000003205066_2", list(PerfSpell_Key, g0_Delinq5, TimeInPerfSpell, StateShift, State_Num, TimeInState,TimeInState)]

# - State specific variables (direction of transitions embedded)
# Indicator for when a shift in the state of g0_Delinq occurs (target event) - special case for the first record of an account
dat_train2[, StateShift2 := ifelse(lag(g0_Delinq5, n=1)==g0_Delinq5, 0, ifelse(lag(g0_Delinq5)<g0_Delinq5, -1, 1)), by=PerfSpell_Key]
dat_train2[, StateShift2 := ifelse(is.na(StateShift2)==T, 1, StateShift2)]
# Direction (Numeric) of the shifts in the state (value aggregated per state)
dat_train2[, StateDir := sum(StateShift2), by=State_Key]
# Direction (Numeric) of the shifts in the state (value aggregated per state)
dat_train2[, StateDir2 := as.factor(StateDir), by=State_Key]
# lookup <- cbind(lookup, dat_train2[PerfSpell_Key=="3000003205066_2", list(StateShift2, StateDir)])

# - Volatility of delinquency transitions
# Spell level volatility
dat_train2[, StateVol := sd(State_Num), by = PerfSpell_Key] # Spell level
# 3-month volatility
dat_train2[, StateVol_3_Months := imputeFirstKnown(frollapply(x=State_Num, n=3, align="right", FUN=sd)), by = PerfSpell_Key] # Period level
# 6-month volatility
dat_train2[, StateVol_6_Months := imputeFirstKnown(frollapply(x=State_Num, n=6, align="right", FUN=sd)), by = PerfSpell_Key] # Period level
# lookup <- cbind(lookup, dat_train2[PerfSpell_Key == "3000003205066_2", list(StateVol, StateVol_3_Months, StateVol_6_Months)])


# --- Experimentation 5
# - Lagging delinquency
# 1 Month lag
dat_train2[, g0_Delinq5_lag_1 := lag(g0_Delinq5, n=1), by=PerfSpell_Key]
# 2 Month lag
dat_train2[, g0_Delinq5_lag_2 := lag(g0_Delinq5, n=2), by=PerfSpell_Key]
# 3 Month lag
dat_train2[, g0_Delinq5_lag_3 := lag(g0_Delinq5, n=3), by=PerfSpell_Key]

# lookup <- cbind(lookup, dat_train2[PerfSpell_Key=="3000003205066_2", list(g0_Delinq5_lag_1, g0_Delinq5_lag_2, g0_Delinq5_lag_3)])

# - Lagging delinquency factorised
# 1 Month lag
dat_train2[, g0_Delinq5_lag_1_f := ifelse(is.na(g0_Delinq5_lag_1), "MISSING", g0_Delinq5_lag_1)]
dat_train2[, g0_Delinq5_lag_1_f := factor(g0_Delinq5_lag_1_f), by=PerfSpell_Key]
# 2 Month lag
dat_train2[, g0_Delinq5_lag_2_f := ifelse(is.na(g0_Delinq5_lag_2), "MISSING", g0_Delinq5_lag_2)]
dat_train2[, g0_Delinq5_lag_2_f := factor(g0_Delinq5_lag_2_f), by=PerfSpell_Key]
# 3 Month lag
dat_train2[, g0_Delinq5_lag_3_f := ifelse(is.na(g0_Delinq5_lag_3), "MISSING", g0_Delinq5_lag_3)]
dat_train2[, g0_Delinq5_lag_3_f := factor(g0_Delinq5_lag_3_f), by=PerfSpell_Key]

# lookup <- cbind(lookup, dat_train2[PerfSpell_Key=="3000003205066_2", list(g0_Delinq5_lag_1_f, g0_Delinq5_lag_2_f, g0_Delinq5_lag_3_f)])


# --- Final Preparation
# - Transforming the Default status into a numeric variable to ensure compatibility with time-dependent AUC/ROC functions
dat_train2[, DefaultStatus1 := as.numeric(DefaultStatus1)]; dat_train2[,DefaultStatus1 := ifelse(DefaultStatus1==1,0,1)]
# - Grouping the data according to PerfSpell_Key and TimeInPerfSpell
dat_train2 <- dat_train2 %>% group_by(PerfSpell_Key, TimeInPerfSpell)




# ------ 2. Cox PH Models - Experimentation 1
# --- Raw variable
cph_def_exp1_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq,
                           data=dat_train2, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq: coef = 28.35; exp(coef) = 5046920172101.94; se(coef) = 9655.60; z = 0.003; p = 0.998

# --- Binning Schemes
# - Variable where all g0_Delinq=3 are assigned values of two.
cph_def_exp1_g0_5 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5,
                           data=dat_train2, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq5: coef = 22.23; exp(coef) = 4520192000.62; se(coef) = 478.31; z = 0.046; p = 0.963

# ---  Indicator variables
# - Indicators for each level of g0_Delinq3
cph_def_exp1_g0_7a <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq7_1 + g0_Delinq7_2 + g0_Delinq7_3,
                            data=dat_train2, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  3 ; beta may be infinite.
# g0_Delinq7_1: coef = 0.05598; exp(coef) = 1.05758; se(coef) = 5.17113; z = 0.011; p = 0.991
# g0_Delinq7_2: coef = 0.18407; exp(coef) = 1.20210; se(coef) = 13.75389; z = 0.013; p = 0.989
# g0_Delinq7_3: coef = 408.94868; exp(coef) = 401933698870632771984282400084...; se(coef) = 1.48428; z = 275.520; p < 0.0000000000000002

# - Indacator for each level of g0_Delinq3, where levels 2 and 3 are grouped.
cph_def_exp1_g0_7b <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq7_1 + g0_Delinq7_23,
                            data=dat_train2, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  2 ; beta may be infinite.
# g0_Delinq7_1 : coef = 0.01629; exp(coef) = 1.01642; se(coef) = 1730.46647; z = 0.000; p  = 1.000
# g0_Delinq7_23: coef = 24.76415; exp(coef) = 56876784248.29054; se(coef) = 402.84324; z = 0.061; p = 0.951




# ------ 3. Cox PH Models - Experimentation 2
# --- Combinig g0_Delinq5 (binned variable) with other covariates in attempting to achieve convergence
# - Fitting a Cox model with the binned delinquency variable and an additional covariate in an attempt to achieve convergence
cph_def_exp2_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5 + Principal,
                           data=dat_train2, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq5: coef = 25.36318817884; exp(coef) = 103536306837.56080627441; se(coef) = 2291.50561647697; z = 0.011; p = 0.991
# Principal : coef = 0.00000006210; exp(coef) = 1.00000006210; se(coef) = 0.00000003878; z = 1.602; p = 0.109

# - Fitting a Cox model with the binned delinquency variable and two additional covariates in an attempt to achieve convergence
cph_def_exp2_g0_2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5 + Principal + Instalment,
                           data=dat_train2, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq5: coef = 25.3616350926; exp(coef) = 103375630825.7816467285; se(coef) = 2294.5502774228; z = 0.011; p = 0.991
# Principal : coef = 0.0000004070; exp(coef) = 1.0000004070; se(coef) = 0.0000000684; z = 5.951; p = 0.00000000266
# Instalment: coef = -0.0000401470; exp(coef) = 0.9999598538; se(coef) = 0.0000074914; z = -5.359; p = 0.00000008366




# ------ 3. Cox PH Models - Experimentation 3
# --- Combinig g0_Delinq8a (factorised/ dummy encoded variable) with other covariates in attempting to achieve convergence
#     Because of the non-linear hazards expected by different levels of g0_Delinq, fitting a model with a facotrised variable
#     will solve the non-linearity and might aid convergence.
# - Fitting a Cox model with the factorised variable and an additional covariate in an attempt to achieve convergence
cph_def_exp3_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq8b + Principal,
                           data=dat_train2, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8b_1: coef = 0.02848757612; exp(coef) = 1.02889722785; se(coef) = 1731.62955126728; z = 0.000; p = 1.000
# g0_Delinq8b_2: coef = 24.77558329928; exp(coef) = 57530627739.88069152832; se(coef) = 0.061; z = 0.377; p = 0.951
# Principal    : coef = 0.00000006210; exp(coef) = 1.00000006210; se(coef) = 0.00000003878; z = 1.602; p = 0.109

# - Fitting a Cox model with the factorised variable and two additional covariates in an attempt to achieve convergence
cph_def_exp3_g0_2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq8b + Principal + Instalment,
                           data=dat_train2, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8b_1: coef = 0.0448162562; exp(coef) = 1.0458356765; se(coef) = 1733.8218679946; z = 0; p = 1.000
# g0_Delinq8b_2: coef = 24.7907217431; exp(coef) = 58408177524.0331039429; se(coef) = 402.2316254501; z = 0.062; p = 0.951
# Principal    : coef = 0.0000004070; exp(coef) = 1.0000004070; se(coef) = 0.0000000684; z = 5.951; p = 0.00000000266
# Instalment   : coef =  -0.0000401470; exp(coef) = 0.9999598538 ; se(coef) = 0.0000074914; z = -5.359; p = 0.00000008366




# ------ 5. Cox PH Models - Experimentation 5
### Binned delinquency variable lagged by different number of months
# --- 1-Month lag (Numerical)
### NOTE:~ 457 6245 observations deleted due to missingness ~ 457624/29071039 = 1.574158% of observations
# - 1-Month lag
cph_def_exp5_g0_1_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1,
                             data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState
cph_def_exp5_g0_1_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                g0_Delinq5_lag_1 + TimeInState,
                              data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState and StateVol
cph_def_exp5_g0_1_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                 g0_Delinq5_lag_1 + TimeInState + StateVol,
                               data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# --- 1-Month lag (categorical)
# - 1-Month lag
cph_def_exp5_g0_2_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1_f,
                             data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_f2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_fMISSING : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag and TimeInState
cph_def_exp5_g0_2_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                g0_Delinq5_lag_1_f + TimeInState,
                              data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_f2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_fMISSING : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# --- 1-Month lag (numerical) and binned delinquency
cph_def_exp5_g0_3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1 + g0_Delinq5,
                             data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# --- 2-Month lag (numerical) and binned delinquency
# - 2-Month lag with binned delinquency
cph_def_exp5_g0_4 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_2 + g0_Delinq5,
                             data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_2 : coef = 0.05689; exp(coef) = 1.05854; se(coef) = 0.03645; z = 1.561; p = 0.119
# g0_Delinq5       : coef = 23.57859; exp(coef) = 17380030625.37346; se(coef) = 963.98020; z = 0.024; p = 0.980

# --- Combining lags
# - 1-Month and 2-month lag (numerical)
cph_def_exp5_g0_5_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                g0_Delinq5_lag_1 + g0_Delinq5_lag_2,
                              data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1       : coef = 4.31432; exp(coef) = 74.76279; se(coef) = 0.04607; z = 93.640; p < 0.0000000000000002
# g0_Delinq5_lag_2       : coef = -0.30560; exp(coef) = 0.73668; se(coef) = 0.04172; z = -7.324; p = 0.00000000000024

# - 1-Month and 2-month lag (categorical)
cph_def_exp5_g0_5_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                 g0_Delinq5_lag_1_f + g0_Delinq5_lag_2_f,
                               data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_2_f       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# --- Combined lags with increased input space
# - 1-Month and 2-month lag (numerical) with Principal
cph_def_exp5_g0_6_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1 + g0_Delinq5_lag_2 + Principal,
                             data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_2         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal                : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month and 2-month lag (numerical) with TimeInState
cph_def_exp5_g0_6_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                g0_Delinq5_lag_1 + g0_Delinq5_lag_2 + TimeInState,
                              data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_2         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState              : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month and 2-month lag (numerical) with Principal and TimeInState
cph_def_exp5_g0_6_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                 g0_Delinq5_lag_1 + g0_Delinq5_lag_2 + Principal + TimeInState,
                               data=dat_train2, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_2         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal                : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState              : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
