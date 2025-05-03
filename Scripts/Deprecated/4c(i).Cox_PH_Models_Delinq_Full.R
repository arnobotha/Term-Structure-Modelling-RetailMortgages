# ====================================== Cox PH Model (Quasi-Complete Separation) ======================================
# Fitting Cox Proportional Hazard (PH) models using g0_Delinq, which is a variable known to cause quasi-complete
# separation issues.
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

# --- Setting the seed to ensure reporducability
set.seed(123)

# --- Experiments 1
# - Random assignment of levels whete g0_Delinq=3
#  Applying random sampling to instances where g0_Delinq==3 with values c(0,1,2) with uniform probabilities (1/3)
dat_train[, g0_Delinq2a := ifelse(g0_Delinq!=3, g0_Delinq, sample(x = 0:2, size=1, prob=rep(1/3, times = 3)))]
#  Applying random sampling to instances where g0_Delinq==3 with values c(0,1,2) with prior probabilities of g0_Delinq
g0_Outcomes <- table(dat_train$g0_Delinq) %>% prop.table() # vector containnig the prior probabilities
dat_train[, g0_Delinq2b := ifelse(g0_Delinq!=3, g0_Delinq, sample(0:2, size=1, prob=g0_Outcomes[1:3]))]


# - Assignment of g0_Delinq=3 with mean and median imputation
# Assigning the median value of g0_Delinq to instances where g0_Delinq==3
g0_Median <- median(dat_train[g0_Delinq != 3, g0_Delinq])
dat_train[ , g0_Delinq3 := ifelse(g0_Delinq!=3, g0_Delinq, g0_Median)]
# Assiging the mean value of g0_Delinq to instances where g0_Delinq==3
g0_Mean <- mean(dat_train[g0_Delinq != 3, g0_Delinq])
dat_train[ , g0_Delinq4 := ifelse(g0_Delinq!=3, g0_Delinq, g0_Mean)]


# - Binning Schemes
# Binning values where g0_Delinq==3 with g0_Delinq==2
dat_train[, g0_Delinq5 := ifelse(g0_Delinq!=3, g0_Delinq, 2)]
# Binning values where g0_Delinq==2 or where g0_Delinq==3 with g0_Delinq==1
dat_train[, g0_Delinq6 := ifelse(g0_Delinq %in% c(1,2,3), 1, 0)]


# - Indicator variables
# Indicator variable for each level
dat_train[, g0_Delinq7_1 := ifelse(g0_Delinq == 1, 1, 0)]
dat_train[, g0_Delinq7_2 := ifelse(g0_Delinq == 2, 1, 0)]
dat_train[, g0_Delinq7_3 := ifelse(g0_Delinq == 3, 1, 0)]
# Indicator variable for each level, where levels 2 and 3 are grouped together
dat_train[, g0_Delinq7_23 := ifelse(g0_Delinq %in% c(2,3), 1, 0)]
# Indicator variables for levels 1, 2, and 3 grouped and 0 alone
dat_train[, g0_Delinq7_123 := ifelse(g0_Delinq %in% c(1,2,3), 1, 0)]


# - Facotorise g0_Delinq
# Factorise the raw variable
dat_train[, g0_Delinq8a := factor(g0_Delinq)]
# Factorise the binned variable
dat_train[, g0_Delinq8b := factor(g0_Delinq5)]


# --- Experiments 4
# - State specific variables (raw variables)
# Indicator for when a shift in the state of g0_Delinq occurs (target event) - special case for the first record of an account
dat_train[, StateShift := ifelse(lag(g0_Delinq5, n=1)==g0_Delinq5,0,1), by=PerfSpell_Key]
dat_train[, StateShift := ifelse(is.na(StateShift)==T, 1, StateShift)]
# State number
dat_train[, State_Num := cumsum(StateShift), by=PerfSpell_Key]
# State key (unique ID per state)
dat_train[, State_Key := paste0(PerfSpell_Key, "_", State_Num)] # Exercise caution using this variable, since the spells are not correctly defined: There is no overlap between the end of an old spell and the begining of a new spell
# State counter (Time in state)
dat_train[ , TimeInState := 1:.N, by=State_Key] # Exercise caution using this variable, since the spells are not correctly defined: There is no overlap between the end of an old spell and the begining of a new spell
# Performance spell example: 3000003205066_2
# lookup <- dat_train[PerfSpell_Key=="3000003205066_2", list(PerfSpell_Key, g0_Delinq5, TimeInPerfSpell, StateShift, State_Num, TimeInState)]

# - State specific variables (direction of transitions embedded)
# Indicator for when a shift in the state of g0_Delinq occurs (target event) - special case for the first record of an account
dat_train[, StateShift2 := ifelse(lag(g0_Delinq5, n=1)==g0_Delinq5, 0, ifelse(lag(g0_Delinq5)<g0_Delinq5, -1, 1)), by=PerfSpell_Key]
dat_train[, StateShift2 := ifelse(is.na(StateShift2)==T, 1, StateShift2)]
# Direction (Numeric) of the shifts in the state (value aggregated per state)
dat_train[, StateDir := sum(StateShift2), by=State_Key]
# Direction (Numeric) of the shifts in the state (value aggregated per state)
dat_train[, StateDir2 := as.factor(StateDir), by=State_Key]
# lookup <- cbind(lookup, dat_train[PerfSpell_Key=="3000003205066_2", list(StateShift2, StateDir)])

# - Volatility of delinquency transitions
# Spell level volatility
dat_train[, StateVol := sd(g0_Delinq), by = PerfSpell_Key] # Spell level
# 3-month volatility
dat_train[, StateVol_3_Months := imputeFirstKnown(frollapply(x=g0_Delinq, n=3, align="right", FUN=sd)), by = PerfSpell_Key] # Period level
# 6-month volatility
dat_train[, StateVol_6_Months := imputeFirstKnown(frollapply(x=g0_Delinq, n=6, align="right", FUN=sd)), by = PerfSpell_Key] # Period level
# lookup <- cbind(lookup, dat_train[PerfSpell_Key == "3000003205066_2", list(StateVol, StateVol_3_Months, StateVol_6_Months)])


# --- Experimentation 5
# - Lagging delinquency
# 1 Month lag
dat_train[, g0_Delinq5_lag_1 := lag(g0_Delinq5, n=1), by=PerfSpell_Key]
# 2 Month lag
dat_train[, g0_Delinq5_lag_2 := lag(g0_Delinq5, n=2), by=PerfSpell_Key]
# 3 Month lag
dat_train[, g0_Delinq5_lag_3 := lag(g0_Delinq5, n=3), by=PerfSpell_Key]

# lookup <- cbind(lookup, dat_train[PerfSpell_Key=="3000003205066_2", list(g0_Delinq5_lag_1, g0_Delinq5_lag_2, g0_Delinq5_lag_3)])

# - Lagging delinquency factorised
# 1 Month lag
dat_train[, g0_Delinq5_lag_1_f := ifelse(is.na(g0_Delinq5_lag_1), "MISSING", g0_Delinq5_lag_1)]
dat_train[, g0_Delinq5_lag_1_f := factor(g0_Delinq5_lag_1_f), by=PerfSpell_Key]
# 2 Month lag
dat_train[, g0_Delinq5_lag_2_f := ifelse(is.na(g0_Delinq5_lag_2), "MISSING", g0_Delinq5_lag_2)]
dat_train[, g0_Delinq5_lag_2_f := factor(g0_Delinq5_lag_2_f), by=PerfSpell_Key]
# 3 Month lag
dat_train[, g0_Delinq5_lag_3_f := ifelse(is.na(g0_Delinq5_lag_3), "MISSING", g0_Delinq5_lag_3)]
dat_train[, g0_Delinq5_lag_3_f := factor(g0_Delinq5_lag_3_f), by=PerfSpell_Key]

# lookup <- cbind(lookup, dat_train[PerfSpell_Key=="3000003205066_2", list(g0_Delinq5_lag_1_f, g0_Delinq5_lag_2_f, g0_Delinq5_lag_3_f)])


# --- Final Preparation
# - Transforming the Default status into a numeric variable to ensure compatibility with time-dependent AUC/ ROC functions
dat_train[, DefaultStatus1 := as.numeric(DefaultStatus1)]; dat_train[,DefaultStatus1 := ifelse(DefaultStatus1==1,0,1)]
# - Grouping the data according to PerfSpell_Key and TimeInPerfSpell
dat_train <- dat_train %>% group_by(PerfSpell_Key, TimeInPerfSpell)




# ------ 2. Cox PH Models - Experimentation 1
# --- Raw variable
cph_def_exp1_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq: coef = 29.25; exp(coef) = 5046920172101.94; se(coef) = 1935.38; z = 0.015; p = 0.988


# --- Random assignment of levels whete g0_Delinq=3
# - Variable where random levels assigned to g0_Delinq=3; uniform probabilities (1/3).
cph_def_exp1_g0_2a <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq2a,
                            data=dat_train, id=PerfSpell_Key)
# g0_Delinq: coef = 2.425755; exp(coef) = 11.310767; se(coef) = 0.004012; z = 604.7; p < 0.0000000000000002

# - Variable where random levels assigned to g0_Delinq=3; prior probabilities of g0_Delinq
cph_def_exp1_g0_2b <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq2b,
                            data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq2b: coef = -16.77133865913; exp(coef) = 0.00000005204; se(coef) = 58.41565152545; z = -0.287; p = 0.774


# --- Assignment of g0_Delinq=3 with mean and median imputation
# - Variable where all g0_Delinq=3 are assigned the value of the median
cph_def_exp1_g0_3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq3,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq3: coef = -16.77133865913; exp(coef) = 0.00000005204; se(coef) = 58.41565152545; z = -0.287; p = 0.774

# - Variable where all g0_Delinq=3 are assigned the value of the mean
cph_def_exp1_g0_4 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq4,
                           data=dat_train, id=PerfSpell_Key)
# g0_Delinq4: coef = -0.1193; exp(coef) = 0.8875; se(coef) = 0.0124; z = -9.62; p < 0.0000000000000002

# --- Binning Schemes
# - Variable where all g0_Delinq=3 are assigned values of two.
cph_def_exp1_g0_5 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq5,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq5: coef = 21.60; exp(coef) = 2401330681.67; se(coef) = 64.49; z = 0.335; p = 0.738

cph_def_exp1_g0_6 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq6,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq5: coef = 22.40; exp(coef) = 5339121628.21; se(coef) = 64.79; z = 0.346; p = 0.73


# ---  Indicator variables
# - Indicators for each level of g0_Delinq3
cph_def_exp1_g0_7a <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq7_1 + g0_Delinq7_2 + g0_Delinq7_3,
                            data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  3 ; beta may be infinite.
# g0_Delinq7_1: coef = 0.004050; exp(coef) = 1.004059; se(coef) = 427.374545; z = 0; p = 1
# g0_Delinq7_2: coef = 0.009081; exp(coef) = 1.009123; se(coef) = 1178.554568; z = 0; p = 1
# g0_Delinq7_3: coef = 28.202998; exp(coef) = 1771765913517.191406; se(coef) = 100.667289; z = 0.28; p = 0.779

# - Indacator for each level of g0_Delinq3, where levels 2 and 3 are grouped.
cph_def_exp1_g0_7b <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq7_1 + g0_Delinq7_23,
                            data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  2 ; beta may be infinite.
# g0_Delinq7_1 : coef = 0.01646; exp(coef) = 1.01660; se(coef) = 286.48729; z = 0; p  = 1
# g0_Delinq7_23: coef = 24.54190; exp(coef) = 45541827534.90582; se(coef) = 67.17892; z = 0.365; p = 0.715

# - Indacator for levels 1, 2, and 3 grouped and 0 as another group.
cph_def_exp1_g0_7c <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq7_123,
                            data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  1 ; beta may be infinite.
# g0_Delinq7_123: coef = 22.40; exp(coef) = 5339121628.21; se(coef) = 64.79; z = 0.346; p = 0.73


# --- Factorised variable
# - Variable is factorised
cph_def_exp1_g0_8 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                          g0_Delinq8,
                          data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8_1: coef = 0.004050; exp(coef) = 1.004059; se(coef) = 427.374545; z = 0; p = 1
# g0_Delinq8_2: coef = 0.009081; exp(coef) = 1.009123; se(coef) = 1178.554568; z = 0; p = 1
# g0_Delinq8_3: coef = 28.202998; exp(coef) = 1771765913517.191406; se(coef) = 100.667289; z = 0.28; p = 0.779


# --- Penalised Cox Regression
# WARNING: VERY LONG RUNTIME (~ 1.5 hours +)
cph_def_exp1_g0_9 <- coxphf::coxphf(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                    g0_Delinq,
                                    data=dat_train)
# NOTE: No errors or warnings (run time is substantially longer than normal coxph())
# g0_Delinq: coef = 12.83037; exp(coef) = 373387.4; se(coef) = 1.414128; z = ?; p = 0




# ------ 3. Cox PH Models - Experimentation 2
# --- Combinig g0_Delinq5 (binned variable) with other covariates in attempting to achieve convergence
# - Fitting a Cox model with the binned delinquency variable and an additional covariate in an attempt to achieve convergence
cph_def_exp2_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq5 + Principal,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq5: coef = 16.633446891487; exp(coef) = 16742264.547202795744; se(coef) = 5.388288086162; z = 3.087; p = 0.00202
# Principal : coef = 0.000000066238; exp(coef) = 1.000000066238; se(coef) = 0.000000006491; z = 10.205; p < 0.0000000000000002

# - Fitting a Cox model with the binned delinquency variable and two additional covariates in an attempt to achieve convergence
cph_def_exp2_g0_2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5 + Principal + Instalment,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq5: coef = 12.055777558857; exp(coef) = 172090.806088245095; se(coef) = 0.529045029487; z = 22.788; p < 0.0000000000000002
# Principal : coef = 0.000000043040; exp(coef) = 1.000000043040; se(coef) = 0.000000006576; z = 6.545; p = 0.0000000000596
# Instalment: coef = 0.000002458192; exp(coef) = 1.000002458195; se(coef) = 0.000000096423; z = 25.494; p < 0.0000000000000002

# - Fitting a Cox model with the binned delinquency variable and three additional covariates in an attempt to achieve convergence
cph_def_exp2_g0_3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq5 + Principal + Instalment + Balance,
                            data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq5: coef = 9.837667177790; exp(coef) = 18725.980642775663; se(coef) = 0.178817524689; z = 55.015; p < 0.0000000000000002
# Principal : coef = 0.000000104316; exp(coef) = 1.000000104316; se(coef) = 0.000000003311; z = 6.545; p < 0.0000000000000002
# Instalment: coef =  0.000001766905; exp(coef) = 1.000001766906; se(coef) = 0.000000107777; z = 16.394; p < 0.0000000000000002
# Balance   : coef = -0.000000049523; exp(coef) =  0.999999950477; se(coef) = 0.000000007928; z = 25.494; p = 0.00000000042

# - Fitting a Cox model with the binned delinquency variable and four additional covariates in an attempt to achieve convergence
cph_def_exp2_g0_4 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                            g0_Delinq5 + Principal + Instalment + Balance + InterestRate_Margin,
                            data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq5          : coef = 6.123756168489; exp(coef) = 456.576455751139; se(coef) = 0.026363059684; z = 232.285; p < 0.0000000000000002
# Principal           : coef = 0.000000072884; exp(coef) = 1.000000072884; se(coef) = 0.000000001831; z = 39.810; p < 0.0000000000000002
# Instalment          : coef =  0.000001302194; exp(coef) = 1.000001302195; se(coef) = 0.000000021809; z = 59.708; p < 0.0000000000000002
# Balance             : coef = 0.000000152576; exp(coef) =  1.000000152576; se(coef) = 0.000000002275; z = 67.076; p < 0.0000000000000002
# InterestRate_Margin : coef = -0.342860721455; exp(coef) =  0.709737055821; se(coef) = 0.258871766074; z = -1.324; p = 0.185

# - Fitting a Cox model with the binned delinquency variable and four additional covariates in an attempt to achieve convergence
cph_def_exp2_g0_5 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5 + Principal + Instalment + Balance + InterestRate_Nom,
                           data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5          : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal           : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Instalment          : coef =  NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Balance             : coef = NA; exp(coef) =  NA; se(coef) = NA; z = NA; p = NA
# InterestRate_Nom    : coef = NA; exp(coef) =  NA; se(coef) = NA; z = NA; p = NA


# --- Fitting the same models as above, without the binned delinquency variable to confirm that non-convergence is caused by it
# - Fitting the model with the largest input space from above
cph_def_exp2_g0_6 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           Principal + Instalment + Balance + InterestRate_Margin,
                           data=dat_train, id=PerfSpell_Key)
# Principal           : coef = -0.00000217116; exp(coef) = 0.99999782884; se(coef) = 0.00000003097; z = -70.105; p < 0.0000000000000002
# Instalment          : coef =  0.00000063399; exp(coef) = 1.00000063399; se(coef) = 0.00000014260; z = 4.446; p = 0.00000875
# Balance             : coef = 0.00000161609; exp(coef) =  1.00000161609; se(coef) = 0.00000003192; z = 50.636; p < 0.0000000000000002
# InterestRate_Margin : coef = 5.67114559608; exp(coef) =  290.36698721056; se(coef) = 0.28579364069; z = 19.843; p < 0.0000000000000002

cph_def_exp2_g0_7 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           Principal + Instalment + Balance,
                           data=dat_train, id=PerfSpell_Key)
# Principal  : coef = -0.00000229664; exp(coef) = 0.99999770336; se(coef) = 0.000000001831; z = -74.844; p < 0.0000000000000002
# Instalment : coef =  0.00000065919; exp(coef) = 1.00000065919; se(coef) = 0.000000021809; z = 4.964; p = 0.00000069
# Balance    : coef = 0.00000173957; exp(coef) =  1.00000173957; se(coef) = 0.000000002275; z = 55.005; p < 0.0000000000000002

cph_def_exp2_g0_8 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           Principal + Instalment,
                           data=dat_train, id=PerfSpell_Key)
# Principal  : coef = -0.000000758056; exp(coef) = 0.999999241944; se(coef) = 0.000000009095; z = -83.35; p < 0.0000000000000002
# Instalment : coef =  0.000000867540; exp(coef) = 1.000000867540; se(coef) = 0.000000072233; z = 12.01; p < 0.0000000000000002

cph_def_exp2_g0_9 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           Principal,
                           data=dat_train, id=PerfSpell_Key)
# Principal : coef = -0.000000749695; exp(coef) = 0.999999250306; se(coef) = 0.000000009029; z = -83.03; p < 0.0000000000000002




# ------ 3. Cox PH Models - Experimentation 3
# --- Combinig g0_Delinq8a (factorised/ dummy encoded variable) with other covariates in attempting to achieve convergence
#     Because of the non-linear hazards expected by different levels of g0_Delinq, fitting a model with a facotrised variable
#     will solve the non-linearity and might aid convergence.
# - Fitting a Cox model with the factorised variable and an additional covariate in an attempt to achieve convergence
cph_def_exp3_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq8b + Principal,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8b_1: coef = -0.002260622394; exp(coef) = 0.997741930889; se(coef) = 286.060062719290; z = 0.000; p = 1.000
# g0_Delinq8b_2: coef = 24.522033930499; exp(coef) = 44646150721.102226257324; se(coef) = 65.053053796357; z = 0.377; p = 0.706
# Principal    : coef = 0.000000066238; exp(coef) = 1.000000066238; se(coef) = 0.000000006491; z = 10.205; p < 0.0000000000000002

# - Fitting a Cox model with the factorised variable and two additional covariates in an attempt to achieve convergence
cph_def_exp3_g0_2 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq8b + Principal + Instalment,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8b_1: coef = 55.89130632944; exp(coef) = 1876229349380512165220206.00; se(coef) = 1470.34007999316; z = 0.038; p = 0.970
# g0_Delinq8b_2: coef = 83.74133803887; exp(coef) = 2335613718562459164260488840884460448.00; se(coef) = 48.92421052491; z = 1.712; p = 0.087
# Principal    : coef = 0.00000036971; exp(coef) = 1.00000036971; se(coef) = 0.00000001528; z = 24.198; p < 0.0000000000000002
# Instalment   : coef =  -0.00003365831; exp(coef) = 0.99996634226 ; se(coef) = 0.00000163210; z = -20.623; p < 0.0000000000000002

# - Fitting a Cox model with the factorised variable and three additional covariates in an attempt to achieve convergence
cph_def_exp3_g0_3 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq8b + Principal + Instalment + Balance,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8b_1: coef = 46.98688069854; exp(coef) = 254766903872029917184.00; se(coef) = 895.87602633580; z = 0.052; p = 0.958
# g0_Delinq8b_2: coef = 73.84152873669; exp(coef) = 117211025489903045658844622086226.00; se(coef) = 54.17291319420; z = 1.363; p = 0.173
# Principal    : coef = 0.00000032113; exp(coef) = 1.00000032113; se(coef) = 0.00000002022; z = 15.881; p < 0.0000000000000002
# Instalment   : coef =  -0.00003776544; exp(coef) = 0.99996223527 ; se(coef) = 0.00000192801; z = -19.588; p < 0.0000000000000002
# Balance      : coef =  0.00000009695; exp(coef) = 1.00000009695 ; se(coef) = 0.00000002491; z = 3.892; p = 0.0000995

# - Fitting a Cox model with the factorised variable and three additional covariates in an attempt to achieve convergence
cph_def_exp3_g0_4 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           g0_Delinq8b + Principal + Instalment + Balance + InterestRate_Margin,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# g0_Delinq8b_1      : coef = 45.59171945973; exp(coef) = 63129475610463158272.00; se(coef) = 1476.99139095606; z = 0.031; p = 0.975
# g0_Delinq8b_2      : coef = 73.44675133997; exp(coef) = 78980307571572965464662460284648.00; se(coef) = 81.20233625496; z = 0.904; p = 0.366
# Principal          : coef = 0.00000031532; exp(coef) = 1.00000031532; se(coef) = 0.00000002067; z = 15.256; p < 0.0000000000000002
# Instalment         : coef =  -0.00003759308; exp(coef) = 0.99996240762 ; se(coef) = 0.00000193563; z = -19.422; p < 0.0000000000000002
# Balance            : coef =  0.00000010038; exp(coef) = 1.00000010038 ; se(coef) = 0.00000002509; z = 4.001; p = 0.000063
# InterestRate_Margin: coef =  -0.39966782697; exp(coef) = 0.67054274526 ; se(coef) = 0.26154427388; z = -1.528; p = 0.126




# ------ 4. Cox PH Models - Experimentation 4
### Time in delinquency state
# --- Raw time variable
cph_def_exp4_g0_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           TimeInState,
                           data=dat_train, id=PerfSpell_Key)
# TimeInState: coef = -0.490850; exp(coef) = 0.612106; se(coef) = 0.001995; z = -246; p < 0.0000000000000002


# --- State Direction
# - Numerical variable
cph_def_exp4_g0_2_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             StateDir,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# StateDir: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - Categorical variable
cph_def_exp4_g0_2_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              StateDir2,
                              data=dat_train, id=PerfSpell_Key)
# TimeInState: coef = -0.490850; exp(coef) = 0.612106; se(coef) = 0.001995; z = -246; p < 0.0000000000000002


# --- Raw time variable with the direction of the last transition in the delinquency state
# - Numerical state direction variable
cph_def_exp4_g0_3_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             TimeInState + StateDir,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# TimeInState: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateDir: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - Categorical state direction variable
cph_def_exp4_g0_3_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              TimeInState + StateDir2,
                              data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# TimeInState: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateDir2: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - Raw time variable with an interaction with the direction of the last transition in the delinquency state
cph_def_exp4_g0_4 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           TimeInState*StateDir,
                           data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# TimeInState*StateDir: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

cph_def_exp4_g0_5 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           TimeInState + g0_Delinq5,
                           data=dat_train, id=PerfSpell_Key)
# WARNING:  Ran out of iterations and did not converge.
# TimeInState: coef = 0.092421; exp(coef) = 1.096827; se(coef) = 79.44; z = -246; p < 0.0000000000000002
# g0_Delinq5 : coef = 24.847552; exp(coef) = 61823654252.605461; se(coef) = 191.441681; z = 0.13; p = 0.897
# TimeInState : coef = 24.847552; exp(coef) = 61823654252.605461; se(coef) = 191.441681; z = 0.13; p = 0.897


# --- Raw time variable with state number (counter)
cph_def_exp4_g0_6 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           TimeInState + State_Num,
                           data=dat_train, id=PerfSpell_Key)
# TimeInState: coef = -0.4850934; exp(coef) = 0.6156397; se(coef) = 0.0020025; z = -242.25; p < 0.0000000000000002
# g0_Delinq5 : coef = 0.0144557; exp(coef) = 1.0145607; se(coef) = 0.0005299; z = 27.28; p = 0.897


# --- Raw time variable with state volatility
# - Spell level volatility
cph_def_exp4_g0_7_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             TimeInState + StateVol,
                             data=dat_train, id=PerfSpell_Key)
# TimeInState: coef = -0.533032; exp(coef) = 0.586823; se(coef) = 0.002055; z = -259.4; p < 0.0000000000000002
# StateVol   : coef = -0.224422; exp(coef) = 0.798978; se(coef) = 0.001881; z = -119.3; p < 0.0000000000000002

# - 3-Month rolling volatility
cph_def_exp4_g0_7_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              TimeInState + StateVol_3_Months,
                              data=dat_train, id=PerfSpell_Key)
# TimeInState       : coef = -0.390937; exp(coef) = 0.676423; se(coef) = 0.002269; z = -172.30; p < 0.0000000000000002
# StateVol_3_Months : coef = 0.855842; exp(coef) = 2.353355; se(coef) = 0.012246; z = 69.89; p < 0.0000000000000002

# - 6-Month rolling volatility
cph_def_exp4_g0_7_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               TimeInState + StateVol_6_Months,
                               data=dat_train, id=PerfSpell_Key)
# TimeInState       : coef = -0.349910; exp(coef) = 0.704752; se(coef) = 0.002098; z = -166.8; p < 0.0000000000000002
# StateVol_6_Months : coef = 1.125403; exp(coef) = 3.081457; se(coef) = 0.008651; z = 130.1; p < 0.0000000000000002


# --- Raw time variable with state direction and state volatility
# - Spell level volatility
cph_def_exp4_g0_8 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                           TimeInState + StateDir + StateVol,
                           data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# TimeInState: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateDir: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol: coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA




# ------ 5. Cox PH Models - Experimentation 5
### Binned delinquency variable lagged by different number of months
# --- 1-Month lag (Numerical)
### NOTE:~ 457 6245 observations deleted due to missingness ~ 457624/29071039 = 1.574158% of observations
# - 1-Month lag
cph_def_exp5_g0_1_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_1,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState
cph_def_exp5_g0_1_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_1 + TimeInState,
                              data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState and StateVol
cph_def_exp5_g0_1_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1 + TimeInState + StateVol,
                               data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState, StateVol, and PerfSpell_Num
cph_def_exp5_g0_1_iv <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_1 + TimeInState + StateVol + PerfSpell_Num,
                              data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# PerfSpell_Num    : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState, StateVol, PerfSpell_Num and Principal
cph_def_exp5_g0_1_v <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_1 + TimeInState + StateVol + PerfSpell_Num + Principal,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# PerfSpell_Num    : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal        : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag Principal
cph_def_exp5_g0_1_vi <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_1 + Principal,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal        : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState and Principal
cph_def_exp5_g0_1_vii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1 + TimeInState + Principal,
                               data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_1 : coef = 3.646152248864; exp(coef) = 38.326909556434; se(coef) = 0.008014793261; z = 454.93; p < 0.0000000000000002
# TimeInState      : coef = -0.116441213637; exp(coef) = 0.890082420147; se(coef) = 0.001821092475; z = -63.94; p < 0.0000000000000002
# Principal        : coef = 0.000000053357; exp(coef) = 1.000000053357; se(coef) = 0.000000002135; z = 24.99; p < 0.0000000000000002

# - 1-Month lag with TimeInState, StateVol, and Principal
cph_def_exp5_g0_1_viii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                g0_Delinq5_lag_1 + TimeInState + Principal + Instalment,
                                data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal        : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# --- 1-Month lag (categorical)
# - 1-Month lag
cph_def_exp5_g0_2_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_1_f,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_f2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_fMISSING : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag and TimeInState
cph_def_exp5_g0_2_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_1_f + TimeInState,
                              data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_f2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_fMISSING : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState      : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState and Principal
cph_def_exp5_g0_2_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1_f + TimeInState + Principal,
                               data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_1_f1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_f2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_fMISSING : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState               : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal                 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month lag with TimeInState, StateVol, and Principal
cph_def_exp5_g0_2_iv <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_1_f + TimeInState + StateVol + Principal,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_1_f1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_f2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_1_fMISSING : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# TimeInState               : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# StateVol                  : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# Principal                 : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# --- 1-Month lag (numerical) and binned delinquency
cph_def_exp5_g0_3_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_1 + g0_Delinq5,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5         : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA


# --- 2-Month lag (numerical)
### NOTE:~ 913 375 observations deleted due to missingness ~ 913375/29071039 = 3.141873% of observations

# - 2-Month lag
cph_def_exp5_g0_4_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_2,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2 : coef = -0.2602589; exp(coef) = 0.7708519; se(coef) = 0.0007645; z = -286; p < 0.0000000000000002

# - 2-Month lag with TimeInState
cph_def_exp5_g0_4_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_2 + TimeInState,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2 : coef = 1.819207; exp(coef) = 6.166963; se(coef) = 0.005427; z = 335.2; p < 0.0000000000000002
# TimeInState      : coef = -0.359717; exp(coef) = 0.697874; se(coef) = 0.002074; z = -173.5; p < 0.0000000000000002

# - 2-Month lag with Principal
cph_def_exp5_g0_4_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_2 + Principal,
                               data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_2 : coef = 2.920419631755; exp(coef) = 18.549069605440; se(coef) = 0.004294035234; z = 680.111; p < 0.0000000000000002
# Principal        : coef = -0.000000051973; exp(coef) = 0.999999948027; se(coef) = 0.000000007587; z = -6.851; p = 0.00000000000735

# - 2-Month lag with TimeInState and Principal
cph_def_exp5_g0_4_iv <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_2 + TimeInState + Principal,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2 : coef = 1.820581267907; exp(coef) = 6.175446995980; se(coef) = 0.005427165639; z = 335.46; p < 0.0000000000000002
# TimeInState      : coef = -0.360237191132; exp(coef) = 0.697510863058; se(coef) = 0.002075670219; z = -173.55; p < 0.0000000000000002
# Principal        : coef = 0.000000049881; exp(coef) = 1.000000049881; se(coef) = 0.000000002846; z = 17.52; p < 0.0000000000000002

# - 2-Month lag with significantly increased input space
cph_def_exp5_g0_4_v <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_2 + TimeInState + Principal + Instalment + Balance + PerfSpell_Num +
                             InterestRate_Margin + M_Repo_Rate + slc_pmnt_method + LN_TPE,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2    : coef = 1.675767427583; exp(coef) = 5.342893860261; se(coef) = 0.005608726197; z = 298.779; p < 0.0000000000000002
# TimeInState         : coef = -0.359349954488; exp(coef) = 0.698129994872; se(coef) = 0.002086045718; z = -172.264; p < 0.0000000000000002
# Principal           : coef = 0.000000050064; exp(coef) = 1.000000050064; se(coef) = 0.000000003231; z = 15.493; p < 0.0000000000000002
# Instalment          : coef = -0.000045458345; exp(coef) = 0.999954542688; se(coef) = 0.000002120740; z = -21.435; p < 0.0000000000000002
# Balance             : coef = 0.000000463613; exp(coef) = 1.000000463613; se(coef) = 0.000000021364; z = 21.701; p < 0.0000000000000002
# PerfSpell_Num       : coef = -0.013283091987; exp(coef) = 0.986804738961; se(coef) = 0.005884816951; z = -2.257; p = 0.0240
# InterestRate_Margin : coef = 3.831339725528; exp(coef) = 46.124290749059; se(coef) = 0.235063885110; z = 16.299; p < 0.0000000000000002
# M_Repo_Rate         : coef = -0.332205613990; exp(coef) = 0.717339812589; se(coef) = 0.168187610901; z = -1.975; p = 0.0482
# slc_pmnt_methodDebit Order Other Bank : coef = 0.021322129550; exp(coef) = 1.021551070428; se(coef) = 0.014471726896; z = 1.473; p = 0.1407
# slc_pmnt_methodDebit MISSING_DATA     : coef = 0.056730007105; exp(coef) = 1.058370019431; se(coef) = 0.012438387509; z = 4.561; p = 0.00000509
# slc_pmnt_methodDebit Salary           : coef = 0.013143498502; exp(coef) = 1.013230253952; se(coef) = 0.018766669762; z = 0.700; p = 0.4837
# slc_pmnt_methodDebit Statement        : coef = 0.310227371309; exp(coef) = 1.363735153130; se(coef) = 0.011240589962; z = 27.599; p < 0.0000000000000002
# slc_pmnt_methodDebit Suspense         : coef = 1.987697172039; exp(coef) = 7.298706728930; se(coef) = 0.013879375163; z = 143.212; p < 0.0000000000000002
# LN_TPE WHL                            : coef = 0.202480634959; exp(coef) = 1.224436373632; se(coef) = 0.016571442554; z = 12.219; p < 0.0000000000000002


# --- 2-Month lag (categorical)
# - 2-Month lag
cph_def_exp5_g0_5_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_2_f,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2_f1       : coef = 5.13299; exp(coef) = 169.52349; se(coef) = 0.01343; z = 382.12; p < 0.0000000000000002
# g0_Delinq5_lag_2_f2       : coef = 6.52259; exp(coef) = 680.34011; se(coef) = 0.01488; z = 438.26; p < 0.0000000000000002
# g0_Delinq5_lag_2_fMISSING : coef = 2.19768; exp(coef) = 9.00412; se(coef) = 0.03628; z = 60.57; p < 0.0000000000000002

# - 2-Month lag with TimeInState
cph_def_exp5_g0_5_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_2_f + TimeInState,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2_f1       : coef = 3.067218; exp(coef) = 21.482049; se(coef) = 0.014179; z = 216.325; p < 0.0000000000000002
# g0_Delinq5_lag_2_f2       : coef = 4.231954; exp(coef) = 68.851642; se(coef) = 0.015636; z = 270.654; p < 0.0000000000000002
# g0_Delinq5_lag_2_fMISSING : coef = -0.267590; exp(coef) = 0.765221; se(coef) = 0.036677; z = -7.296; p < 0.0000000000000002
# TimeInState               : coef = -0.307532; exp(coef) = 0.735260; se(coef) = 0.002069; z = -148.654; p < 0.0000000000000002

# - 2-Month lag with Principal
cph_def_exp5_g0_5_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_2_f + Principal,
                               data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2_f1       : coef = 5.142027598737; exp(coef) = 171.062262554406; se(coef) = 0.013451599754; z = 382.26; p < 0.0000000000000002
# g0_Delinq5_lag_2_f2       : coef = 6.531187320775; exp(coef) = 686.212482332343; se(coef) = 0.014898341555; z = 438.38; p < 0.0000000000000002
# g0_Delinq5_lag_2_fMISSING : coef = 2.210558119060; exp(coef) = 9.120805468131; se(coef) = 0.036293301636; z = 60.91; p < 0.0000000000000002
# Principal                 : coef = 0.000000045189; exp(coef) = 1.000000045189; se(coef) = 0.000000002402; z = 18.81; p < 0.0000000000000002

# - 2-Month lag with TimeInState and Principal
cph_def_exp5_g0_5_iv <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_2_f + TimeInState + Principal,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2_f1       : coef = 3.072494463436; exp(coef) = 21.595705239930; se(coef) = 0.014195439951; z = 216.442; p < 0.0000000000000002
# g0_Delinq5_lag_2_f2       : coef = 4.237216086427; exp(coef) = 69.214895090047; se(coef) = 0.015650761432; z = 270.735; p < 0.0000000000000002
# g0_Delinq5_lag_2_fMISSING : coef = -0.257227458239; exp(coef) = 0.773192324798; se(coef) = 0.036685814575; z = -7.012; p < 0.0000000000000002
# TimeInState               : coef = -0.307707911825; exp(coef) = 0.735130009436; se(coef) = 0.002070247878; z = -148.633; p < 0.0000000000000002
# Principal                 : coef = 0.000000049669; exp(coef) = 1.000000049669; se(coef) = 0.000000003145; z = 15.792; p < 0.0000000000000002

# - 2-Month lag with significantly increased input space
cph_def_exp5_g0_5_v <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_2_f + TimeInState + Principal + Instalment + Balance + PerfSpell_Num +
                             InterestRate_Margin + M_Repo_Rate + slc_pmnt_method + LN_TPE,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2_f1                   : coef = 6.743713954117; exp(coef) = 848.706948455372; se(coef) = 0.038300175676; z = 176.075; p < 0.0000000000000002
# g0_Delinq5_lag_2_f2                   : coef = 6.757360047683; exp(coef) = 860.367864908867; se(coef) = 0.040552939945; z = 166.631; p < 0.0000000000000002
# g0_Delinq5_lag_2_fMISSING             : coef = 2.298728778009; exp(coef) = 9.961511110490; se(coef) = 0.080917634117; z = 28.408; p < 0.0000000000000002
# TimeInState         : coef = -0.022787288839; exp(coef) = 0.977470380521; se(coef) = 0.000765830857; z = -29.755; p < 0.0000000000000002
# Principal           : coef = 0.000000050064; exp(coef) = 1.000000050064; se(coef) = 0.000000003231; z = 256.879; p < 0.0000000000000002
# Instalment          : coef = 0.000001315258 exp(coef) = __________; se(coef) = 0.000000006178; z = 212.885; p < 0.0000000000000002
# Balance             : coef = -0.000000012798; exp(coef) = 0.999999987202; se(coef) = 0.000000003644; z = -3.513; p = 0.000444
# PerfSpell_Num       : coef = -0.138866005962; exp(coef) = 0.870344641638; se(coef) = 0.005873059552; z = -23.645; p < 0.0000000000000002
# InterestRate_Margin : coef = 1.381035578159; exp(coef) = 3.979020080926; se(coef) = 0.229689979501; z = 6.013; p = 0.00000000183
# M_Repo_Rate         : coef = -7.401417596628; exp(coef) = 0.000610386865; se(coef) = 0.167204349584; z = -44.266; p < 0.0000000000000002
# slc_pmnt_methodDebit Order Other Bank : coef = -0.035055617976; exp(coef) = 0.965551712732; se(coef) = 0.014120618562; z = -2.483; p = 0.013043
# slc_pmnt_methodDebit MISSING_DATA     : coef = 0.050876223534; exp(coef) = 1.052192648535; se(coef) = 0.012354470954; z = 4.118; p = 0.00003821060
# slc_pmnt_methodDebit Salary           : coef = -0.362598510313; exp(coef) = 0.695865760348; se(coef) = 0.019586984645; z = -18.512; p = 0.4837
# slc_pmnt_methodDebit Statement        : coef = 0.459243002451; exp(coef) = 1.582875298629; se(coef) = 0.010684964633; z = 42.980; p < 0.0000000000000002
# slc_pmnt_methodDebit Suspense         : coef = 1.203797971692; exp(coef) = 3.332750608828; se(coef) = 0.018230648912; z = 66.032; p < 0.0000000000000002
# LN_TPE WHL                            : coef = 0.165624077402; exp(coef) = 1.180129381023; se(coef) = 0.015522927120; z =  10.670; p < 0.0000000000000002


# --- 2-Month lag (numerical) and binned delinquency
# - 2-Month lag with binned delinquency
cph_def_exp5_g0_6_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_2 + g0_Delinq5,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_2 : coef = 0.074905; exp(coef) = 1.077782; se(coef) = 0.006637; z = 11.29; p < 0.0000000000000002
# g0_Delinq5       : coef = 24.395131; exp(coef) = 39325173003.143692; se(coef) = 270.290858; z = 0.09; p = 0.928


# --- 3-Month lag (numerical)
### NOTE:~ 1 361 499 observations deleted due to missingness ~ 1361499/29071039 = 4.683352% of observations

# - 3-Month lag
cph_def_exp5_g0_7_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_3,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_3 : coef = -0.2312096; exp(coef) = 0.7935731; se(coef) = 0.0007966; z = -252.9 p < 0.0000000000000002

# - 3-Month lag with TimeInState
cph_def_exp5_g0_7_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_3 + TimeInState,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_3 : coef = 1.412604; exp(coef) = 4.106636; se(coef) = 0.005467; z = 335.2; p < 0.0000000000000002
# TimeInState      : coef = -0.429079; exp(coef) = 0.651108; se(coef) = 0.002094; z = -173.5; p < 0.0000000000000002

# - 3-Month lag with Principal
cph_def_exp5_g0_7_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_3 + Principal,
                               data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_3 : coef = 2.698286484475; exp(coef) = 14.854256905532; se(coef) = 0.004355018710; z = 619.58; p < 0.0000000000000002
# Principal        : coef = -0.000000111776; exp(coef) = 0.999999888224; se(coef) = 0.000000007935; z = -14.09; p = 0.00000000000735

# - 3-Month lag with TimeInState and Principal
cph_def_exp5_g0_7_iv <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_3 + TimeInState + Principal,
                              data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_3 : coef = 1.414419344683; exp(coef) = 4.114096900044; se(coef) = 0.005466750633; z = 258.73; p < 0.0000000000000002
# TimeInState      : coef = -0.429709407275; exp(coef) = 0.650698155402; se(coef) = 0.002096250088; z = -204.99; p < 0.0000000000000002
# Principal        : coef = 0.000000051743; exp(coef) = 1.000000051743; se(coef) = 0.000000002754; z = 18.79; p < 0.0000000000000002

# - 3-Month lag with significantly increased input space
cph_def_exp5_g0_7_v <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_3 + TimeInState + Principal + Instalment + Balance + PerfSpell_Num +
                             InterestRate_Margin + M_Repo_Rate + slc_pmnt_method + LN_TPE,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_2    : coef = 1.282784486386; exp(coef) = 3.606668475901; se(coef) = 0.005616000930; z = 228.416; p < 0.0000000000000002
# TimeInState         : coef = -0.417514093389; exp(coef) = 0.658682208718; se(coef) = 0.002103070187; z = -198.526; p < 0.0000000000000002
# Principal           : coef = 0.000000051193; exp(coef) = 1.000000051193; se(coef) = 0.002103070187; z = 15.745; p < 0.0000000000000002
# Instalment          : coef = -0.000045595922; exp(coef) = 0.999954405118; se(coef) = 0.000002011418; z = -22.669; p < 0.0000000000000002
# Balance             : coef = 0.000000458003; exp(coef) = 1.000000458003; se(coef) = 0.000000020128; z = 22.755; p < 0.0000000000000002
# PerfSpell_Num       : coef = -0.032517009623; exp(coef) = 0.968005984274; se(coef) = 0.006145446079; z = -5.291; p = 0.000000121
# InterestRate_Margin : coef = 3.277916832407; exp(coef) = 26.520468535358; se(coef) = 0.241851062956; z = 13.553; p < 0.0000000000000002
# M_Repo_Rate         : coef = -0.073248235426; exp(coef) = 0.929370098807; se(coef) = 0.169566215342; z = -0.432; p = 0.6658
# slc_pmnt_methodDebit Order Other Bank : coef = 0.037264233465; exp(coef) = 1.037967250290; se(coef) = 0.014552715642; z = 2.561; p = 0.0104
# slc_pmnt_methodDebit MISSING_DATA     : coef = 0.108174979284 ; exp(coef) = 1.114242697719; se(coef) = 0.012636810664; z = 8.560; p < 0.0000000000000002
# slc_pmnt_methodDebit Salary           : coef = -0.020225353046; exp(coef) = 0.979977807437; se(coef) = 0.018871947674; z = -1.072; p = 0.2838
# slc_pmnt_methodDebit Statement        : coef = 0.340186916408; exp(coef) = 1.405210222865; se(coef) = 0.011342802659; z = 29.991; p < 0.0000000000000002
# slc_pmnt_methodDebit Suspense         : coef = 2.274527673792; exp(coef) = 9.723325346402; se(coef) = 0.013933257739; z = 163.244; p < 0.0000000000000002
# LN_TPE WHL                            : coef = 0.208531220813; exp(coef) = 1.231867389362; se(coef) = 0.016774333099; z = 12.432; p < 0.0000000000000002


# --- 3-Month lag (categorical)
# - 3-Month lag
cph_def_exp5_g0_8_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_3_f,
                             data=dat_train, id=PerfSpell_Key)
# g0_Delinq5_lag_3_f1       : coef = 4.299900; exp(coef) = 73.692413; se(coef) = 0.009408; z = 457.02; p < 0.0000000000000002
# g0_Delinq5_lag_3_f2       : coef = 5.194348; exp(coef) = 180.250497; se(coef) = 0.013339; z = 389.42; p < 0.0000000000000002
# g0_Delinq5_lag_3_fMISSING : coef = 1.444620; exp(coef) = 4.240239; se(coef) = 0.027687; z = 52.18; p < 0.0000000000000002


# --- 3-Month lag (numerical) and binned delinquency
# - 3-Month lag with binned delinquency
cph_def_exp5_g0_9_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                              g0_Delinq5_lag_3 + g0_Delinq5,
                              data=dat_train, id=PerfSpell_Key)
# WARNING:  Loglik converged before variable  2 ; beta may be infinite.
# g0_Delinq5_lag_3 : coef =  0.208827; exp(coef) = 1.232231; se(coef) = 0.005592; z = 37.343; p < 0.0000000000000002
# g0_Delinq5       : coef = 24.653248; exp(coef) = 50906070735.790077; se(coef) = 315.113143; z = 0.078; p = 0.938


# --- Combining lags
# - 1-Month and 2-month lag (numerical)
cph_def_exp5_g0_10_i <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                             g0_Delinq5_lag_1 + g0_Delinq5_lag_2,
                             data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month and 3-month lag (numerical)
cph_def_exp5_g0_10_ii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1 + g0_Delinq5_lag_3,
                               data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_3       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month, 2-month lag, and 3-month lag (numerical)
cph_def_exp5_g0_10_iii <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                                g0_Delinq5_lag_1 + g0_Delinq5_lag_2 + g0_Delinq5_lag3,
                                data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA
# g0_Delinq5_lag_2       : coef = NA; exp(coef) = NA; se(coef) = NA; z = NA; p = NA

# - 1-Month and 2-month lag (categorical)
cph_def_exp5_g0_10_iv <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1)~
                               g0_Delinq5_lag_1_f + g0_Delinq5_lag_2_f,
                               data=dat_train, id=PerfSpell_Key)
# Error in agreg.fit(X, Y, istrat, offset, init, control, weights = weights,  : exp overflow due to covariates
# g0_Delinq5_lag_1_f       : coef = 4.299900; exp(coef) = 73.692413; se(coef) = 0.009408; z = 457.02; p < 0.0000000000000002
# g0_Delinq5_lag_2_f       : coef = 5.194348; exp(coef) = 180.250497; se(coef) = 0.013339; z = 389.42; p < 0.0000000000000002



