# ======================================= INPUT SPACE: Cox PH (PWP) BASIC ==============================
# Divide data into thematic groups and perform data analysis on them to compile an input space for a
# basic Cox Proportional Hazards (CPH) model, having used the PWPST-time defintiion
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





# ------ 1. Final model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# - remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Initialize variables
# Selection based on expert judgement alone, where the initial list being based on thematic selection
vars2 <- c("Arrears", "InterestRate_Nom")

# - Build model based on variables
cox_PWPST_basic <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", paste(vars2,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST_basic); AIC(cox_PWPST_basic); concordance(cox_PWPST_basic)
### RESULTS: AIC: 221736.4 Harrel's C: 0.891   

c <- coefficients(cox_PWPST_basic)
(c <- data.table(Variable=names(c),Coefficient=c))
#1:          Arrears  0.00005137344
#2: InterestRate_Nom 14.87640164574

# - Test Goodness of fit using bootstrapped B-statistics (1-KS statistic) over single-factor models
csTable_PWPST <- csTable(datCredit_train, vars2, TimeDef="PWPST", strataVar="PerfSpell_Num_binned", seedVal=1, numIt=10,
                         fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="DefaultStatus1", genPath=genPath)
### RESULTS: Top single-factor models: Arrears, InterestRate_Nom
# Results do not vary much from each other, at least not meaningfully

# - Test goodness-of-fit using AIC-measure
aicTable_PWPST <- aicTable(datCredit_train, variables=vars2, fldSpellID="PerfSpell_Key",
                           TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genPath)
### RESULTS: Top single-factor models: Arrears, InterestRate_Nom
# Where the first 3 results have AIC values significantly different from the rest.

# Test accuracy using Harrell's c-statistic over single-factor models
concTable_PWPST <- concTable(datCredit_train, datCredit_valid, variables=vars2, 
                             fldSpellID="PerfSpell_Key", TimeDef="PWPST", strataVar="PerfSpell_Num_binned", genPath=genPath)
### RESULTS: Top single-factor models (>85%): Arrears, InterestRate_Nom

# - Combine results into a single object
Table_PWPST_basic <- concTable_PWPST[,1:2] %>% left_join(aicTable_PWPST, by ="Variable") %>% 
  left_join(data.table(csTable_PWPST$Results), by="Variable")

# - Test Goodnes-of-fit using Cox-Snell, having measured distance between residual distribution and unit exponential using KS-statistic
GoF_CoxSnell_KS(cox_PWPST_basic, datCredit_train, GraphInd=TRUE, legPos=c(0.6,0.4), panelTitle="Prentice-Williams-Peterson (PWP) model: Basic",
                fileName = paste0(genFigPath, "KS_Test_CoxSnellResiduals_Exp_PWPST_basic", ".png"), dpi=280) # 0.6167
### RESULTS: Goodness of fit for the model seems to be a bit low.



# - Manual test of Cox-Snell residuals
# 1. Get baseline cumulative hazard
datBaseHaz <- as.data.table(basehaz(cox_PWPST_basic, centered=F))
# 2. Compute linear predictor
datCredit_train[, LP := predict(cox_PWPST_basic, type="lp")]
# 3a. Iterate per stratum
for (i in 1:length(unique(datBaseHaz$strata))) {
  # i <- 1
  # 3b. Approximate H_0(T_i) for each subject using interpolation since "basehaz" returns discrete values
  getH0 <- approxfun(datBaseHaz[strata==paste0("PerfSpell_Num_binned=",i), time], 
                     datBaseHaz[strata==paste0("PerfSpell_Num_binned=",i), hazard],
                     rule = 2)
  datCredit_train[PerfSpell_Num_binned==i, H0 := getH0(TimeInPerfSpell)]
}
# 4. Compute Cox-Snell (CS) residuals as per definition H_0(t).exp(X^T \beta )
datCredit_train[, CS_residual := H0 * exp(LP)]
# 5. Treat CS-residuals as survival times towards computing cumulative hazard using Kaplan-Meier estimate of S(t)
r_KM <- survfit(Surv(CS_residual, DefaultStatus1) ~ 1, data=datCredit_train)
# 6. Plot diagnostic
datPlot <- data.table(time=r_KM$time, surv=r_KM$surv); sLimit <- 1000
plot(datPlot[1:sLimit,time], -log(datPlot[1:sLimit, surv]), type="s")
abline(0,1, col="red", lty=2)
### RESULTS: Deviates significantly from the red-line.



# - Save objects
pack.ffdf(paste0(genObjPath,"PWPST_Univariate_Models_basic"), Table_PWPST_basic)
pack.ffdf(paste0(genPath,"PWPST_Cox_Model_basic"), cox_PWPST_basic)

# - Cleanup
rm(datCredit_train, datCredit_valid,cox_PWPST_basic, aicTable_PWPST, concTable_PWPST, csTable_PWPST, c); gc()
