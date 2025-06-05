# ======================================= INPUT SPACE: Cox PO (PWP) ====================================
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





# ------ 1. Preliminaries

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# - Create start point variable
datCredit_train[, Start := TimeInPerfSpell - 1]
datCredit_valid[, Start := TimeInPerfSpell - 1]

# - remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Bespoke model evaluation function
evalLR <- function(model, model_base, datGiven, targetFld, predClass) {
  require(data.table); require(scales)
  # - Test conditions
  # model <- modLR; model_base <- modLR_base; datGiven <- datCredit_train
  # targetFld = "PerfSpell_Event"; predClass <- 1
  result1 <- AIC(model) # 1164537 
  result2 <- coefDeter_glm(model, model_base) # 0.29%
  matPred <- predict(model, newdata=datGiven, type="response")
  actuals <- ifelse(datGiven[[targetFld]] == predClass, 1,0)
  result3 <- roc(response=actuals, predictor = matPred)
  objResults <- data.table(AIC=comma(result1), result2, AUC=percent(result3$auc,accuracy=0.01))
  return(objResults)
  # - Cleanup, if run interactively
  rm(result1, result2, matPred, actuals, result3, objResults, model, model_base, datGiven, targetFld, predClass)
}

# - Fit an empty model as a performance gain, used within some diagnostic functions
modLR_base <- glm(PerfSpell_Event ~ 1, 
                  data=datCredit_train, family="binomial")



# ------ 2a. Modelling theme: Handling time towards embedding baseline hazard h_0(t)
# Which way in embedding time is best?

# --- Single-factor models: Single time variable
modLR <- glm(PerfSpell_Event ~  TimeInPerfSpell,
                   data = datCredit_train, family="binomial")
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  177,477; McFadden R^2:  0.26%; AUC:  59.12%

# --- Single-factor models: Single time variable (transform)
modLR <- glm(PerfSpell_Event ~  log(TimeInPerfSpell),
             data = datCredit_train, family="binomial")
evalLR(modLR, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  176,610; McFadden R^2:  0.75%; AUC:  59.12%
# Log-transformation affords a better goodness-of-fit, though AUC understandably remains unchanged





##### SCRATCH

# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC_PWPST_CoxPO_bas <- tROC.multi(datGiven=datCredit_valid, modGiven=modLR, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                                fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                                graphName="ROC_CoxPO_Basic_TimeVar", genFigPath=paste0("tROC-Analyses/",genFigPath), 
                                caseStudyName=paste0("CoxPO_PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC_PWPST_CoxPO_bas$AUC; objROC_PWPST_CoxPO_bas$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.01%, achieved in 201.11 secs