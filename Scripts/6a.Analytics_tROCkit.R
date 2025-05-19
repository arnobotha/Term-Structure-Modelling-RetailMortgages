# ========================= TIME-DEPENDENT ROC-ANALYSIS ==========================
# Compare various functions from various packages in conducting time-dependent 
# ROC-analyses on the same fitted Cox regression model, having used the 
# prepared credit data
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# --------------------------------------------------------------------------------
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

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b
#
# -- Outputs:
#   - <Analytics> | tROC-graphs
# ================================================================================





# ----------------- 1. Load & prepare data for tROC-analyses

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# Create start and stop columns
datCredit_train[, Start := TimeInPerfSpell - 1]
datCredit_valid[, Start := TimeInPerfSpell - 1]




# ----------------- 2. Fit a Cox regression model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic
# - Initialize variables
vecVars_PWPST_bas <- c("Arrears", "InterestRate_Nom")

# - Build model based on variables
cox_PWPST_basic <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                           paste(vecVars_PWPST_bas,collapse=" + "), 
                                           " + strata(PerfSpell_Num_binned)")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST_basic); AIC(cox_PWPST_basic); concordance(cox_PWPST_basic)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced
# - Initialize variables
vecVars_PWPST_adv <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_adv <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                         paste(vecVars_PWPST_adv,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, datCredit_train, ties="efron")
summary(cox_PWPST_adv); AIC(cox_PWPST_adv); concordance(cox_PWPST_adv)





# ----------------- 3. Conduct & compare ROC-analyses at fixed prediction time
# NOTE: ROC-analyses are conducted using the same fitted Cox model at a specific prediction time t

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic

# --- Package: tROCkit() | custom "package"/function
# NOTE: Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator and 0/1-kernel


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC1_PWPST_bas <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_basic, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                                fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                                graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_bas", genFigPath=paste0(genFigPath), 
                                caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC1_PWPST_bas$AUC; objROC1_PWPST_bas$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.01%, achieved in 201.11 secs


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 12 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 12
objROC2_PWPST_bas <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_basic, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4,
                                fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                                graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_bas", genFigPath=paste0(genFigPath),
                                caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC2_PWPST_bas$AUC; objROC2_PWPST_bas$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 93.05%, achieved in 346.30  secs


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 24
objROC3_PWPST_bas <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_basic, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                                fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                                graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_bas", genFigPath=paste0(genFigPath), 
                                caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC3_PWPST_bas$AUC; objROC3_PWPST_bas$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 93.33%, achieved in 535.03  secs


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 36
objROC4_PWPST_bas <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_basic, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                                fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                                graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_bas", genFigPath=paste0(genFigPath), 
                                caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC4_PWPST_bas$AUC; objROC4_PWPST_bas$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 93.21%, achieved in 730.53  secs




# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced

# --- Package: tROCkit() | custom "package"/function
# NOTE: Using custom tROC()-function from script 0b(iii) under the CD-approach with an NN-estimator and 0/1-kernel


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 3 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 3
objROC1_PWPST_adv <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_adv, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_adv", genFigPath=paste0(genFigPath), 
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC1_PWPST_adv$AUC; objROC1_PWPST_adv$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.19%, achieved in 712.29 secs ( 11.9 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 12 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 12
objROC2_PWPST_adv <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_adv, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4,
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_adv", genFigPath=paste0(genFigPath),
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC2_PWPST_adv$AUC; objROC2_PWPST_adv$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.24%, achieved in 1370.21 secs (22.8 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 24
objROC3_PWPST_adv <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_adv, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_adv", genFigPath=paste0(genFigPath), 
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC3_PWPST_adv$AUC; objROC3_PWPST_adv$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.33%, achieved in 2264.65 secs (37.7 mins)


# -- Multi-threaded calculation of the # AUC from given start up to given prediction time 36 in following the CD-approach
# NOTE: Uses the superior Nearest Neighbour Estimator (NNE) method for S(t) with a 0/1-kernelNNE-kernel for S(t)
# NOTE2: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
ptm <- proc.time() #IGNORE: for computation time calculation;
predictTime <- 36
objROC4_PWPST_adv <- tROC.multi(datGiven=datCredit_valid, cox=cox_PWPST_adv, month_End=predictTime, sLambda=0.05, estMethod="NN-0/1", numDigits=4, 
                          fld_ID="PerfSpell_Key", fld_Event="PerfSpell_Event", eventVal=1, fld_StartTime="Start", fld_EndTime="TimeInPerfSpell",
                          graphName="DefaultSurvModel-Cox-PWPST-ROC_Depedendence_adv", genFigPath=paste0(genFigPath), 
                          caseStudyName=paste0("PWPST_", predictTime), numThreads=12, logPath=genPath)
objROC4_PWPST_adv$AUC; objROC4_PWPST_adv$ROC_graph
proc.time() - ptm
### RESULTS: AUC up to t: 96.40%, achieved in 3166.55 secs (52.8 mins)


# -- Store experimental objects | Memory optimisation
# PWPST-model: Basicc
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_03_bas"), objROC1_PWPST_bas);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_12_bas"), objROC2_PWPST_bas);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_24_bas"), objROC3_PWPST_bas);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_36_bas"), objROC4_PWPST_bas);

# PWPST-model: Advanced
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_03_adv"), objROC1_PWPST_adv);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_12_adv"), objROC2_PWPST_adv);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_24_adv"), objROC3_PWPST_adv);
pack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_36_adv"), objROC4_PWPST_adv);





# ----------------- 4. Create combined ROC-graph across multiple prediction times

# ------ Prentice-Williams-Peterson (PWP) Total-time definition | Basic-model

# - Ensure required objects exist in memory
if (!exists('objROC1_PWPST_bas')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_03_bas"), tempPath);gc()
if (!exists('objROC2_PWPST_bas')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_12_bas"), tempPath);gc()
if (!exists('objROC3_PWPST_bas')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_24_bas"), tempPath);gc()
if (!exists('objROC4_PWPST_bas')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_36_bas"), tempPath);gc()

# - Set ROC-parameters and initialize data structures
vecPercTimepoint <- c(3,12,24,36)
vecTROC <- list(objROC1_PWPST_bas, objROC2_PWPST_bas, objROC3_PWPST_bas, objROC4_PWPST_bas)
vLabels <- vector("list", length=length(vecPercTimepoint))

# -- Create a combined data object for plotting purposes
for (i in 1:length(vecPercTimepoint)) {
  # i <-1 # testing condition
  
  # datGraph <- data.frame(x = vFPR[-(nThresh+1)], y=vTPR[-1])
  
  # - Create a data object for the current prediction time
  if (i == 1) {
    datGraph <- data.table(PredictTime=paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                           x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR)
    
  } else {
    datGraph <- rbind(datGraph, 
                      data.table(PredictTime= paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                                 x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR))
  }
  vLabels[[i]] <- bquote("Prediction time "*italic(t)==.(vecPercTimepoint[i])*"; AUC: "*.(percent(vecTROC[[i]]$AUC, accuracy=0.01)))
}


# -- Graph a combined ROC-graph across prediction times t
# - Aesthetic parameters
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) model: Basic"]
vCol <- brewer.pal(8,"Set1")
vLabels_F <- setNames(vLabels, paste0(letters[1:length(vecPercTimepoint)],"_", vecPercTimepoint))
chosenFont <- "Cambria"

# - Create ROC-graph
(gg <- ggplot(datGraph, aes(x=x,y=y,group=PredictTime)) + theme_minimal() + 
    theme(text = element_text(family=chosenFont), legend.position="inside", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90),
          legend.position.inside = c(0.55,0.45),
          legend.background = element_rect(fill="snow2", color="black",
                                           linetype="solid", linewidth=0.1)) +
    labs(x = bquote("False Positive Rate "*italic(F^"+")), y = 
           bquote("True Positive Rate "*italic(T^"+"))) + 
    # Add 45-degree line
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "grey", linewidth=0.2) +
    # Main line graph
    geom_step(aes(x=x, y=y, linetype=PredictTime, colour=PredictTime), linewidth=0.05) + 
    geom_point(aes(x=x, y=y, shape=PredictTime, colour=PredictTime), size=0.25) +
    # Facets and scales
    facet_grid(FacetLabel ~ .) +  
    scale_color_manual(name=bquote("ROC"*(italic(t))), values=vCol, labels=vLabels) + 
    scale_linetype_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_shape_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_y_continuous(label=percent) + scale_x_continuous(label=percent))


# - Save graph
dpi <- 300
ggsave(gg, file=paste0(paste0(genFigPath, "DefaultSurvModel-Cox-PWPST-CombinedROC_Depedendence_bas.png")), 
       width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")




# ------ Prentice-Williams-Peterson (PWP) Total-time definition | Advanced-model

# - Ensure required objects exist in memory
if (!exists('objROC1_PWPST_adv')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_03_adv"), tempPath);gc()
if (!exists('objROC2_PWPST_adv')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_12_adv"), tempPath);gc()
if (!exists('objROC3_PWPST_adv')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_24_adv"), tempPath);gc()
if (!exists('objROC4_PWPST_adv')) unpack.ffdf(paste0(genPath,"DefaultSurvModel-Cox-PWPST-ROC_Depedendence_36_adv"), tempPath);gc()

# - Set ROC-parameters and initialize data structures
vecPercTimepoint <- c(3,12,24,36)
vecTROC <- list(objROC1_PWPST_adv, objROC2_PWPST_adv, objROC3_PWPST_adv, objROC4_PWPST_adv)
vLabels <- vector("list", length=length(vecPercTimepoint))

# -- Create a combined data object for plotting purposes
for (i in 1:length(vecPercTimepoint)) {
  # i <-1 # testing condition
  
  # datGraph <- data.frame(x = vFPR[-(nThresh+1)], y=vTPR[-1])
  
  # - Create a data object for the current prediction time
  if (i == 1) {
    datGraph <- data.table(PredictTime=paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                           x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR)
    
  } else {
    datGraph <- rbind(datGraph, 
                      data.table(PredictTime= paste0(letters[i], "_", vecPercTimepoint[i]), Threshold=vecTROC[[i]]$Thresholds, 
                                 x=vecTROC[[i]]$FPR, y=vecTROC[[i]]$TPR))
  }
  vLabels[[i]] <- bquote("Prediction time "*italic(t)==.(vecPercTimepoint[i])*"; AUC: "*.(percent(vecTROC[[i]]$AUC, accuracy=0.01)))
}


# -- Graph a combined ROC-graph across prediction times t
# - Aesthetic parameters
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) model: Advanced"]
vCol <- brewer.pal(8,"Set1")
vLabels_F <- setNames(vLabels, paste0(letters[1:length(vecPercTimepoint)],"_", vecPercTimepoint))
chosenFont <- "Cambria"

# - Create ROC-graph
(gg <- ggplot(datGraph, aes(x=x,y=y,group=PredictTime)) + theme_minimal() + 
    theme(text = element_text(family=chosenFont), legend.position="inside", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90),
          legend.position.inside = c(0.55,0.45),
          legend.background = element_rect(fill="snow2", color="black",
                                           linetype="solid", linewidth=0.1)) +
    labs(x = bquote("False Positive Rate "*italic(F^"+")), y = 
           bquote("True Positive Rate "*italic(T^"+"))) + 
    # Add 45-degree line
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "grey", linewidth=0.2) +
    # Main line graph
    geom_step(aes(x=x, y=y, linetype=PredictTime, colour=PredictTime), linewidth=0.05) + 
    geom_point(aes(x=x, y=y, shape=PredictTime, colour=PredictTime), size=0.25) +
    # Facets and scales
    facet_grid(FacetLabel ~ .) +  
    scale_color_manual(name=bquote("ROC"*(italic(t))), values=vCol, labels=vLabels) + 
    scale_linetype_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_shape_discrete(name=bquote("ROC"*(italic(t))), labels=vLabels) + 
    scale_y_continuous(label=percent) + scale_x_continuous(label=percent))


# - Save graph
dpi <- 300
ggsave(gg, file=paste0(paste0(genFigPath, "DefaultSurvModel-Cox-PWPST-CombinedROC_Depedendence_adv.png")), 
       width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


# - cleanup
suppressWarnings( rm(gg, vLabels, vLabels_F, vecTROC, datGraph, dat, 
                     objROC1_PWPST, objROC2_PWPST, objROC3_PWPST, objROC4_PWPST,
   cox_PWPST_adv,
   datCredit_train_PWPST, datCredit_valid_PWPST, datCredit_train, datCredit_valid
   ) )
