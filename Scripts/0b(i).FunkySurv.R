# ============================== SURVIVAL FUNCTIONS ==============================
# Defining bespoke customs relating to various generic aspects of survival analysis
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Marcel Muller, Dr Arno Botha, Bernard Scheepers
# VERSION: 2.0 (Apr-2025)
# DESCRIPTION: 
# This script defines various functions specific to survival modelling
# that are used elsewhere in this project or, indeed, used across other projects.
# Functions are grouped thematically.
# ================================================================================




# ----------------- 1. Functions related to the modelling process------------------



# --- Function to add and/or remove certain values from a vector (Improves readability of code)
# Input:  [mat]: (2 x n) Matrix on which changes will be made. The first column contains the variable names (vars) and the second their respective variable type (vartypes);
#         [Remove]: Vector containing the entries to be removed.
# Output:        [m], Updated matrix.
vecChange <- function(mat,Remove=FALSE,Add=FALSE){
  m <- mat
  if (all(Remove != FALSE)) {m <- mat[!(mat$vars %in% Remove)]} # Remove values in [Remove]
  if (all(Add != FALSE)) {m <- rbind(m,Add[!(Add$vars %in% mat$vars)])} # Add values in [Add]
  return(m)
}



# --- Function to detect significant correlations (abs(cor) > 0.6) between vectors.
# Input:  [data_train]: Training data; [variables]: List of variables used in correlation analysis;
#         [corrThresh]: Absolute correlation threshold above which correlations are deemed significant;
#         [method]: Correlation method.
# Output: <graph>  Upper half of correlation matrix.
#         <print> Text indicating the variable pairs with high correlation.
corrAnalysis <- function(data_train, variables, corrThresh = 0.6, method = 'spearman') {
  # Compute the correlation matrix
  corrMat <- as.data.table(data_train) %>% subset(select = variables) %>% cor(method = method)
  
  # Visualize the correlation matrix
  if(length(variables) <= 5){
    corrplot(corrMat, type = 'upper', addCoef.col = 'black', tl.col = 'black', diag=FALSE,
             tl.srt = 45, cl.pos="n")
  }
  
  # Find correlation coordinates exceeding the threshold
  corrCoordinates <- which(abs(corrMat) > corrThresh & abs(corrMat) < 1 & upper.tri(corrMat), arr.ind = TRUE)
  
  if(nrow(corrCoordinates) != 0){
    # Create a data table with correlation pairs
    corrProbs <- data.table(x = rownames(corrMat)[corrCoordinates[, 1]], y = colnames(corrMat)[corrCoordinates[, 2]])
    
    # Print the identified correlations
    for (i in 1:nrow(corrProbs)) {
      cat("Absolute correlations of ",percent(corrMat[corrProbs[i, x], corrProbs[i, y]]),
          " found for ", corrProbs[i, x], " and ", corrProbs[i, y],"\n")
    }
  }else{
    cat("No significant correlations were detected")
  }
}



# --- Function to return the appropriate formula object based on the time definition.
#         [TimeDef]: Time definition incorporated;
#         [variables]: List of variables used to build single-factor models;
TimeDef_Form <- function(TimeDef, variables, strataVar=""){
  # Create formula based on time definition of the dataset.
  if(TimeDef=="TFD"){# Formula for time to first default time definition (containing only the fist performance spell).
    formula <- as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ",
                                 paste(variables,collapse=" + ")))
    
  } else if(TimeDef=="AG"){# Formula for Andersen-Gill (AG) time definition
    formula <- as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ PerfSpell_Num + ",
                                 paste(variables,collapse=" + ")))
    
  } else if(TimeDef=="PWPST"){# Formula for Prentice-Williams-Peterson (PWP) Spell time definition
    formula <- as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ strata(", strataVar, ") + ",
                                 paste(variables,collapse=" + ")))
  } else {stop("Unkown time definition")}
  
  return(formula)
}



# --- Function to fit a given formula within a Cox regression model towards extracting Akaike Information Criterion (AIC) and related quantities
#         [formula]: Cox regression formula object; [data_train]: Training data;
#         [data_valid]: Validation data; [variables]: List of variables used to build single-factor models;
#         [it]: Number of variables being compared; [logPath], Optional path for log file for logging purposes;
#         [fldSpellID]: Field name of spell-level ID.
calc_AIC <- function(formula, data_train, variables="", it=NA, logPath="", fldSpellID="PerfSpell_Key") {
  # - Testing conditions
  # j <- 1; formula=TimeDef_Form(TimeDef,variables[j], strataVar=strataVar); 
  
  tryCatch({
    model <- coxph(formula,id=get(fldSpellID), data = data_train) # Fit Cox model
    
    if (!is.na(it)) {# Output the number of models built, where the log is stored in a text file afterwards.
      cat(paste0("\n\t ", it,") Single-factor survival model built. "),
          file=paste0(logPath,"AIC_log.txt"), append=T)
    }
    
    AIC <- AIC(model) # Calculate AIC of the model.
    
    # Return results as a data.table
    return(data.table(Variable = variables, AIC = AIC, pValue=summary(model)$coefficients[5]))
    
  }, error=function(e) {
    AIC <- Inf
    if (!is.na(it)) {
      cat(paste0("\n\t ", it,") Single-factor survival model failed. "),
          file=paste0(logPath,"AIC_log.txt"), append=T)
    }
    return(data.table(Variable = variables, AIC = AIC, pValue=NA)) 
  })
}



# --- Function to extract the Akaike Information Criterion (AIC) from single-factor models
# Input:  [data_train]: Training data; [data_valid]: [variables]: List of variables used to build single-factor models;
#         [fldSpellID]: Field name of spell-level ID; [TimeDef]: Time definition incorporated.
#         [numThreads]: Number of threads used; [genPath]: Optional path for log file. 
# Output: [matResults]: Result matrix.
aicTable <- function(data_train, variables, fldSpellID="PerfSpell_Key",
                      TimeDef, numThreads=6, genPath, strataVar="") {
  # - Testing conditions
   # data_train <- datCredit_train; TimeDef="PWPST"; numThreads=6
   # fldSpellID<-"PerfSpell_Key"; variables<-"g0_Delinq_SD_4"; strataVar="PerfSpell_Num_binned"
  
  # - Iterate across loan space using a multi-threaded setup
  ptm <- proc.time() #IGNORE: for computation time calculation
  cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
  cat("New Job: Estimating AIC for each variable as a single-factor survival model ..",
      file=paste0(genPath,"AIC_log.txt"), append=F)
  
  results <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                     .packages=c('data.table', 'survival'), .export=c('calc_AIC', 'TimeDef_Form')) %dopar%
    { # ----------------- Start of Inner Loop -----------------
      # - Testing conditions
      # j <- 1
      calc_AIC(formula=TimeDef_Form(TimeDef,variables[j], strataVar=strataVar), variables=variables[j],
                    data_train=data_train, it=j, logPath=genPath,  fldSpellID=fldSpellID)
      } # ----------------- End of Inner Loop -----------------
  stopCluster(cl.port); proc.time() - ptm  
  
  # Sort by concordance in ascending order.
  setorder(results, AIC)
  
  # Return resulting table.
  return(results)
}



# --- Function to fit a given formula within a Cox regression model towards extracting Harrell's C-statistic and related quantities
#         [formula]: Cox regression formula object; [data_train]: Training data;
#         [data_valid]: Validation data; [variables]: List of variables used to build single-factor models;
#         [it]: Number of variables being compared; [logPath], Optional path for log file for logging purposes;
#         [fldSpellID]: Field name of spell-level ID.
calc_HarrellC <- function(formula, data_train, data_valid, variables="", it=NA, logPath="", fldSpellID="PerfSpell_Key") {
  # formula <- TimeDef_Form(TimeDef,variables[j], strataVar=strataVar)
  
  tryCatch({
    model <- coxph(formula,id=get(fldSpellID), data = data_train) # Fit Cox model
    
    if (!is.na(it)) {# Output the number of models built, where the log is stored in a text file afterwards.
      cat(paste0("\n\t ", it,") Single-factor survival model built. "),
          file=paste0(logPath,"HarrelsC_log.txt"), append=T)
    }
    
    c <- concordance(model, newdata=data_valid) # Calculate concordance of the model based on the validation set.
    conc <- as.numeric(c[1])# Extract concordance
    sd <- sqrt(c$var)# Extract concordance variability as a standard deviation
    lr_stat <- round(2 * (model$loglik[2] - model$loglik[1]),0)# Extract LRT from the model's log-likelihood
    
    # Return results as a data.table
    return(data.table(Variable = variables, Concordance = conc, SD = sd, LR_Statistic = lr_stat))
  }, error=function(e) {
    conc <- 0
    sd <- NA
    lr_stat <- NA
    if (!is.na(it)) {
      cat(paste0("\n\t ", it,") Single-factor survival model failed. "),
          file=paste0(logPath,"HarrelsC_log.txt"), append=T)
    }
    return(data.table(Variable = variables, Concordance = conc, SD = sd, LR_Statistic = lr_stat)) 
  })
}



# --- Function to extract the concordances (Harrell's C) from single-factor models
# Input:  [data_train]: Training data; [data_valid]: Validation data;
#         [variables]: List of variables used to build single-factor models;
#         [fldSpellID]: Field name of spell-level ID; [TimeDef]: Time definition incorporated.
# Output: [matResults]: Result matrix.
concTable <- function(data_train, data_valid, variables, fldSpellID="PerfSpell_Key",
                      TimeDef, numThreads=6, genPath, strataVar="") {
  # - Testing conditions
  # data_valid <- datCredit_train_PWPST; TimeDef="PWPST"; numThreads=6
  # fldEventInd<-"Default_Ind"
  
  # - Iterate across loan space using a multi-threaded setup
  ptm <- proc.time() #IGNORE: for computation time calculation
  cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
  cat("New Job: Estimating B-statistic (1-KS) for each variable as a single-factor survival model ..",
      file=paste0(genPath,"HarrelsC_log.txt"), append=F)
  
  results <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                                 .packages=c('data.table', 'survival'), .export=c('calc_HarrellC', 'TimeDef_Form')) %dopar%
    { # ----------------- Start of Inner Loop -----------------
      # - Testing conditions
      # j <- 1
      calc_HarrellC(formula=TimeDef_Form(TimeDef,variables[j], strataVar=strataVar), variable=variables[j],
                    data_train=data_train, data_valid=data_valid, it=j, logPath=genPath,  fldSpellID=fldSpellID)
    } # ----------------- End of Inner Loop -----------------
  stopCluster(cl.port); proc.time() - ptm  
  
  # Sort by concordance in descending order.
  setorder(results, -Concordance)
  
  # Return resulting table.
  return(results)
}



# --- Function to calcualte the complement of the KS test statistic "B-statistic"
# Inputs: [formula]: Cox regression formula object; [data_train]: Training data;
#         [fldSpellID]: Field name of spell-level ID; [vEvents]: spell-level vector of event indicators
#         [seedVal]: Seed value for random number generation; [it]: optional iteration parameter for logging purposes;
#         [logPath]: Optional path for log file for logging purposes.
# Outputs: b-statistic (single value)
calcBStat <- function(formula, data_train, fldSpellID="PerfSpell_Key", vEvents, seedVal, it=NA, logPath=NA) {
  # - Testing conditions
  # formula <- TimeDef_Form(TimeDef,variables[j], strataVar=strataVar)
  
  # Fit Model
  tryCatch({
    model <- coxph(formula, data = data_train, id=get(fldSpellID)) 
    
    if (!is.na(it)) {
      cat(paste0("\n\t ", it,") Single-factor survival model built. "),
          file=paste0(logPath,"BStat_log.txt"), append=T)
    }
    
    # Calculate Cox-Snell (adjusted) residuals
    vCS <- calc_CoxSnell_Adj(model, vIDs=data_train[[fldSpellID]], vEvents=vEvents)
    # Initialize a unit exponential distribution
    set.seed(seedVal, kind = "Mersenne-Twister")
    vExp <- rexp(length(vCS),1)
    # Perform the two-sample Kolmogorov-Smirnov test of distribution equality
    #   H_0: vCS and vExp originates from the same distribution
    #   NOTE: We only desire the KS test statistic in measuring distributional dissimilarity
    #   Then, we subtract this from 1 in creating a coherent statistic; greater is better
    bStat <- 1 - round(suppressWarnings(ks.test(vCS,vExp))$statistic,4)
    return(bStat)
    
  }, error=function(e) {
    bStat <- 0
    if (!is.na(it)) {
      cat(paste0("\n\t ", it,") Single-factor survival model failed. "),
          file=paste0(logPath,"BStat_log.txt"), append=T)
    }
    return(bStat) 
  })
}



# --- Function to extract the B-statistic from a range of models built on a list of variables based on a time definition.
# Input:  [data_train]: Training data; [seedVal]: Seed value for random number generation;
#         [numIt]: Number of iterations; [TimeDef], Time definition incorporated.
#         [fldSpellID]: Field name of spell-level ID; [fldLstRowInd]: Indicates the end of a performance spell;
#         [fldEventInd]: Indicates whether the target event occured; [numThreads]: Number of threads;
#         [genPath]: Optional path for log file.
# Output: [Results]:  Results table.
csTable <- function(data_train, variables, TimeDef, seedVal=1, numIt=5, 
                    fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="DefaultStatus1",
                    numThreads=6, genPath=NA, strataVar=""){
  
  # - Testing conditions
  # data_train <- datCredit_train_PWPST; variables<-vars2; TimeDef<-"PWPST"; seedVal<-1; numIt<-5; 
  # fldLstRowInd="PerfSpell_Exit_Ind";  fldSpellID="PerfSpell_Key"; fldEventInd="Default_Ind"; numThreads=6
  
  # - Initialize results
  results <- data.frame(Variable = variables, B_Statistic = NA_real_)
  
  # - Data preparation
  # Subset last row per performing spell for Goodness-of-Fit (GoF) purposes
  datLstRow <- copy(data_train[get(fldLstRowInd)==1,])
  vLstRow_Events <- datLstRow[, get(fldEventInd)]
  
  # - Simulate null distribution if seedVal is not NA
  if (!is.na(seedVal)) {
    
    
    # - Iterate across loan space using a multi-threaded setup
    ptm <- proc.time() #IGNORE: for computation time calculation
    cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
    cat("New Job: Estimating B-statistic (1-KS) for each variable as a single-factor survival model ..",
        file=paste0(genPath,"BStat_log.txt"), append=F)
    
    results$B_Statistic <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                      .packages=c('data.table', 'survival'), .export=c('calc_CoxSnell_Adj', 'calcBStat', 'TimeDef_Form')) %dopar%
      
      { # ----------------- Start of Inner Loop -----------------
        # - Testing conditions
        # var <- variables[1]; j<-4
        calcBStat(formula=TimeDef_Form(TimeDef,variables[j], strataVar=strataVar), data_train=data_train, fldSpellID=fldSpellID, vEvents=vLstRow_Events,
                  seedVal=seedVal, it=j, logPath=genPath)
        
      } # ----------------- End of Inner Loop -----------------
    stopCluster(cl.port); proc.time() - ptm
    
    # Sort results by B statistic in descending order
    results <- results[order(-results$B_Statistic, na.last = TRUE), ]
    
    # Return results and range of B statistics
    return(list(Results = results, Range = diff(range(results$B_Statistic, na.rm = TRUE))))
    
  } else {
    # Perform iterative B calculation when seedVal is NA
    # Initialize Results matrix to contain the number of interations
    matResults <- matrix(NA, nrow = length(variables), ncol = numIt,
                         dimnames = list(variables,
                                         paste0("Iteration_", 1:numIt))) %>%
                          as.data.table()
    
    # - Iterate across loan space using a multi-threaded setup
    ptm <- proc.time() #IGNORE: for computation time calculation
    cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
    cat("New Job: Estimating B-statistics (1-KS) ..",
        file=paste0(genPath,"BStat_log.txt"), append=F)
    
    for (it in seq_len(numIt)) {
      
      cat(paste0("\n Estimating B-statistic (1-KS) for each variable as a single-factor survival model for iteration ", it, " .."),
          file=paste0(genPath,"BStat_log.txt"), append=T)
      
      matResults[, it] <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                                     .packages=c('data.table', 'survival'), .export=c('calc_CoxSnell_Adj', 'calcBStat', 'TimeDef_Form')) %dopar%
        
        { # ----------------- Start of Inner Loop -----------------
          # - Testing conditions
          # var <- variables[1]
          calcBStat(formula=TimeDef_Form(TimeDef,variables[j], strataVar=strataVar), data_train=data_train, fldSpellID=fldSpellID, vEvents=vLstRow_Events,
                    seedVal=seedVal*it, it=j, logPath=genPath)
          
        } # ----------------- End of Inner Loop -----------------
    }
    stopCluster(cl.port); proc.time() - ptm
    
    # Compute additional statistics for the results matrix
    colRanges <- matResults[, lapply(.SD, function(x) diff(range(x, na.rm = TRUE)))] # Calculate the Range of B-statistic value for each iteration
    matResults[, Average := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^Iteration_")]# Calculate the average B-statistic for each variable
    matResults[, Variable := variables] 
    matResults <- matResults %>% relocate(Variable, .before=Iteration_1)
    
    #matResults <- cbind(Variables = c(variables,"Range"),matResults)# Add a column to the Results matrix to cross reference the variables with their respective B-statistics
    setorder(matResults,-Average)# Arrange matrix according to average
    
    # Return matrix of B statistics
    return(list(Results=matResults, IterationRanges=colRanges))
  }
}



# --- Function to calculate various survival-related quantities for a given loan history
# Input:    [datGiven]: given loan history; [coxGiven]: fitted cox PH model; [it]: current iteration index; 
#           [numKeys]: total keys; [genPath]: Optional path for log file.
# Output:   Survival probability, cumulative hazard, hazard, and event probability
survQuants <- function(datGiven, coxGiven, it=1, numKeys, genPath="", timeVar="End", startVar="Start") {
  # datGiven <- subset(datCredit_valid_TFD,PerfSpell_Key == vSpellKeys[j]); coxGiven <- cox_TFD
  # it=1; numKeys <- numSpellKeys; timeVar="End"; startVar="Start"
  
  # - Add a row when scoring S(t), merely to facilitate estimation of f(t)
  if (datGiven[,get(timeVar)][1] > 1) {
    datAdd <- datGiven[1, ]
    datAdd[, (startVar) := get(startVar) - 1]
    datAdd[, (timeVar) := get(timeVar) - 1]
    datGiven <- rbind(datAdd, datGiven)
    addedRow <- T
  } else addedRow <- F
  
  # - Compute individual survival curve from fitted Cox model
  survFit_pred <- survfit(coxGiven, centered=F, newdata=datGiven, id=PerfSpell_Key)
  
  cat("\n\t", it, "of", numKeys, "| Estimation completed for spell key:", unique(datGiven$PerfSpell_Key),
      file=paste0(genPath,"survQuants_log.txt"), append=T)
  
  # - Compile survival table
  datSurv <- data.table(PerfSpell_Key = unique(datGiven$PerfSpell_Key), End=datGiven$End, # composite key
                        CHaz=survFit_pred$cumhaz, #RiskSetSize=survFit_pred$n.risk,
                        #NumEvents=survFit_pred$n.event, NumCensored=survFit_pred$n.censor,
                        Survival=round(survFit_pred$surv,digits=15))
  # plot(survFit_pred)
  # plot(datSurv$Survival, type="b")
  
  # - Due to data irregularities, the survival probability can increase again over certain t,
  # which breaks the axioms of probability since S(t) = 1 - F(t) and F(t) is monotonically increasing by definition
  # Luckily, these cases seem isolated to the later parts of the spell life
  # In these cases, force the S(t)-estimate to equal the previous
  vecErrors <- which(diff(datSurv$Survival) > 0)
  if (length(vecErrors) > 1) {
    
    while (looping==T) {
      datSurv$Survival[vecErrors[1]+1] <- datSurv$Survival[vecErrors[1]]
      # Testing end condition
      vecErrors <- which(diff(datSurv$Survival) > 0)
      if (length(vecErrors) > 1) {looping=T} else {looping=F}
    }
  }
  # plot(datSurv$Survival, type="b")
  
  # - Approximate baseline hazard h_0(t) from cumulative baseline hazard H_0(t)
  #datSurv[, Hazard := c(datSurv$CHaz[1], diff(datSurv$CHaz))]
  datSurv[, Hazard := (shift(Survival,n=1,type="lag",fill=1) - Survival)/shift(Survival,n=1,type="lag",fill=1)]
  datSurv[, EventProb := shift(Survival,n=1,type="lag",fill=1) * Hazard] # f(t|X) = S(t-1|X) . h(t|X)
  
  # - Remove added row (if added)
  if (addedRow) {
    datSurv <- datSurv[2:NROW(datSurv),]
  }
  
  return(datSurv)
}



### AB: Given the work of Bernard, I'm no longer sure of the utility of the below.

# --- function to compute the (unscaled) Schoenfeld residuals for a Cox PH model
#   1) Schoenfeld residuals are computed for each specified variable in the training dataset
#   2) Tests are conducted for significance of the correlation of the Schoenfeld residuals against time
#   3) Graphs are created of the Schoenfeld residuals against time
# Input:  cph - Cox PH model to be assessed
#         dat_train - Dataset used to train the Cox PH model
#         var - Name of the variables for which the Schoenfeld residuals are to be computed
#         id - Name of the column uniquely identifying each row of the dat_train
#         time - Name of the column identifying the associated time in dat_train
#         status - Name of the column identifying the target variable in dat_train
#         verbose - Indicator variable for supressing graphs created by the function
#         max_time - The maximum time for which the graph should be plotted
# Output: A list containing the following:
#         data - A dataset containing the id, time, raw variable's value, and the associated Schoenfeld residual
#         CorTest - A list of the results from test(s) for significance of the correlation of the Schoenfeld residuals against time
#         plots - A list of graphs of the Schoenfeld residuals against time
cph_schoen <- function(cph, var=NULL, dat_train, id, time, status, verbose=T, max_time=NULL){
  # UNIT TEST (VARIABLE INITIALISATION)
  # cph <- cph_Default_PH_test; var <- c("Principal_wins"); id <- "PerfSpell_Key"; time <- "TimeInPerfSpell"; status <- "DefaultStatus1"; max_time <- 240; verbose <- F
  
  # Copying the training dataset to ensure no contamination and ensuring that the dataset is of the correct class (this step usually takes a considerable amount of time for large datasets)
  if (any(class(dat_train) %in% "tbl_df")){
    dat_train2 <- as.data.table(dat_train)
  } else {
    dat_train2 <- copy(dat_train)
  }
  if (any(class(dat_train2) %in% "grouped_df")){
    dat_train2 <- ungroup(dat_train2)
  }
  
  cph_sum <- summary(cph) # Getting a summary of the model
  row_names <- rownames(cph_sum$coefficients) # Getting the names of all the variables in the Cox model (includes levels of categorical variables)
  
  # Object for graphs
  gplots <- list()
  
  # Object for correlation test
  cor_test <- list()
  
  # Getting all the variable names in the model if no names are specified
  if (is.null(var)){
    var <- unlist(strsplit(toString(cph_sum$call$formula[[3]]), '[,+ ]+'))
    var <- var[var!=""]
  }
  
  # Initialising the dataset to be returned
  dat_return <- data.table(ID = numeric(),
                           Time = numeric(),
                           Var_Val = numeric(),
                           Sch_Res = numeric(),
                           Var_Name = as.character(),
                           Var_Name_Base = as.character())
  
  # Initialising the temporary dataset used in the loop for all variables
  col_names <- c(id, time, status, var) # Getting the names of the columns from the training dataset to be subsetted
  dat_temp <- dat_train2[, ..col_names] # Subsetting from the main dataset
  colnames(dat_temp)[1:3] <- c("ID", "Time", "Status") # Renaming the columns for conveinience
  dat_temp <- cbind(dat_temp, predict(cph, dat_train2, type="risk")); colnames(dat_temp)[length(colnames(dat_temp))] <- "RiskScore" # Getting the risk score of each observation
  dat_temp <- merge(dat_temp, dat_temp[, list(SumScore=sum(RiskScore)), by=list(Time)], by="Time", all.x=T) # Getting the total risk score at each time point and merging it back into the dataset
  
  # Computing the Schoenfeld residuals for each selected variable and appending the dataset which to return
  k <- 0 # Counting variable for lists in for loop below
  for (i in 1:length(var)){
    # i <- 1
    var_type <- class(dat_temp[[var[i]]]) # The class of the selected variable
    var_name <- var[i]
    
    col_names2 <- c("ID", "Time", "Status", "RiskScore", "SumScore", var_name)
    dat_temp2 <- dat_temp[, ..col_names2] # Subsetting from the temporary dataset to only include the i'th variable's values
    colnames(dat_temp2)[6] <- "Var_Val"
    dat_temp2[, Var_Name := var_name]; dat_temp2[, Var_Name_Base := var_name] # [Var_Name] is the name of the variable as in the training dataset; [Var_Name_Base] is the level of the categorical variable (equal to Var_Name for numeric variables)
    
    # Computing the Schoenfeld residuals for a numeric variable
    if (var_type == "numeric"){
      # Computing the risk-weighted value of the chosen variable
      dat_temp2[, RW_Val := Var_Val*RiskScore/SumScore]
      # Computing the expected value of the numerical variable at each time point and merging it back into the dataset
      dat_temp2 <- merge(dat_temp2, dat_temp2[, list(Exp_Val=sum(RW_Val)), by=list(Time)], by="Time", all.x = TRUE)
      
      # Subsetting to only include defaulted accounts || Change naming accordingly as structure of data changes
      dat_temp3 <- dat_temp2[Status==1,]; rm(dat_temp2)
      # Computing the Schoenfeld residuals
      dat_temp3[, Sch_Res := Var_Val - Exp_Val] # -447886.85507
      
      # Appending the temporary dataset to the main dataset that is to be returned
      dat_return <- rbind(dat_return, dat_temp3[, list(ID, Time, Var_Val, Sch_Res, Var_Name, Var_Name_Base)])
      
      # Updating the counter variable
      k <- k + 1
      
      # Plotting if verbose = True
      if (verbose==F){
        if (is.null(max_time)){
          max_time <- max(dat_return$Time)
        } # if

        g <- ggplot(dat_temp3, aes(x=Time, y=Sch_Res)) +
             theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
             theme(text=element_text(family=chosenFont),legend.position = "bottom",
                   strip.background=element_rect(colour="grey", fill="#D3D3D3")) +
             facet_wrap(~ Var_Name, strip.position = "right") +
             geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
             scale_x_continuous(limits = c(NA, max_time)) +
             scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
        
        g_name <- paste0("SchRes_",var[i])
        
        gplots[[g_name]] <- g; names(gplots[k])  
      } # if

      # Correlation Test (Formal)
      cor_test[[k]] <- cor.test(dat_temp3$Time, dat_temp3$Sch_Res, method = "spearman"); names(cor_test)[k] <- var[i] # Spearman rank correlation is used for robustness
      
      # Clean up
      rm(dat_temp3)
      
    } else if (var_type %in% c("character", "factor")){
      levels_n <- length(grep(var[i], row_names)) # Computing the number of levels (-1) of the categorical variable | k-1 levels of the categorical variable
      levels <-  substring(row_names[grep(var[i], row_names)], nchar(var[i])+1) # Getting the levels of the categorical variable (-1)
      
      # Computing the Schoenfeld residuals for each level of the categorical variable
      for (j in 1:levels_n){
        # j <- 1
        dat_temp3 <- copy(dat_temp2) # Copying the temporary dataset to increase efficiency with multiple levels of a categorical variable
        dat_temp3[, Var_Name := levels[j]] # Setting [Var_Name] to level j of the categorical variable
        
        # Computing the risk-weighted value of the chosen variable
        dat_temp3[, RW_Val := as.numeric(Var_Val==levels[j])*RiskScore/SumScore]
        # Computing the expected value of the chosen variable at each time point and merging it back into the dataset
        dat_temp3 <- merge(dat_temp3, dat_temp3[, list(Exp_Val=sum(RW_Val)), by=list(Time)], by="Time", all.x = TRUE)
        # Subsetting to only include defaulted accounts
        dat_temp4 <- dat_temp3[Status==1,]; rm(dat_temp3)
        # Computing the Schoenfeld residuals
        dat_temp4[, Sch_Res := as.numeric(Var_Val==levels[j]) - Exp_Val]
        
        # Appending the temporary dataset to the main dataset
        dat_return <- rbind(dat_return, dat_temp4[, list(ID, Time, Var_Val, Sch_Res, Var_Name, Var_Name_Base)])
        
        # Updating the counter variable
        k <- k+1
        
        # Plotting if verbose = True
        if (verbose==F){
          if (is.null(max_time)){
            max_time <- max(dat_return$Time)
          } # if
          
          # Modification for faceting to plotting dataset
          dat_temp4[, fac := paste0(Var_Name_Base, " ~ ", Var_Name)]
          
          g <- ggplot(dat_temp4, aes(x=Time, y=Sch_Res)) +
               theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
               theme(text=element_text(family=chosenFont),legend.position = "bottom",
                     strip.background=element_rect(colour="grey", fill="#D3D3D3")) +
               facet_wrap(~ fac, strip.position = "right") +
               geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
               scale_x_continuous(limits = c(NA, max_time)) +
               scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
          g_name <- paste0("SchRes_",var[i], levels[j])  
          gplots[[g_name]] <- g; names(gplots[k])    
        } # if
        
        # Correlation Test
        cor_test[[k]] <- cor.test(dat_temp4$Time, dat_temp4$Sch_Res, method = "spearman"); names(cor_test)[k] <- paste(var[i], levels[j]) # Spearman rank correlation is used for robustness
      }
    } # else if
    
  } # for
  
  dat_return_wider <- dat_return %>% pivot_wider(names_from = c(Var_Name_Base, Var_Name), values_from = c(Var_Val, Sch_Res)) %>% setDT()
  
  # Small correction in naming
  names(dat_return_wider)[grep("Var_Val*", names(dat_return_wider))[names(dat_return_wider)[grep("Var_Val*", names(dat_return_wider))] == paste0("Var_Val", "_", var, "_", var)]] <- paste0("Var_Val_", var[sapply(dat_train2[,..var], is.numeric)])
  names(dat_return_wider)[grep("Sch_Res*", names(dat_return_wider))[names(dat_return_wider)[grep("Sch_Res*", names(dat_return_wider))] == paste0("Sch_Res", "_", var, "_", var)]] <- paste0("Sch_Res_", var[sapply(dat_train2[,..var], is.numeric)])
  
  # Clean up
  suppressWarnings(rm(dat_train2, dat_temp, dat_temp2, dat_temp3, dat_temp4, var, var_type, var_name, col_names, g, g_name, k, levels, levels_n))

  return(list(data = dat_return_wider, CorTest = cor_test, plots = gplots))
  
  # rm(dat_temp, dat_temp2, dat_return, id, var, time, status, max_time, verbose)
} # function

# --- UNIT TEST (Breslow Approximation)
# - Setup
# method <- "breslow"
# cph <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1) ~
#              Principal + Instalment + LN_TPE + slc_pmnt_method,
#              id= PerfSpell_Key, data = dat_train %>% group_by(PerfSpell_Key, TimeInPerfSpell),
#              ties = method)
# cph_bres_scaled_schoenfeld <- residuals(cph, type="schoenfeld") # | Use residuals() for classical Schoenfeld residuals and not Scaled; Schoefeld residuals computed for each level

# - Function Execution
# return <- cph_schoen(cph=cph, id = "PerfSpell_Key", dat_train = dat_train, time = "TimeInPerfSpell", 
#             status = "DefaultStatus1", verbose = F, max_time = 240)

# - Comparison of Schoenfeld residuals for a numerical variable
# a <- cph_bres_scaled_schoenfeld[,1]
# b <- return$data[, Sch_Res_Principal]
# ab <- data.table(Sch_Res_Principal_Residuals = a,
#                  Sch_Res_Principal_CustFunc = b)
# all.equal(ab$Sch_Res_Principal_Residuals, ab$Sch_Res_Principal_CustFunc)
### RESULTS:~ TRUE
# - Plot comparison for numeric variable
# par(mfcol=c(1,2)); plot(x=rownames(cph_bres_scaled_schoenfeld), y=cph_bres_scaled_schoenfeld[,1], xlim=c(0,240)); plot(x=return$data$Time, y=return$data$Sch_Res_Principal, xlim=c(0,240))

# - Comparison of Schoenfeld residuals for a categorical variable
# c <- cph_bres_scaled_schoenfeld[,4]
# d <- as.numeric(return$data$`Sch_Res_slc_pmnt_method_Debit Order other bank`)
# cd <- data.table(Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_Residuals = c,
#                  Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_CustFunc = d)
# all.equal(cd$Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_Residuals,cd$Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_CustFunc)
### RESULTS:~ TRUE
# - Plot comparison for categorical variable
# par(mfcol=c(1,2)); plot(x=rownames(cph_bres_scaled_schoenfeld), y=cph_bres_scaled_schoenfeld[,8], xlim=c(0,240)); plot(x=return$data$Time, y=return$data$Sch_Res_slc_pmnt_method_Suspense, xlim=c(0,240))
