# ============================== SURVIVAL FUNCTIONS ==============================
# Defining custom functions used across various projects
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Marcel Muller

# VERSION: 1.0 (July-2023)
# DESCRIPTION: 
# This script defines various functions specific to survival modelling
# that are used elsewhere in this project or, indeed, used across other projects.
# Functions are grouped thematically.
# ================================================================================

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
# return <- cph_schoen(cph=cph, id = "PerfSpell_Key", dat_train = dat_train, time = "TimeInPerfSpell", status = "DefaultStatus1", verbose = F, max_time = 240)

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


# --- function to return the time-dependent ROC curve and time-dependent AUC value (Incidence/ Dynamic) of a Cox survival model
# Input:  dat_explain - The dataset containing the required variables for time-dependent ROC and AUC computation
#         predict.time - Time point for valuation
#         Input_Names - A vector containing the names in dat_explain
#             First element is the name of the column containing the event time
#             Second element is the name of the column containing the event indicator
#          marker - A vector containing the predictions of the Cox Model on the records in dat_explain
#             Third element is the name of the column containing the marker values (predictions); we work specifically with the linear prediction
# Output: A list containing time-dependent marker values, TPs, FPs, and AUCs.
risksetROC_helper <- function(dat_explain, Input_Names = c("TimeInPerfSpell", "DefaultStatus1"), marker = cph_lp, predict.time=12) {
  risksetROC::risksetROC(Stime        = dat_explain[[Input_Names[1]]],  
                         status       = dat_explain[[Input_Names[2]]],                 
                         marker       = marker,                             
                         entry        = dat_explain[[Input_Names[1]]]-1,                              
                         predict.time = predict.time,
                         method       = "Cox",
                         plot = FALSE)
}


# --- function for computing various assessment measures for a given Cox Proportional Hazards model.
#   1) A test can be conducted for proportional hazards and the Schoenfeld residuals for each covariate is plotted.
#   2) The VIF can be computed for each covariate.
#   3) An analysis of the reduction of deviance (ANOVA) by sequnetial addition of covariates to the model can be conducted and a visual aid provided.
#   4) The time-depdendent ROC-curves can be plotted and the time-dependent AUC-values computed for various time-points.
# Input: dat_train - The dataset used to train the Cox PH model
#        dat_explain - The dataset with which the Cox PH model needs to be assessed.
#        model - The Cox PH model to be assessed.
#        input.l - list containing the relevant variable names in dat_train (& dat_explain) to enable the various assessments
#             input.l <- list(Time of Observation Variable Name, Event Indicator Variable Name, c(Variable Name 1,..., Variable Name n))
#        PH_assumption - Indicator for conducting the PH assumption test and displaying the Schoenfeld residuals (only if verbose  = FALSE).
#        VIF - Indicator for computing the VIF of the covariates.
#        anova_explain - Indicator for conducting an analysis of the reduction of the deviance of the model by sequential addition of variables and displaying the results (if verbos =FALSE)
#        predict.time - Vector containg the desired times to compute the AUC values and plot the corresponding ROC curves (if verbose = FALSE)
#        verbose - Indicator variable used to supress visual graphs created by some of the performance measures and analysis.
# Output: List of the various test and analysis resutls, the elements of this list are as follows:
#         Model_VIF - The VIF of each covariate
#         Model_ANOVA - The results of a reducion in deviance analysis by the sequnetial addition of variables
#         Model_AIC - The AIC of the fitted Cox PH model
#         AUC - The corresponding time-dependent AUC values of the given prediction times
#         PH_Test - Graphical- and statistical tests for proportional hazards of each of the given covariates (defaults to all covariates in the model). Also included is a dataset of the original values and the residuals.
swiss_model <- function(dat_train, dat_explain, model, input.l = NULL, PH_assumption = TRUE, VIF = TRUE, anova_explain = FALSE, AUC_explain = FALSE, predict.time = 12, verbose = FALSE, max_time = NULL){
  # dat_explain<-dat_valid; input.l<-list("PerfSpell_Key", "TimeInPerfSpell", "DefaultStatus1"); model<-cph; PH_assumption<-TRUE; VIF<-TRUE; anova_explain<-FALSE; AUC_explain<-TRUE; predict.time<- c(12); verbose<-FALSE; max_time<-240
  # - Initialise the output vector
  output <- NULL
  
  if (length(input.l) < 4){
    input.l[[4]] <- as.vector(unlist(strsplit(toString(summary(model)$call$formula[[3]]), '[,+ ]+')))
    input.l[[4]] <- input.l[[4]][input.l[[4]]!=""]
  }
  
  # --- Model Assumptions
  # - Multicollinearity
  if (VIF){
    Model_VIF <- vif(model)
    output[["VIF"]] = Model_VIF
  }
  
  # - Proportional Hazards Assumption
  # Proportional Hazards Test
  if (PH_assumption){
    PH_Test <- cph_schoen(cph=model, dat_train = dat_train, var = input.l[[4]], id = input.l[[1]], time = input.l[[2]], status = input.l[[3]], verbose = F, max_time = max_time)
    output[["PH_Test"]] <- PH_Test
  }
  
  # --- Variable Importance Measures
  # - ANOVA Analysis - Reduction in Deviance
  # Conduct the ANOVA analysis, conditional on if it is required
  if (anova_explain){
    Model_ANOVA <- anova(model, test = 'chisq')
    ANOVA <- list(anova = Model_ANOVA)
    
    # Create the ANOVA plot
    if (!verbose){
      # Getting the names of the relevant variables in the ANOVA object.
      col_names <- input.l[[4]]
      
      # Creating a dataset from the ANOVA object for graphing.
      datPlot <- data.frame(Variable_Name = col_names, ChiSq_Value = round(Model_ANOVA$Chisq[2:(length(col_names)+1)]), Var_Order = 1:length(col_names))
      
      # Creating a graph to visually present the contribution of each variable in the Cox model.
      g_anova <- ggplot(datPlot, aes(x=reorder(Variable_Name, as.numeric(ChiSq_Value)), y=as.numeric(ChiSq_Value))) +
        geom_col(col='blue', fill='blue') + 
        theme_minimal() + theme(plot.title = element_text(hjust=0.5)) + coord_flip() +
        labs(x="Variable", y="Reduction in Model's Deviance") +
        geom_label(aes(label=Var_Order), position = position_stack(vjust=0.5))
      
      ANOVA[["plots"]] <- g_anova
    }
    
    output[["anova"]] = ANOVA
  }
  
  # --- Model Assessment
  # - AIC Value
  output[["AIC"]] <- extractAIC(model)[2]
  
  # - Time-dependent AUC Values and ROC Curves
  if (AUC_explain){
    # Computing the linear predictions for each recored in the given dataset and amending the dataset to include those predictions
    cph_lp <- predict(model, dat_explain, type="lp") # linear predictions are used for assessment (X*Beta)
    
    # Creating a dataset containing the FPs and TPs and AUC values for each time in the vector
    risksetROC_data <- data.frame(Predict_Time = NULL, AUC = NULL, FP = NULL, TP = NULL, Marker = NULL)
    for (i in 1:length(predict.time)){
      temp <- risksetROC_helper(dat_explain = dat_explain, Input_Names = (c(input.l[[2]], input.l[[3]])), marker = cph_lp, predict.time = predict.time[i])
      risksetROC_data <- rbind(risksetROC_data, data.frame(Predict_Time = predict.time[i], AUC = temp$AUC, FP = temp$FP, TP = temp$TP))
    }
    rm(temp)
    risksetROC_data <- risksetROC_data %>% arrange(Predict_Time, FP, TP) # Arranging the dataset according to the given vector of prediction times (and then by FPs and TPs)
    
    AUC <- list(data = data.table(predict_times = predict.time,
                                  auc = unique(risksetROC_data$AUC)))
    
    # Plotting the ROC curves for the desired prediction times (if verbose = FALSE)
    if (!verbose){
      g_auc <- ggplot(data = risksetROC_data, mapping = aes(x = FP, y = TP)) +
        geom_line(col = 'blue') +
        geom_abline(col='red') +
        geom_label(data = risksetROC_data %>% dplyr::select(Predict_Time,AUC) %>% unique,
                   mapping = aes(label = sprintf("%.3f", AUC)), x = 0.5, y = 0.5) +
        facet_wrap(vars(Predict_Time)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5),
              strip.background = element_blank()) +
        xlab("FP") + ggtitle('ROC Curve(s) for the given CPH Model')
      
      AUC[["plots"]] <- g_auc
      output[["AUC"]] <- AUC
    }
  }
  
  return(output)
  # rm(dat_explain, model, PH_assumption, VIF, anova_explain, AUC_explain, input.l, predict.time, verbose, risksetROC_data, datPlot, Model_VIF, PH_Test, output, temp)
}

# --- Unit Test
# swiss_cph <- swiss_model(dat_train = dat_train, dat_explain = dat_valid, model = cph, input.l = list("PerfSpell_Key", "TimeInPerfSpell", "DefaultStatus1"), PH_assumption = T, VIF = T, anova_explain = F, AUC_explain = T, predict.time = 12, verbose = F, max_time = 240)  




# --- function for identifying FALSE default spells in a given dataset
# Input: dat_given - The dataset in which false default spells should be identified
#        LoanID - Name of the account ID
#        PerfSpellID - Name of the performance spell ID of an account
#        DefSpellID - Name of the default spell ID of an account
#        Counter - Name of the counter variable (counts the observation of a account)
#        PerfSpell_Counter - Name of the counter variable for an assocaited performance spell
#        DefSpell_Counter - Name of the counter variable for an assocaited default spell
#        DefSpellResol_Type_Hist - Name of variable indicating how the default spell was resolved
#        PerfSpell_Max_Date - Name of variable indicating the last observed date of the assocaited performance spell (optional)
#        DefSpell_Max_Date - Name of variable indicating the last observed date of the assocaited performance spell (optional)
# Output: dat_given - The dataset conatining the variable identifying FALSE default spells
False_Perf_Def <- function(dat_given, LoanID=NA, Date=NA, PerfSpellID=NA, DefSpellID=NA, Counter=NA, PerfSpell_Counter=NA, DefSpell_Counter=NA,
                           PerfSpell_Max_Date=NA, DefSpell_Max_Date=NA){
  # --- Unit test parameters:
  # dat_given <- copy(datCredit_smp); LoanID <- "LoanID"; Date <- "Date"; PerfSpellID <- "PerfSpell_Key"; DefSpellID <- "DefSpell_Key"
  # Counter <- "Counter"; PerfSpell_Counter <- "PerfSpell_Counter"; DefSpell_Counter <- "DefSpell_Counter"; DefSpell_Max_Date <- "DefSpell_Max_Date"; PerfSpell_Max_Date <- "PerfSpell_Max_Date"
  
  # --- Arranging the dataset according to the LoanID and Date
  dat_given <- arrange(dat_given, get(deparse(substitute(LoanID))), get(deparse(substitute(Date)))) %>% setDT(key=c(get(deparse(substitute(LoanID))), get(deparse(substitute(Date)))))
  
  # --- Creating a subset containing only the required column names
  colnames <- c(LoanID, Date, PerfSpellID, DefSpellID, Counter, PerfSpell_Counter, DefSpell_Counter,
                ifelse(!is.na(DefSpell_Max_Date), DefSpell_Max_Date, NA), ifelse(!is.na(PerfSpell_Max_Date), PerfSpell_Max_Date, NA))
  dat_sub <- subset(dat_given, select=colnames[!is.na(colnames)])
  
  # --- Renaming the columns to enable easier coding
  colnames_new <- c("LoanID", "Date", "PerfSpellID", "DefSpellID", "Counter", "PerfSpell_Counter", "DefSpell_Counter",
                    ifelse(!is.na(DefSpell_Max_Date), "DefSpell_Max_Date", NA), ifelse(!is.na(PerfSpell_Max_Date), "PerfSpell_Max_Date", NA))
  colnames(dat_sub) <- colnames_new[!is.na(colnames_new)]
  
  # --- Creating a variable showing the previous/ next counter value
  dat_sub[, Counter_Prev := as.numeric(shift(x=Counter, n=1, type="lag")), by=LoanID] # Creating a temporary variable for checking if the next observation of the account is in the dataset (exists)
  dat_sub[, Counter_Next := as.numeric(shift(x=Counter, n=1, type="lead")), by=LoanID] # Creating a temporary variable for checking if the next observation of the account is in the dataset (exists)
  
  # --- Creating variables for indicating whether the previous/ next record of an account exists
  dat_sub[, Prev_Exist := Counter_Prev==Counter-1]
  dat_sub[is.na(Prev_Exist), Prev_Exist := F] # Checking whether this variable is missing (should then be FALSE)
  
  dat_sub[, Next_Exist := Counter_Next==Counter+1]
  dat_sub[is.na(Next_Exist), Next_Exist := F] # Checking whether this variable is missing (should then be FALSE)
  
  # --- Creating a variable for identifying whether a performance/ default spell should be included in an subsequent analysis
  dat_sub[!is.na(DefSpellID) & !(is.na(PerfSpellID)) & Prev_Exist==F, PerfSpell_F := T] # Identifying all instances where a loan is in default (and performance) and the previous observation doesn't exist
  dat_sub[!is.na(DefSpellID) & !is.na(PerfSpellID) & Next_Exist==F, DefSpell_F := T] # Identifying all instances where a loan is in default (and performance) and the next observation doesn't exist
  
  # - Amending the variables for identifying FALSE performance/ default spells
  dat_sub[is.na(PerfSpell_F), PerfSpell_F := F] # Correcting for all other instances of performance spells (since the FALSE ones have already been identified)
  dat_sub[is.na(DefSpell_F), DefSpell_F := F]
  
  # - Amending the default variable for single observation default spells
  dat_sub[DefSpell_Counter==1 & Date==DefSpell_Max_Date & Prev_Exist==T, DefSpell_F := F] #  Because of the overlap, this one observation will be in both the training and validation dataset. This ensures that only one of the observations is chosen to take into account.
  
  # --- Subsetting dat_sub to only include the required variables (that which is to be returned)
  dat_sub <- subset(dat_sub, select=c("LoanID", "Date", "PerfSpell_F", "DefSpell_F")) %>% setDT(key=c("LoanID", "Date"))
  
  # --- Adding the false performance- and default variables to the given dataset
  dat_given <- cbind(dat_given, dat_sub[, list(PerfSpell_F, DefSpell_F)])
  
  # --- Amending the dataset so that FALSE performance- and default spells do not have an associated spell key, counter, and max date (enables easier subsettin)
  dat_given[PerfSpell_F==T, c((PerfSpellID), (PerfSpell_Counter), (PerfSpell_Max_Date)) := lapply(.SD, function(x) {x=NA}), .SDcols = c((PerfSpellID), (PerfSpell_Counter), (PerfSpell_Max_Date))]
  dat_given[DefSpell_F==T, c((DefSpellID), (DefSpell_Counter), (DefSpell_Max_Date)) := lapply(.SD, function(x) {x=NA}), .SDcols = c((DefSpellID), (DefSpell_Counter), (DefSpell_Max_Date))]

  # --- Returning the dataset
  return(dat_given)

}

# --- Checks
# dat_train1 <- False_Pef_Def(dat_train1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# dat_valid1 <- False_Pef_Def(dat_valid1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# 
# lookup_IDs <- unique(datCredit_real[PerfSpell_Num>=5, LoanID])
# 
# lookup <- datCredit_real[LoanID=="3000002499333", ]
# lookup_t <- dat_train1[LoanID=="3000002499333", ]
# lookup_v <- dat_valid1[LoanID=="3000002499333", ]
# lookup_sub <- dat_sub[LoanID==lookup_IDs[200], ]
# 
# # - Overlaps
# overlaps <- which(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key] %in% dat_valid2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key])
# 
# lookup_IDs <- str_sub(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key][overlaps], start=1, end=-3)
# 
# lookup <- datCredit_real[LoanID==lookup_IDs[10], ]
# lookup_t <- dat_train2[LoanID==lookup_IDs[10], ]
# lookup_v <- dat_valid2[LoanID==lookup_IDs[10], ]


# --- Checks
# dat_train1 <- False_Pef_Def(dat_train1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# dat_valid1 <- False_Pef_Def(dat_valid1, LoanID="LoanID", Date="Date", PerfSpellID="PerfSpell_Key", DefSpellID="DefSpell_Key",
#                             Counter="Counter", PerfSpell_Counter="PerfSpell_Counter", DefSpell_Counter="DefSpell_Counter",
#                             DefSpell_Max_Date="DefSpell_Max_Date", PerfSpell_Max_Date="PerfSpell_Max_Date", DefSpellResol_Type_Hist="DefSpellResol_Type_Hist")
# 
# lookup_IDs <- unique(datCredit_real[PerfSpell_Num>=5, LoanID])
# 
# lookup <- datCredit_real[LoanID=="3000002499333", ]
# lookup_t <- dat_train1[LoanID=="3000002499333", ]
# lookup_v <- dat_valid1[LoanID=="3000002499333", ]
# lookup_sub <- dat_sub[LoanID==lookup_IDs[200], ]
# 
# # - Performance spells
# # The start of default spells overlaps with the end of performance spells. When sub-setting for performance spells, those overlaps are including (although the entire performance spell is not included, only the last observation).
# # Need to identify the observations with no further performance spell history as they are false.
# 
#
# # - Overlaps
# overlaps <- which(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key] %in% dat_valid2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key])
# 
# lookup_IDs <- str_sub(dat_train2[DefSpell_Counter==1 & DefSpell_F==F, DefSpell_Key][overlaps], start=1, end=-3)
# 
# lookup <- datCredit_real[LoanID==lookup_IDs[10], ]
# lookup_t <- dat_train2[LoanID==lookup_IDs[10], ]
# lookup_v <- dat_valid2[LoanID==lookup_IDs[10], ]


