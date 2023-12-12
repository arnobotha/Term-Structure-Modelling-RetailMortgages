# ======================== OPTIMAL SUBSAMPLE SIZES FOR SURVIVAL MODELS USING sCS ========================
# Determining the effect of a wide range of subsample sizes within a simple clustered sampling resampling
# scheme amidst a survival modelling setup (PD-modelling).
# --------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Marcel Muller

# DESCRIPTION:
# This ancillary & exploratory script iterates across a given vector of subsample sizes towards
# resampling data into a basic cross-validation setup (training:validation), using either random sampling
# or n-way stratified sampling. For the latter, a frequency analysis can be conducted to exclude
# observations/ records from strata with sizes below a specified threshold.
# Each chosen subsample size is executed multiple times, using simple clustered sampling for both the
# subsampling and resampling, over various seed values to account for randomness in following a broader
# Monte Carlo setup. Thereafter, various error measures are calculated within each iteration, whereupon
# these error values are appropriately aggregated and graphed into a single cohesive graph.
# --------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup
#   - 2f.Data_Fusion1

# -- Inputs:
#   - datCredit_real | Prepared from script 2d.
#
# -- Outputs:
#   - Difference graph across subsample sizes, where 'difference' is the value of an error measure 
#       between prior probabilities across resampled sets (full:training & training:validation)
# --------------------------------------------------------------------------------------------------------




# ------ 1. Preliminaries
# --- Load in Dataset
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

# --- Some feature engineering
# - Creating a variable for the maximum performance spells per loan (used as stratification variable)
datCredit_real[, PerfSpell_Num_Max := ifelse(all(is.na(PerfSpell_Num)), 0, max(PerfSpell_Num, na.rm=T)), by=LoanID]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date := NA] # Ensure that this performance spell level variable is undefined for default spells.

# - Creating a variable for the first observation of a loan (used as stratification variable)
datCredit_real[, Date_First := Date[1], by=LoanID]
# [SANITY CHECK] Checking if the variable was created correctly
(check.1 <- datCredit_real[is.na(Date_First),.N] == 0) # Should be TRUE
cat(check.1 %?% 'SAFE: variable [Date_First] was successfully created.\n' %:% 'WARNING: variable [Date_First] was not successfully created!\n')

# - Max date of each default spell ( required for curing indicator)
datCredit_real[, DefSpell_Max_Date := max(Date), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), DefSpell_Max_Date := NA] # Ensure that this default spell level variable is undefined for performance spells.

# - Creating a curing indicator for default spells (enables graphing of curing events in default spells)
datCredit_real[, Cured_Ind := ifelse(!is.na(DefSpell_Key) & DefSpellResol_Type_Hist=="Cured" & Date==DefSpell_Max_Date,1,0)] # | Reduce conditions to only the second and third (check feasibility)

# - Loan ever defaulted
datCredit_real[, HasDefault_Ever := ifelse(max(DefSpell_Num)>0, "Has default(s)", "Has no default(s)"), by=list(LoanID)]

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "Settled & Other",
                                                                                  TRUE ~ NA))
datCredit_real[is.na(PerfSpell_Key),.N] == datCredit_real[is.na(PerfSpellResol_Type_Hist2),.N] # TRUE, field created successfully


# - Required field names
targetVar <- c("DefaultStatus1_lead_12_max") # Field name of the main target (i.e., the 12-month default rate)
CurStatus <- "DefaultStatus1" # Field name of the current status of an account (default vs non-default)
resolPerf <- "PerfSpellResol_Type_Hist2" # Field name of performance spell resolution types - first level should be the target event (default)
resolPerf_targetVar <- "Defaulted" # Reference level in the performance spell resolution type (as specified by [resolPerf_start]) for the target variable
clusVar <- "LoanID"
clusVar_Perf <- "PerfSpell_Key"
timeVar <- "Date"
counter <- "Counter"
# - Optional field names
stratifiers <- c("Date_First") # First variable should be of type "date" | Assign "NA" for no stratifiers | Other good stratifier candidates are "Event_Type", "LN_TPE", and HasDefaulted_Ever
# - Final selection
selectionVar <- unique(c(clusVar, clusVar_Perf, timeVar, counter,
                         resolPerf, CurStatus, stratifiers, targetVar)) # Variables to subset
selectionVar <- selectionVar[!is.na(selectionVar)] # Facilitating cases where the variables are left unspecified (specifically for use of no stratifiers)

# --- Subset given dataset accordingly; an efficiency enhancement
datCredit <- subset(datCredit_real, select=selectionVar)
rm(datCredit_real); gc()


# --- Feature engineering spell level date variables
# - Max- and min date of each performance spell (used for computing the prior probabilities)
datCredit[!is.na(get(clusVar_Perf)), c("timeVar_Perf_Min","timeVar_Perf_Max") := as.list(range(get(timeVar), na.rm=TRUE)), by=list(get(clusVar_Perf))]
datCredit[is.na(get(clusVar_Perf)), c("timeVar_Perf_Min","timeVar_Perf_Max") := as.list(c(NA,NA)), by=list(get(clusVar_Perf))]
datCredit[!is.na(get(clusVar_Perf)),.N] == datCredit[!is.na(timeVar_Perf_Max),.N] # TRUE, field created successfully
datCredit[!is.na(get(clusVar_Perf)),.N] == datCredit[!is.na(timeVar_Perf_Min),.N] # TRUE, field created successfully

# --- Calculate prior probability of the (conditional) target event. 
# NOTE: Precalculating the prior probability of the target variable - This is merely a coding optimisation
resolPerf_targetVar <- "Defaulted" # Reference level in the performance spell resolution corresponding to the target variable"
prior_prop <- datCredit[get(timeVar)==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
(prior_prop <- prior_prop[names(prior_prop)[names(prior_prop) == resolPerf_targetVar]][[1]]) # Subsetting for only the target variable
# Prior = 0.1880574

# --- Calculate 12-month conditional default rate on population
# NOTE: Precalculating this is merely a coding optimisation
StartDte <- min(datCredit[,get(timeVar)], na.rm=T) # First date in sampling window
EndDte <- max(datCredit[,get(timeVar)], na.rm=T) # Last date in sampling window
maxDate <- EndDte - years(1) # Post-hoc filter for a conditional 12-month default rate
eventRate_prop <- datCredit[get(CurStatus)==0, list(EventRate = sum(get(targetVar), na.rm=T)/.N),
                            by=timeVar][get(timeVar) >= StartDte & get(timeVar) <= maxDate, ][order(get(timeVar)), ] # Computing the event rate by date
colnames(eventRate_prop) <- c("timeVar", "EventRate") # Renaming columns for easier subsetting
plot(x=eventRate_prop$timeVar, y=eventRate_prop$EventRate, type="b") # Quick visual check

# --- Some parameters
# - Confidence interval parameter
confLevel <- 0.95
# - General parameters
cpu.threads <- 6
confLevel <- 0.95
# - Iteration parameters
smp_size_v <- c(seq(from=2500, to=25000, by=2500), seq(from=30000, to=75000, by=5000),
                seq(from=82500, to=150000, by=7500), seq(from=165000, to=255000, by=15000),
                seq(from=285000, to=465000, by=30000), seq(from=505000, to=655000, by=50000))
seed_v <- c(1:100)
# - Sampling fraction
smp_frac <- 0.7




# ------ 2. Subsampled resampling scheme: basic cross-validation with random sampling

# --- Define function for applying a subsampled resampling scheme given parameters on given data
# This function serves as an "outer job" to be called within a multithreaded environment
# - Inputs: smp_size: Subsample size; smp_frac: sampling fraction for resmpling scheme;
# stratifiers: vector of stratification field names for n-way stratified sampling inner technique;
# targetVar: outcome field name within performance spells;
# clusVar: unique account ID name
# CurStatus: Current default status of an acount (used to compute the conditional 12-month default rate for performance spells)
# Counter: name of field counting the sequential observation of an account
# timeVar: field name of date for event rate calculations;
# seed: specific seed value
# prior_prop: pre-calculated prior probability within population for error measurement
# eventRate_prop: pre-calculated performance spell event rates over time within population for error measurement
# minStrata_size: Minimum strata size which should be enforced in the resulting subsample
# datGiven: given dataset from which to subsample and resample
# EventRate: indicator for computing the event rate as the second error measure
subSmp_strat <- function(smp_size, smp_frac, stratifiers=NA, targetVar=NA, clusVar=NA, CurStatus=NA, Counter=NA, timeVar=NA, timeVar_Min=NA,
                         seed=123, prior_prop=NA, eventRate_prop=NA, minStrata_size=0, datGiven, EventRate=F) {
  # --- UNIT TEST
  # datGiven <- copy(datCredit); smp_size <- 5000; seed <- 1; minStrata_size <- 0; timeVar_Min<-"timeVar_Perf_Min"; EventRate <- T
  
  # --- Preliminaries
  # - Error Checks
  if (any(is.na(stratifiers)) & is.na(targetVar)) { stop("Stratifiers and target variable are unspecified! Must at least include the target variable")}
  # - Implied sampling fraction for downsampling step
  smp_perc <- smp_size/length(unique(datGiven[,get(clusVar)]))
  
  
  # --- Downsample data into a set with a fixed size (using stratified sampling) before implementing resampling scheme
  # - Set seed
  set.seed(seed)
  # - Conditional loop for stratifiers
  if (all(is.na(stratifiers))){ # - No stratifiers
    # Get unique loan account IDs from the full dataset
    dat_keys <- unique(datGiven[, mget(c(clusVar))])
    # Use simple random sampling to select the loan IDs that ought to be in the subsampled dataset
    dat_smp_keys <- dat_keys %>% slice_sample(prop=smp_perc) %>% as.data.table()
  } else { # - Stratifiers
    # Get unique loan account IDs from the full dataset
    dat_keys <- unique(datGiven[, mget(c(clusVar, stratifiers))])
    # Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the subsampled dataset
    dat_smp_keys <- dat_keys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
  }
  # - Obtain the associated loan records as to create the subsampled dataset
  datGiven_smp <- copy(datGiven[get(clusVar) %in% dat_smp_keys[, get(clusVar),]])
  
  # --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold
  if (all(!is.na(stratifiers))){ # - Conditional loop for strata
    selectionVar_smp <- c(clusVar, timeVar, stratifiers)
    datStrata_smp_min <- datGiven_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
    
    # - Conditionally applying the exclusions
    if (sum(datStrata_smp_min[,Freq]) > 0){
      # Initiating a vector which will contain the exclusion IDs
      dat_keys_exc <- NA
      # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
      for (i in 1:datStrata_smp_min[,.N]){
        class_type <- sapply(datStrata_smp_min[,1:length(stratifiers)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
        
        excCond <- datStrata_smp_min[i,1:length(stratifiers)] # Getting the values of the ith minimum strata
        excCond <- data.table(Stratifier = colnames(excCond), # Building a dataset
                              Value = unname(t(excCond)), # Ensure that the column name is Value instead of Value.V1
                              Class = class_type)
        excCond[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
        excCond[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
        excCond2 <- parse(text = paste0(paste0(excCond$Condition, collapse = ""), counter,"==1")) # Compiling the ith condition
        
        dat_keys_exc <- c(dat_keys_exc, as.vector(datGiven_smp[eval(excCond2), get(clusVar)]))
      }
      dat_keys_exc <- dat_keys_exc[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
      
      # Applying the exclusions to the subsampled dataset
      datGiven_smp <- copy(datGiven_smp[!(get(clusVar) %in% dat_keys_exc),])
      
      # - Obtaining the stratum that are below the minimum
      # datStrata_smp_min <- datGiven_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
      # cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datGiven_smp[get(counter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datGiven_smp[get(counter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
    }
  }  
  
  
  # --- Implementing the resampling scheme
  # - Set seed
  set.seed(1)
  # - Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the training dataset
  if (all(!is.na(stratifiers))){
    dat_train_keys <- dat_smp_keys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_frac) %>% as.data.table() 
  } else {
    dat_train_keys <- dat_smp_keys %>% slice_sample(prop=smp_frac) %>% as.data.table()
  }
  # - Obtain the associated loan records as to create the training dataset
  datGiven_train <- copy(datGiven_smp[get(clusVar) %in% dat_train_keys[, get(clusVar)],])
  # - Obtain the associated loan records of the validation dataset
  datGiven_valid <- copy(datGiven_smp[!(get(clusVar) %in% dat_train_keys[, get(clusVar)]),])
  # - [SANITY CHECK]
  # (check.2 <- datGiven_smp[,.N] == datGiven_train[,.N] + datGiven_valid[,.N]) # Should be TRUE

  
  # --- Create the output table
  if (all(!is.na(stratifiers))){ # Stratifiers
    datTemp <- data.table("SubSample_Size"=datGiven_smp[get(counter)==1,.N], "SubSample_FullSize"=datGiven_smp[,.N], "Exclusions"=sum(datStrata_smp_min[,Freq]))
  } else { # No stratifiers
    datTemp <- data.table("SubSample_Size"=datGiven_smp[get(counter)==1,.N], "SubSample_FullSize"=datGiven_smp[,.N])
  }
  
  # --- Calculate error measure 1: Difference in prior probabilities between population and training (as subsampled + resampled)
  # Calculate prior probabilities within each relevant dataset, e.g., proportion of defaults across all time
  
  # - Population
  if (is.na(prior_prop)) {
    prior_prop <- datGiven[get(timeVar)==get(timeVar_Min), get(resolPerf), with=T] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
    prior_prop <- prior_prop[names(prior_prop)[names(prior_prop) == resolPerf_targetVar]][[1]] # Subsetting for only the target variable
  }
  
  # - Subsampled & resampled training set
  prior_prop_train <- datGiven_train[get(timeVar)==get(timeVar_Min), get(resolPerf)] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
  prior_prop_train <- prior_prop_train[names(prior_prop_train)[names(prior_prop_train) == resolPerf_targetVar]][[1]] # Subsetting for only the target variable
  
  # - Compare population with training set using chosen error measure
  err_priorProb_AE <- abs(prior_prop - prior_prop_train) # absolute error
  err_priorProb_SqrdErr <- (prior_prop - prior_prop_train)^2 # squared error
  
  # - Append output table (preliminary)
  datTemp <- data.table(datTemp, "SampleFrac"=smp_frac, "Stratifiers"=paste(stratifiers, collapse="; "),
                        "Seed"=seed, "Err_PriorProb_AE" = err_priorProb_AE, "Err_PriorProb_SqrdErr" = err_priorProb_SqrdErr)
  
  
  # --- Calculate error measure 2: MAE between 2 time series of the event rate between population/training and training/validation (as subsampled + resampled)
  # NOTE: performance spell event rate is a 12-month conditional event rate, e.g., k-month default rate at t+k given that event has not happened at t
  # This is an optional error measure
  if (EventRate) {
    # - Population
    if (any(is.na(eventRate_prop))) {
      StartDte <- min(datGiven[,get(timeVar)], na.rm=T) # First date in sampling window
      EndDte <- max(datGiven[,get(timeVar)], na.rm=T) # Last date in sampling window
      maxDate <- EndDte - years(1) # Post-hoc filter for a conditional 12-month default rate
      eventRate_prop <- datGiven[get(CurStatus)==0, list(EventRate = sum(get(targetVar), na.rm=T)/.N),
                                  by=timeVar][get(timeVar) >= StartDte & get(timeVar) <= maxDate, ][order(get(timeVar)), ] # Computing the event rate by date
      colnames(eventRate_prop) <- c("timeVar", "EventRate") # Renaming columns for easier subsetting
    }
    
    # - Subsampled-resampled training set
    StartDte <- min(datGiven_train[,get(timeVar)], na.rm=T) # First date in sampling window
    EndDte <- max(datGiven_train[,get(timeVar)], na.rm=T) # Last date in sampling window
    maxDate <- EndDte - years(1) # Post-hoc filter for a conditional 12-month default rate
    eventRate_prop_train <- datGiven_train[get(CurStatus)==0, list(EventRate = sum(get(targetVar), na.rm=T)/.N),
                                            by=timeVar][get(timeVar) >= StartDte & get(timeVar) <= maxDate, ][order(get(timeVar)), ] # Computing the event rate by date
    colnames(eventRate_prop_train) <- c("timeVar", "EventRate") # Renaming columns for easier subsetting
    
    # - Subsampled-resampled validation set
    eventRate_prop_valid <- datGiven_valid[get(CurStatus)==0, list(EventRate = sum(get(targetVar), na.rm=T)/.N),
                                           by=timeVar][get(timeVar) >= StartDte & get(timeVar) <= maxDate, ][order(get(timeVar)), ] # Computing the event rate by date
    colnames(eventRate_prop_valid) <- c("timeVar", "EventRate") # Renaming columns for easier subsetting
    
    # - Compare event rates across different sets using chosen error measure
    err_eventRate_MAE_train <- mean(abs(eventRate_prop$EventRate - eventRate_prop_train$EventRate), na.rm=T) # mean absolute error
    err_eventRate_MAE_valid <- mean(abs(eventRate_prop$EventRate - eventRate_prop_valid$EventRate), na.rm=T) # mean absolute error
    err_eventRate_MAE_trainvalid <- mean(abs(eventRate_prop_train$EventRate - eventRate_prop_valid$EventRate), na.rm=T) # mean absolute error
    
    # - Append error value to output table
    datTemp <- data.table(datTemp, "Err_EventRate_PopTrain_MAE" = err_eventRate_MAE_train, 
                          "Err_EventRate_PopValid_MAE" = err_eventRate_MAE_valid,
                          "Err_EventRate_TrainValid_MAE" = err_eventRate_MAE_trainvalid)
  }
  
  
  # --- Return value of chosen error measure
  return(datTemp)
} # end of function


# - Testing function call
ptm <- proc.time() #IGNORE: for computation time calculation
test <- subSmp_strat(smp_size=15000, smp_frac=smp_frac, stratifiers=stratifiers, targetVar=targetVar, clusVar=clusVar, CurStatus=CurStatus, Counter=counter,
                     timeVar=timeVar, timeVar_Min="timeVar_Perf_Min", seed=1, prior_prop=NA, eventRate_prop=NA, minStrata_size=0, datGiven=datCredit, EventRate=T)
proc.time() - ptm  #IGNORE: for computation time calculation
# Stratification with smp_size=150000 and seed=1:
#     Optimal Sampling : SubSample_Size = 14905; Exclusions = 0; Err_EventRate_PopTrain_MAE = 0.001966296; Err_EventRate_PopValid_MAE = 0.003110211; Err_EventRate_TrainValid_MAE = 0.003564881
#     Resampling Tool  : SubSample_Size = 14905 ; Exclusions = 0; Err_EventRate_PopTrain_MAE = 0.001966; Err_EventRate_PopValid_MAE = 0.003110; Err_EventRate_TrainValid_MAE = 0.003565
# No stratification with smp_size=150000 and seed=1:
#     Optimal Sampling : SubSample_Size = 14999; Exclusions = 0; Err_EventRate_PopTrain_MAE = 0.001869999; Err_EventRate_PopValid_MAE = 0.003527054; Err_EventRate_TrainValid_MAE = 0.003610729
#     Resampling Tool  : SubSample_Size = 64858 ; Exclusions = 0; Err_EventRate_PopTrain_MAE = 0.001869999; Err_EventRate_PopValid_MAE = 0.003527054; Err_EventRate_TrainValid_MAE = 0.003610729

# --- Cleanup
rm(test, ptm)


# --- Main Loop (outer function call)
cl.port <- makeCluster(cpu.threads-1)
registerDoParallel(cl.port)

cat(paste0("1 (", Sys.time(),"). Iterating across subsample sizes ..."),
    file="subsampleLoop", append=T)

ptm <- proc.time() #IGNORE: for computation time calculation

# - Multithreaded looping procedure using the foreach-package
datResults <- foreach(it=1:(length(seed_v)*length(smp_size_v)), .combine='rbind', .verbose=F, .inorder=T, 
                      .packages=c("dplyr","data.table","tidyselect","lubridate", "scales")) %dopar%
  {
    # - Testing 
    # it <- 3
    
    # - Set indices
    iSeed <- (it-1) %% length(seed_v) + 1 # modulo operation
    iSize <- (it-1) %/% length(seed_v) + 1 # integer-valued division
    
    # - Iterate 
    temp <- subSmp_strat(smp_size=smp_size_v[iSize], smp_frac=smp_frac, stratifiers=stratifiers, targetVar=targetVar, clusVar=clusVar, CurStatus=CurStatus, Counter=counter,
                         timeVar=timeVar, timeVar_Min="timeVar_Perf_Min", seed=seed_v[iSeed], prior_prop=prior_prop, eventRate_prop=eventRate_prop, minStrata_size=0, datGiven=datCredit, EventRate=T)
    
    # - Reporting
    if (iSeed == length(seed_v)) {
      cat(paste0("\n2 (", Sys.time(),"). Subsample size: ", comma(smp_size_v[iSize]), " tested ",length(seed_v), " times."),
          file="subsampleLoop.txt", append=T) 
    }
    
    return(temp)
  }  

t <- proc.time() - ptm  #IGNORE: for computation time calculation
cat(paste0("\n3 (", Sys.time(),"). ForEach-loop done. Elapsed time: ", round(t[3]/60), " minutes."),
    file="subsampleLoop.txt", append=T)

# - Save to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genObjPath, "subSampleSizes_alt"), datResults); gc()
stopCluster(cl.port)




# ------ 3. Graphing
# --- Load in Dataset
if (!exists('datResults')) unpack.ffdf(paste0(genObjPath,"subSampleSizes_alt"), tempPath)


# --- Aggregate to subsample size level
datGraph <- datResults[, list(PriorProb_MAE = mean(Err_PriorProb_AE , na.rm=T), PriorProb_MAE_SD = sd(Err_PriorProb_AE , na.rm=T),
                              PriorProb_RMSE = sqrt(sum(Err_PriorProb_SqrdErr, na.rm=T)/.N), PriorProb_RMSE_SE = sd(Err_PriorProb_SqrdErr, na.rm=T),
                              EventRate_PopTrain_MAE_Mean = mean(Err_EventRate_PopTrain_MAE, na.rm=T), EventRate_PopTrain_MAE_SD = sd(Err_EventRate_PopTrain_MAE, na.rm=T),
                              EventRate_TrainValid_MAE_Mean = mean(Err_EventRate_TrainValid_MAE, na.rm=T), EventRate_TrainValid_MAE_SD = sd(Err_EventRate_TrainValid_MAE, na.rm=T),
                              SubSample_FullSize_Mean = round(mean(SubSample_FullSize, na.rm=T)), N=.N),
                       by=list(SubSample_Size)]

# --- Create 95% confidence interval for point estimate (mean) : Population-training set comparison
datGraph[, EventRate_MAE_PopTrain_Mean_ErrMargin := (qnorm(1-(1-confLevel)/2)*EventRate_PopTrain_MAE_SD/sqrt(N))]
datGraph[, EventRate_MAE_PopTrain_Mean_lower := EventRate_PopTrain_MAE_Mean - EventRate_MAE_PopTrain_Mean_ErrMargin]
datGraph[, EventRate_MAE_PopTrain_Mean_upper := EventRate_PopTrain_MAE_Mean + EventRate_MAE_PopTrain_Mean_ErrMargin]


# --- Create 95% confidence interval for point estimate (mean) : Training-Validation set comparison
datGraph[, EventRate_MAE_TrainValid_Mean_ErrMargin := (qnorm(1-(1-confLevel)/2)*EventRate_TrainValid_MAE_SD/sqrt(N))]
datGraph[, EventRate_MAE_TrainValid_Mean_lower := EventRate_TrainValid_MAE_Mean - EventRate_MAE_TrainValid_Mean_ErrMargin]
datGraph[, EventRate_MAE_TrainValid_Mean_upper := EventRate_TrainValid_MAE_Mean + EventRate_MAE_TrainValid_Mean_ErrMargin]


# --- Create summary table for annotations within graph - Performance Spells
# - Rounding the sample sizes to the nearest 50000 to help aesthetics
datGraph[,SubSample_FullSize_Mean := plyr::round_any(SubSample_FullSize_Mean, 50000)]
datAnnotate <- datGraph[SubSample_FullSize_Mean %in% c(200000, 550000, 900000, 1300000, 1650000, 2200000, 2950000, 3650000, 4400000, 6600000, 12100000, 20900000),
                        list(`"Size "*italic(s)`=comma(SubSample_FullSize_Mean), 
                               `italic(E)(epsilon(italic(s)))*" for "*italic(D):italic(D[T])`=paste0(sprintf("%.3f", EventRate_PopTrain_MAE_Mean*100),"%"),
                               `95% CI`=paste0("± ",sprintf("%.4f",EventRate_MAE_PopTrain_Mean_ErrMargin*100),"%"),
                               `italic(E)(epsilon(italic(s)))*" for "*italic(D[T]):italic(D[V])`=paste0(sprintf("%.3f", EventRate_TrainValid_MAE_Mean*100),"%"),
                               `95% CI`=paste0("± ",sprintf("%.4f",EventRate_MAE_TrainValid_Mean_ErrMargin*100),"%"))]


# - SCRATCH
plot(x=datGraph$SubSample_FullSize_Mean, y=datGraph$PriorProb_MAE, type="b")
plot(x=datGraph$SubSample_FullSize_Mean, y=datGraph$PriorProb_RMSE, type="b")
plot(x=datGraph$SubSample_FullSize_Mean, y=datGraph$EventRate_PopTrain_MAE_Mean, type="b")
plot(x=datGraph$SubSample_FullSize_Mean, y=datGraph$EventRate_TrainValid_MAE_Mean, type="b")
plot(x=datGraph$SubSample_FullSize_Mean, y=datGraph$EventRate_MAE_PopTrain_Mean_ErrMargin, type="b")


# --- Pivot for graphing purposes: 2 different set comparisons using single error measure
datGraph2 <- pivot_longer(datGraph[,list(SubSample_FullSize_Mean, a_EventRate_PopTrain=EventRate_PopTrain_MAE_Mean, b_EventRate_TrainValid=EventRate_TrainValid_MAE_Mean)],
                               cols=a_EventRate_PopTrain:b_EventRate_TrainValid, names_to = "Set", values_to = "Value") %>% as.data.table()

datGraph2_lower <- pivot_longer(datGraph[,list(SubSample_FullSize_Mean, a_EventRate_PopTrain=EventRate_MAE_PopTrain_Mean_lower, b_EventRate_TrainValid=EventRate_MAE_TrainValid_Mean_lower)],
                                cols=a_EventRate_PopTrain:b_EventRate_TrainValid, names_to = "Set", values_to = "Value_Lower") %>% as.data.table()

datGraph2_upper <- pivot_longer(datGraph[,list(SubSample_FullSize_Mean, a_EventRate_PopTrain=EventRate_MAE_PopTrain_Mean_upper, b_EventRate_TrainValid=EventRate_MAE_TrainValid_Mean_upper)],
                                cols=a_EventRate_PopTrain:b_EventRate_TrainValid, names_to = "Set", values_to = "Value_Upper") %>% as.data.table()

datGraph2_margins <- merge(datGraph2_lower, datGraph2_upper, by=c("SubSample_FullSize_Mean", "Set"))

datGraph3 <- merge(datGraph2, datGraph2_margins, by=c("SubSample_FullSize_Mean", "Set"))


# --- Find elbow point where differential between subsequent error values becomes negligible
### I.e., find x where 2nd derivative of f(x) is near zero
# - Compute gradients
datEventRate_PopTrain_1st <- datGraph[order(SubSample_Size), list(Gradient = diff(EventRate_PopTrain_MAE_Mean) / diff(SubSample_Size))]$Gradient
datEventRate_TrainValid_1st <- datGraph[order(SubSample_Size), list(Gradient = diff(EventRate_TrainValid_MAE_Mean) / diff(SubSample_Size))]$Gradient
# plot(datEventRate_PopTrain_1st, type="b", main="1st derivative") # 1st derivative neither smooth nor monotonic
# plot(datEventRate_TrainValid_1st, type="b", main="1st derivative") # 1st derivative neither smooth nor monotonic
# plot(diff(datEventRate_PopTrain_1st), type="b", main="2nd derivative") # 2nd derivative neither smooth nor monotonic
# plot(diff(datEventRate_TrainValid_1st), type="b", main="2nd derivative") # 2nd derivative is smooth but not monotonic
# - Find index of stationary points x such that f''(x) <= epislon
statPoint_PopTrain <- which(abs(diff(datEventRate_PopTrain_1st)) < 10^(-10.25))[1] # Matter investigated; there is "jaggedness" in the plot of the second derivative, and that is why we take the second point
statPoint_TrainValid <- which(abs(diff(datEventRate_TrainValid_1st)) < 10^(-10.5))[1] # 2nd derivative is neither smooth nor monotonic
# statPoint_TrainValid <- which(abs(datEventRate_TrainValid_1st) < 10^(-9))[1]
# - Find corresponding x-values at index
SubSamp_PopTrain <- datGraph$SubSample_FullSize_Mean[statPoint_PopTrain+1] # (x + 1 to correspond to indices of original vector)
SubSamp_TrainValid <- datGraph$SubSample_FullSize_Mean[statPoint_TrainValid+1] # (x + 1 to correspond to indices of original vector)
# - Find corresponding y-values at index
Err_PopTrain <- datGraph$EventRate_PopTrain_MAE_Mean[statPoint_PopTrain+1] # (x + 1 to correspond to indices of original vector)
Err_TrainValid <- datGraph$EventRate_TrainValid_MAE_Mean[statPoint_TrainValid+1] # (x + 1 to correspond to indices of original vector)

# - Create elbow annotation object
datAnnotate_elbow <- data.table(Set=c("a_EventRate_PopTrain", "b_EventRate_TrainValid"), 
                                SubSample_Size=c(SubSamp_PopTrain, SubSamp_TrainValid), Value=c(Err_PopTrain,Err_TrainValid),
                                Label=paste0( c(comma(SubSamp_PopTrain/1000000), comma(SubSamp_TrainValid/1000000)), "m") )


# --- Graphing
# - Aesthetic engineering
datGraph3[, Facet_label := factor("'Comparison of 12-month default rate series across sets: '*(italic(D):italic(D[T]))*' ; '*(italic(D[T]):italic(D[V]))")]

# - Graphing parameters
chosenFont <- "Cambria"; dpi <- 180
col.v <- brewer.pal(10, "Paired")[c(10,8)]
fill.v <- brewer.pal(10, "Paired")[c(9,7)]
linetype.v <- c("solid", "dotted")
label.v <- list(expression(italic(D)*" vs "*italic(D[T])),
                expression(italic(D[T])*" vs "*italic(D[V])) )
label.v2 <- list(expression(italic(D)*" vs "*italic(D[T])),
                 expression(italic(D[T])*" vs "*italic(D[V])) )

# --- Create main graph for performance spells
(g1 <- ggplot(datGraph3[SubSample_FullSize_Mean<=25000000,], aes(x=SubSample_FullSize_Mean, y=Value, group=Set)) + theme_minimal() + 
    labs(x=bquote("Subsample size "*italic(s)*" = |"*italic(D[S])*"|"), y=bquote("Error measure value "*italic(E)*'['*epsilon(italic(s))*"] (%)")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # annotate elbow/stationary point beyond which error measure decreases markedly slower
    geom_point(data=datAnnotate_elbow, aes(x=SubSample_Size, y=Value, colour=Set), shape=1, size=5, show.legend=F) +
    geom_segment(data=datAnnotate_elbow, aes(x=SubSample_Size, xend=SubSample_Size, y=0, yend=Value, colour=Set),
                 show.legend=F, linewidth=0.3, linetype="dashed") + 
    geom_label(data=datAnnotate_elbow, aes(x=SubSample_Size, label=Label, y=Value, colour=Set), show.legend = F, 
               nudge_y=0.0002, nudge_x=1450000, size=3) + 
    # main line graph with overlaid points
    geom_ribbon(aes(x=SubSample_FullSize_Mean, ymin=Value_Lower, ymax=Value_Upper, fill=Set), alpha=0.5) + 
    geom_line(aes(colour=Set, linetype=Set), linewidth=0.5) + 
    geom_point(aes(x=SubSample_FullSize_Mean, y=Value, colour=Set, shape=Set), size=1.3) + 
    # annotate data table
    annotate(geom="table", x=4000000, y=0.011, family=chosenFont, size=2.9,
            label=datAnnotate, parse=T) +
    # facets & scale options
    facet_grid(Facet_label ~ ., labeller=label_parsed) + 
    scale_colour_manual(name="Mean MAE", values=col.v, label=label.v) + 
    scale_shape_discrete(name="Mean MAE", label=label.v) + 
    scale_linetype_manual(name="Mean MAE", values=linetype.v, label=label.v) + 
    scale_fill_manual(name="95% CI for mean", values=fill.v, label=label.v2) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_continuous(breaks=pretty_breaks(), label=label_comma(scale=0.000001, suffix="m"))
)
# - Save graph
ggsave(g1, file=paste0(genFigPath_Res, "DefaultRates_SubSampleRates_Experiment.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


# --- Cleanup
rm(datCredit, g1); gc()




