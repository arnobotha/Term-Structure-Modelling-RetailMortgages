# ======================== OPTIMAL SUBSAMPLE SIZES FOR SURVIVAL MODELS USING sCS ========================
# Determining the effect of a wide range of subsample sizes within a simple clustered & 1-way stratified 
# resampling scheme amidst a survival modelling setup (PD-modelling).
# --------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB), Marcel Muller (MM)

# DESCRIPTION:
# This ancillary & exploratory script iterates across a given vector of subsample sizes towards
# resampling data into a basic cross-validation setup (training:validation).
# Each chosen subsample size is executed multiple times within a Monte Carlo setup, using simple 
# clustered sampling for both the subsampling and resampling schemes. 
# 1-way stratification on last spell date [Date_First] is necessarily forced to ensure the 
# eventual time series of spell resolution rates can be created over time
# Thereafter, various error measures are calculated within each iteration, whereupon
# these error values are appropriately aggregated and graphed into a single cohesive graph.
# --------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R

# -- Inputs:
#   - datCredit_real | Prepared from script 2d.
#
# -- Outputs:
#   - Difference graph across subsample sizes, where 'difference' is the value of an error measure 
#       between prior probabilities across resampled sets (full:training & training:validation)
# --------------------------------------------------------------------------------------------------------




# ------ 1. Preliminaries
# --- Load in Dataset
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)


# --- Binned resolution type creation and field specification

# - Creating new spell resolution types for analytical purposes
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "Settled & Other",
                                                                                  TRUE ~ NA))
# [SANITY CHECKS]: Any missing spells keys or mismatch between spell keys and associated resolution types?
datCredit_real[is.na(PerfSpell_Key),.N] == datCredit_real[is.na(PerfSpellResol_Type_Hist2),.N] # TRUE, field created successfully


# - Required field names
clusVar <- "LoanID" # for clustered random sampling on subject-level
clusVar_Spell <- "PerfSpell_Key" # for clustered random sampling on spell-level within a subject's history
timeVar <- "Date" # calendar time variable over which main event is observed & predicted, used in strata analysis
counter <- "PerfSpell_Counter" # subject-level row index variable, used in strata analysis
# - Facet specification field names (for graphing purposes of the resolution rates)
resolType <- "PerfSpellResol_Type_Hist2"  # Performance spell resolution types
resolType_Val <- "Defaulted" # Reference value in the performance spell resolution type field [resolType]

# - Final selection
selectionVar <- unique(c(clusVar, clusVar_Spell, timeVar, counter,resolType)) # Variables to subset
selectionVar <- selectionVar[!is.na(selectionVar)] # Facilitating cases where the variables are left unspecified (specifically for use of no stratifiers)


# --- Subset given dataset accordingly; an efficiency enhancement
datCredit <- subset(datCredit_real, select=selectionVar)
rm(datCredit_real); gc()


# --- Feature engineering spell level date variables

# - Max- and min date of each performance spell (used for computing the prior probabilities)
datCredit[!is.na(get(clusVar_Spell)), c("timeVar_SpellMin","timeVar_SpellMax") := as.list(range(get(timeVar), na.rm=TRUE)), by=list(get(clusVar_Spell))]
datCredit[is.na(get(clusVar_Spell)), c("timeVar_SpellMin","timeVar_SpellMax") := as.list(c(NA,NA)), by=list(get(clusVar_Spell))]
# [SANITY CHECKS]: Any missing spell dates or mismatch between spell dates and associated resolution types?
datCredit[!is.na(get(clusVar_Spell)),.N] == datCredit[!is.na(timeVar_SpellMax),.N] # TRUE, field created successfully
datCredit[!is.na(get(clusVar_Spell)),.N] == datCredit[!is.na(timeVar_SpellMin),.N] # TRUE, field created successfully

# - Creating a variable for the first observation date of a loan (used as a necessary stratification variable)
# This is a necessary stratification variable for strata analysis over time, merely graphing purposes therefore
datCredit[, Date_First := get(timeVar)[1], by=get(clusVar)]
# [SANITY CHECK] Checking if the variable was created correctly
(check.1 <- datCredit[is.na(Date_First),.N] == 0) # Should be TRUE
cat(check.1 %?% 'SAFE: variable [Date_First] was successfully created.\n' %:% 'WARNING: variable [Date_First] was not successfully created!\n')


# - Rename variable fields into fixed names for easier reference throughout script
colnames(datCredit) <- c("clusVar", "clusVar_Spell", "timeVar", "counter", "resolType", 
                         "timeVar_SpellMin", "timeVar_SpellMax", "Date_First")



# --- Calculate prior probability of the (conditional) target event. 
# NOTE: Precalculating the prior probability of the target variable - This is merely a coding optimisation
prior_prop <- datCredit[timeVar==timeVar_SpellMin, resolType] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
(prior_prop <- prior_prop[names(prior_prop)[names(prior_prop) == resolType_Val]][[1]]) # Subsetting for only the target variable
# Prior = 0.1880574

# --- Calculate Spell resolution rates over time on population | Spell stop time (t_s)
# NOTE: Precalculating this is merely a coding optimisation

# - set dates
StartDte <- min(datCredit[,timeVar], na.rm=T) # First date in sampling window
EndDte <- max(datCredit[,timeVar], na.rm=T) # Last date in sampling window
maxDate <- EndDte %m+% months(-1) # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
minDate <- StartDte # %m+% month(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window

# - Fixing to spell stop-time, we aggregate to monthly level and observe the time series up to given point
datAggr_cohorts <- merge(datCredit[timeVar==timeVar_SpellMax, list(Sum_Total = .N), by=list(timeVar)],
                         datCredit[timeVar==timeVar_SpellMax, list(Sum_Resol = .N), by=list(timeVar,resolType)],
                         by="timeVar")[timeVar >= minDate & timeVar <= maxDate,]
datAggr_cohorts[, EventRate := Sum_Resol/Sum_Total]

# - Subset only to main resolution rate amongst subtypes (kappa)
eventRate_prop <- subset(datAggr_cohorts, resolType == resolType_Val)
plot(x=eventRate_prop$timeVar, y=eventRate_prop$EventRate, type="b") # Quick visual check

# --- Some parameters
smp_frac <- 0.7 # Sampling fraction
confLevel <- 0.95 # Confidence interval parameter
cpu.threads <- 10
# - Iteration parameters
smp_size_v <- c(seq(from=2500, to=25000, by=2500), seq(from=30000, to=75000, by=5000),
                seq(from=82500, to=150000, by=7500), seq(from=165000, to=255000, by=15000),
                seq(from=285000, to=465000, by=30000), seq(from=505000, to=655000, by=50000))
seed_v <- c(1:100) # vector of seeds (also number of iterations per sample size)





# ------ 2. Subsampled resampling scheme: basic cross-validation with random sampling

# --- Define function for applying a subsampled resampling scheme given parameters on given data
# This function serves as an "outer job" to be called within a multithreaded environment
# - Inputs: smp_size: Subsample size; smp_frac: sampling fraction for resmpling scheme;
# seed: specific seed value
# prior_prop: pre-calculated prior probability within population for error measurement
# eventRate_prop: pre-calculated performance spell event rates over time within population for error measurement
# datGiven: given dataset from which to subsample and resample
# EventRate: indicator for computing the event rate as the second error measure
# sDateMax_Adj: Number of months by which to adjust the maximum of the eventual rate time series when calculating the MAE
subSmp_strat <- function(smp_size, smp_frac, seed=123, prior_prop=NA, eventRate_prop=NA, datGiven, EventRate=F,
                         sDateMax_Adj=-1) {
  # --- UNIT TEST
  # datGiven <- copy(datCredit); smp_size <- 90000; targetVar <- "resolType"
  # seed <- 1; minStrata_size <- 0; EventRate <- T; sDateMax_Adj=-1
  ptm <- proc.time()
  
  # --- Preliminaries
  # - Implied sampling fraction for downsampling step
  smp_perc <- smp_size/length(unique(datGiven[!is.na(clusVar_Spell), clusVar]))
  
  
  # --- Downsample data into a set with a fixed size (using stratified sampling) before implementing resampling scheme
  # - Set seed
  set.seed(seed, kind="Mersenne-Twister")
  # Get unique loan account IDs from the full dataset, with spell info
  dat_keys <- datGiven[!is.na(clusVar_Spell), list(Date_First = Date_First[1]), by=list(clusVar)]
  # Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the subsampled dataset
  dat_smp_keys <- dat_keys %>% group_by(Date_First) %>% slice_sample(prop=smp_perc) %>% as.data.table()
  # - Obtain the associated loan records as to create the subsampled dataset
  datGiven_smp <- copy(datGiven[clusVar %in% dat_smp_keys[, clusVar,]])
  
  
  # --- Implementing the resampling scheme
  # - Set seed
  set.seed(1, kind="Mersenne-Twister")
  # - Use simple random sampling to select the loan IDs that ought to be in the training dataset
  dat_train_keys <- dat_smp_keys %>% group_by(Date_First) %>% slice_sample(prop=smp_frac) %>% as.data.table()
  # - Obtain the associated loan records as to create the training dataset
  datGiven_train <- copy(datGiven_smp[clusVar %in% dat_train_keys[, clusVar],])
  # - Obtain the associated loan records of the validation dataset
  datGiven_valid <- copy(datGiven_smp[!(clusVar %in% dat_train_keys[, clusVar]),])
  # - [SANITY CHECK]
  # (check.2 <- datGiven_smp[,.N] == datGiven_train[,.N] + datGiven_valid[,.N]) # Should be TRUE

  
  # --- Create the output table
  datTemp <- data.table("SubSample_Size"=length(unique(datGiven_smp$clusVar)), "SubSample_FullSize"=datGiven_smp[,.N])
  
  # --- Calculate error measure 1: Difference in prior probabilities between population and training (as subsampled + resampled)
  # Calculate prior probabilities within each relevant dataset, e.g., proportion of defaults across all time
  
  # - Population
  if (is.na(prior_prop)) {
    prior_prop <- datGiven[timeVar==timeVar_SpellMin, resolType] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
    prior_prop <- prior_prop[names(prior_prop)[names(prior_prop) == resolType_Val]][[1]] # Subsetting for only the target variable
  }
  
  # - Subsampled & resampled training set
  prior_prop_train <- datGiven_train[timeVar==timeVar_SpellMin, resolType] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
  prior_prop_train <- prior_prop_train[names(prior_prop_train)[names(prior_prop_train) == resolType_Val]][[1]] # Subsetting for only the target variable
  
  # - Compare population with training set using chosen error measure
  err_priorProb_AE <- abs(prior_prop - prior_prop_train) # absolute error
  err_priorProb_SqrdErr <- (prior_prop - prior_prop_train)^2 # squared error
  
  # - Append output table (preliminary)
  datTemp <- data.table(datTemp, "SampleFrac"=smp_frac, "Seed"=seed, "Err_PriorProb_AE" = err_priorProb_AE, 
                        "Err_PriorProb_SqrdErr" = err_priorProb_SqrdErr)
  
  
  # --- Calculate error measure 2: MAE between 2 time series of the resolution rate between population/training and training/validation (as subsampled + resampled)
  # NOTE: This is an optional error measure
  if (EventRate) {
    # - Population
    if (any(is.na(eventRate_prop))) {
      # - set dates
      StartDte <- min(datGiven[,timeVar], na.rm=T) # First date in sampling window
      EndDte <- max(datGiven[,timeVar], na.rm=T) # Last date in sampling window
      maxDate <- EndDte  %m+% months(sDateMax_Adj) # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
      minDate <- StartDte # %m+% month(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window
      
      # - Fixing to spell stop-time, we aggregate to monthly level and observe the time series up to given point
      datAggr_cohorts <- merge(datGiven[timeVar==timeVar_SpellMax, list(Sum_Total = .N), by=list(timeVar)],
                               datGiven[timeVar==timeVar_SpellMax, list(Sum_Resol = .N), by=list(timeVar,resolType)],
                               by="timeVar")[timeVar >= minDate & timeVar <= maxDate,]
      datAggr_cohorts[, EventRate := Sum_Resol/Sum_Total]
      
      # - Subset only to main resolution rate amongst subtypes (kappa)
      eventRate_prop <- subset(datAggr_cohorts, resolType == resolType_Val)
      # plot(x=eventRate_prop$timeVar, y=eventRate_prop$EventRate, type="b") # Quick visual check
    }
    
    # - Subsampled-resampled training set
    # set dates
    StartDte <- min(datGiven_train[,timeVar], na.rm=T) # First date in sampling window
    EndDte <- max(datGiven_train[,timeVar], na.rm=T) # Last date in sampling window
    maxDate <- EndDte %m+% months(sDateMax_Adj) # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
    minDate <- StartDte # %m+% month(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window
    # Fixing to spell stop-time, we aggregate to monthly level and observe the time series up to given point
    datAggr_cohorts <- merge(datGiven_train[timeVar==timeVar_SpellMax, list(Sum_Total = .N), by=list(timeVar)],
                             datGiven_train[timeVar==timeVar_SpellMax, list(Sum_Resol = .N), by=list(timeVar,resolType)],
                             by="timeVar")[timeVar >= minDate & timeVar <= maxDate,]
    datAggr_cohorts[, EventRate := Sum_Resol/Sum_Total]
    # Subset only to main resolution rate amongst subtypes (kappa)
    eventRate_prop_train <- subset(datAggr_cohorts, resolType == resolType_Val)
    # plot(x=eventRate_prop_train$timeVar, y=eventRate_prop_train$EventRate, type="b") # Quick visual check
    
    
    # - Subsampled-resampled validation set
    # set dates
    StartDte <- min(datGiven_valid[,timeVar], na.rm=T) # First date in sampling window
    EndDte <- max(datGiven_valid[,timeVar], na.rm=T) # Last date in sampling window
    maxDate <- EndDte %m+% months(sDateMax_Adj) # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
    minDate <- StartDte # %m+% month(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window
    # Fixing to spell stop-time, we aggregate to monthly level and observe the time series up to given point
    datAggr_cohorts <- merge(datGiven_valid[timeVar==timeVar_SpellMax, list(Sum_Total = .N), by=list(timeVar)],
                             datGiven_valid[timeVar==timeVar_SpellMax, list(Sum_Resol = .N), by=list(timeVar,resolType)],
                             by="timeVar")[timeVar >= minDate & timeVar <= maxDate,]
    datAggr_cohorts[, EventRate := Sum_Resol/Sum_Total]
    # Subset only to main resolution rate amongst subtypes (kappa)
    eventRate_prop_valid <- subset(datAggr_cohorts, resolType == resolType_Val)
    # plot(x=eventRate_prop_valid$timeVar, y=eventRate_prop_valid$EventRate, type="b") # Quick visual check
    
    # - Merge resolution rate datasets 
    # Merge on date so that the MAE can be calculated only across those
    # time points that coincide in both the population and sample (train/validation)
    datMrgd_Train <- merge(eventRate_prop[,list(timeVar,EventRate=EventRate)], 
                           eventRate_prop_train[,list(timeVar,EventRate_Samp=EventRate)], by="timeVar")
    datMrgd_Valid <- merge(eventRate_prop[,list(timeVar,EventRate=EventRate)], 
                           eventRate_prop_valid[,list(timeVar,EventRate_Samp=EventRate)], by="timeVar")
    datMrgd_TrainValid <- merge(eventRate_prop_train[,list(timeVar,EventRate=EventRate)], 
                           eventRate_prop_valid[,list(timeVar,EventRate_Samp=EventRate)], by="timeVar")
    
    # - Compare event rates across different sets using chosen error measure
    err_eventRate_MAE_train <- mean(abs(datMrgd_Train$EventRate - datMrgd_Train$EventRate_Samp), na.rm=T) # mean absolute error
    err_eventRate_MAE_valid <- mean(abs(datMrgd_Valid$EventRate - datMrgd_Valid$EventRate_Samp), na.rm=T) # mean absolute error
    err_eventRate_MAE_trainvalid <- mean(abs(datMrgd_TrainValid$EventRate - datMrgd_TrainValid$EventRate_Samp), na.rm=T) # mean absolute error
    
    # record elapsed time
    tme <- (proc.time() - ptm)
    
    # - Append error value to output table
    datTemp <- data.table(datTemp, "Err_EventRate_PopTrain_MAE" = err_eventRate_MAE_train, 
                          "Err_EventRate_PopValid_MAE" = err_eventRate_MAE_valid,
                          "Err_EventRate_TrainValid_MAE" = err_eventRate_MAE_trainvalid,
                          "ElapsedTime"=tme[3])
  }
  
  
  # --- Return value of chosen error measure
  return(datTemp)
} # end of function


# - Testing function call
ptm <- proc.time() #IGNORE: for computation time calculation
test <- subSmp_strat(smp_size=90000, smp_frac=smp_frac, seed=1, prior_prop=NA, eventRate_prop=NA, 
                     datGiven=datCredit, EventRate=T)
proc.time() - ptm  #IGNORE: for computation time calculation
# Stratification with smp_size=15000 and seed=1:
#     Optimal Sampling : SubSample_Size = 89903; Err_EventRate_PopTrain_MAE = 0.02073174; Err_EventRate_PopValid_MAE = 0.03206947; Err_EventRate_TrainValid_MAE = 0.03962074
#     Resampling Tool achieves very similar result

# --- Cleanup
rm(test, ptm)


# --- Main Loop (outer function call)
cl.port <- makeCluster(cpu.threads-1)
registerDoParallel(cl.port)

cat(paste0("1 (", Sys.time(),"). Iterating across subsample sizes ..."),
    file="subsampleLoop.txt", append=F)

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
    temp <- subSmp_strat(smp_size=smp_size_v[iSize], smp_frac=smp_frac, seed=seed_v[iSeed], 
                         prior_prop=prior_prop, eventRate_prop=eventRate_prop, datGiven=datCredit, EventRate=T)
    
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
statPoint_PopTrain <- 18 # Matter investigated; there is "jaggedness" in the plot of the second derivative, and that is why we take the second point
statPoint_TrainValid <- 22 # 2nd derivative is neither smooth nor monotonic
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
    annotate(geom="table", x=4000000, y=0.2, family=chosenFont, size=2.9,
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
ggsave(g1, file=paste0(genFigPath, "DefaultRates_SubSampleRates_Experiment.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


# --- Cleanup
rm(datCredit, g1); gc()




