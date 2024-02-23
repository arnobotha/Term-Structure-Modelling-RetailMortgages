# ======================== OPTIMAL SUBSAMPLE SIZES FOR SURVIVAL MODELS USING ACS ========================
# Determining the effect of a wide range of subsample sizes within a simple advanced sampling resampling
# scheme amidst a survival modelling setup (PD-modelling).
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Marcel Muller

# DESCRIPTION:
# This ancillary & exploratory script iterates across a given vector of subsample sizes towards
# resampling data into a basic cross-validation setup (training:validation), using 2-way stratified 
# sampling. Each chosen size is executed multiple times over various seed values to account for randomness
# in following a broader Monte Carlo setup. Thereafter, various error measures are calculated within 
# each iteration, whereupon these error values are appropriately aggregated and graphed into a 
# single cohesive graph.
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup
#   - 2d(i).Data_Fusion

# -- Inputs:
#   - datCredit_real | Prepared from script 2d.
#
# -- Outputs:
#   - Difference graph across subsample sizes, where 'difference' is the value of an error measure 
#       between prior probabilities across resampled sets (full:training and training:validation)
# ------------------------------------------------------------------------------------------------------





# ------ 1. Preliminaries
# --- Load in Dataset
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

# --- Some feature engineering
# - Max date of each performance spell (used as a stratifier)
# - Max- and min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), c("PerfSpell_Min_Date","PerfSpell_Max_Date") := as.list(range(Date, na.rm=TRUE)), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), c("PerfSpell_Min_Date","PerfSpell_Max_Date") := as.list(c(NA,NA)), by=list(PerfSpell_Key)]
# Sanity check - Should be TRUE
datCredit_real[!is.na(PerfSpell_Key),.N] == datCredit_real[!is.na(PerfSpell_Max_Date),.N] # TRUE, field created successfully
datCredit_real[!is.na(PerfSpell_Key),.N] == datCredit_real[!is.na(PerfSpell_Min_Date),.N] # TRUE, field created successfully
# - Max- and min date of each default spell (used as a stratifier and for identifying FALSE default spells)
datCredit_real[!is.na(DefSpell_Key), c("DefSpell_Min_Date","DefSpell_Max_Date") := as.list(range(Date, na.rm=TRUE)), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), c("DefSpell_Min_Date","DefSpell_Max_Date") := as.list(c(NA,NA)), by=list(DefSpell_Key)]
# Sanity check - Should be TRUE
datCredit_real[!is.na(DefSpell_Key),.N] == datCredit_real[!is.na(DefSpell_Max_Date),.N] # TRUE, field created successfully
datCredit_real[!is.na(DefSpell_Key),.N] == datCredit_real[!is.na(DefSpell_Min_Date),.N] # TRUE, field created successfully

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "Settled & Other",
                                                                                  TRUE ~ NA))
# Checking the proportions of the newly created variable
datCredit_real$PerfSpellResol_Type_Hist2 %>% table() %>% prop.table()

# --- Identifying single observation default spells as they are not taken into account directly through the subsampling/ resampling scheme
datCredit_real[, DefSpell_Exc := F] # Creating a variable for identifying these observations
datCredit_real[DefSpell_Counter==1 & Date==DefSpell_Max_Date, DefSpell_Exc := T] # Identifying FALSE default spells
# - Checking how many single-observation default spells exist and assessing their impact
(check.1 <- datCredit_real[DefSpell_Exc==T, .N] / datCredit_real[DefSpell_Counter==1,.N])
cat(sprintf("%.4f", check.1), "% (", datCredit_real[DefSpell_Exc==T, .N], " of ", datCredit_real[DefSpell_Counter==1,.N], ")", "of default spells are to be excluded and are indirectly taken in to account via the subsampling/ resampling scheme.")


# --- Field names (REMOVE STRATIFICATION TO ALLIGN TO THE SIMPLE CLUSTERED SAMPLING TECHNIQUE || )
# Required field names
targetVar <- "DefaultStatus1_lead_12_max" # Field name of the main target (i.e., the 12-month default rate)
CurStatus <- "DefaultStatus1"
resolPerf <- "PerfSpellResol_Type_Hist2" # Field name of performance spell resolution types - first level should be the target event (default)
resolDef <- "DefSpellResol_Type_Hist" # Field name of default spell resolution types
excDef <- "DefSpell_Exc" # Field for identifying single-observation default spells (they are indirectly excluded from the resampling scheme)
clusVar <- "LoanID" # Field name of unique account identifier
clusVar_Perf <- "PerfSpell_Key" # Field name of unique performance spell identifier
clusVar_Def <- "DefSpell_Key" # Field name of unique default spell identifier
timeVar <- "Date" # Field name of time variable
counter <- "Counter" # Field name of counter for loan observations
perfCounter <- "PerfSpell_Counter" # Field name of counter for performance spell observations
defCounter <- "DefSpell_Counter" # Field name of counter for default spell observations
# Names of the stratifiers (optional)
stratifiers_Perf <- c("PerfSpell_Max_Date", "PerfSpellResol_Type_Hist2") # c("PerfSpell_Min_Date", "PerfSpellResol_Type_Hist2") #First variable should be of type "date"
stratifiers_Def <- c("DefSpell_Max_Date", "DefSpellResol_Type_Hist") # c("DefSpell_Min_Date", "DefSpellResol_Type_Hist")
# Final Selection
selectionVar <- unique(c(clusVar, clusVar_Perf, clusVar_Def, timeVar, CurStatus, counter, perfCounter, defCounter, resolPerf, resolDef,
                         excDef, stratifiers_Perf, stratifiers_Def, targetVar)) # Variables to subset
selectionVar <- selectionVar[!is.na(selectionVar)] # Facilitating cases where the variables are left unspecified (specifically for use of no stratifiers)

# --- Subset given dataset accordingly; an efficiency enhancement
datCredit <- subset(datCredit_real, select=selectionVar)
rm(datCredit_real); gc()

# --- Spell-level feature engineering
# - Max- and min date of each performance spell (used as a stratifier)
datCredit[!is.na(get(clusVar_Perf)), c("timeVar_Perf_Min","timeVar_Perf_Max") := as.list(range(get(timeVar), na.rm=TRUE)), by=list(get(clusVar_Perf))]
datCredit[is.na(get(clusVar_Perf)), c("timeVar_Perf_Min","timeVar_Perf_Max") := as.list(c(NA,NA)), by=list(get(clusVar_Perf))]
# - Max- and min date of each default spell (used as a stratifier and for identifying FALSE default spells)
datCredit[!is.na(get(clusVar_Def)), c("timeVar_Def_Min","timeVar_Def_Max") := as.list(range(get(timeVar), na.rm=TRUE)), by=list(get(clusVar_Def))]
datCredit[is.na(get(clusVar_Def)), c("timeVar_Def_Min","timeVar_Def_Max") := as.list(c(NA,NA)), by=list(get(clusVar_Def))]


# --- Calculate prior probability of the (conditional) target event. 
# NOTE: Precalculating the prior probability of the target variable - This is merely a coding optimisation
resolPerf_targetVar <- "Defaulted" # Reference level in the performance spell resolution corresponding to the target variable"
prior_prop <- datCredit[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
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
smp_size_v <- c(# seq(from=2500, to=25000, by=2500), seq(from=30000, to=75000, by=5000))
                seq(from=82500, to=105000, by=7500)) #, seq(from=165000, to=255000, by=15000)),
                # seq(from=285000, to=465000, by=30000), seq(from=505000, to=655000, by=50000))
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
# Counter_Perf: name of field counting the sequential observaiton of an performance spell
# timeVar: field name of date for event rate calculations;
# seed: specific seed value
# prior_prop: pre-calculated prior probability within population for error measurement
# eventRate_prop: pre-calculated performance spell event rates over time within population for error measurement
# minStrata_size: Minimum strata size which should be enforced in the resulting subsample
# datGiven: given dataset from which to subsample and resample
# EventRate: indicator for computing the event rate as the second error measure
subSmp_strat_adv <- function(smp_size, smp_frac, stratifiers_Perf=NA, stratifiers_Def=NA, targetVar=NA, clusVar=NA, clusVar_Perf, clusVar_Def, CurStatus=NA, Counter=NA,
                             perfCounter=NA, defCounter=NA, timeVar=NA, seed=123, prior_prop=NA, eventRate_prop=NA, minStrata_size_Perf=0, minStrata_size_Def=0, datGiven, EventRate=F) {
  # --- UNIT TEST
  # datGiven <- copy(datCredit); smp_size <- 15000; Counter_Perf <- "PerfSpell_Counter"; seed <- 1; smp_frac=0.7; minStrata_size_Perf <- 0; minStrata_size_Def <- 0; EventRate <- T
  
  # --- Preliminaries
  # - Error Checks
  if (any(is.na(stratifiers_Perf) | is.na(stratifiers_Def)) & is.na(targetVar)) { stop("Stratifiers and target variable are unspecified! Must at least include the target variable")}
  # - Implied sampling fraction for downsampling step
  smp_perc <- smp_size/length(unique(datGiven[,get(clusVar_Perf)]))
  
  
  # --- Downsample data into a set with a fixed size (using stratified sampling) before implementing resampling scheme
  # - Set seed
  set.seed(seed)
  # - Conditional loop for stratifiers
  if ((all(is.na(stratifiers_Perf)) | all(is.na(stratifiers_Def)))){ # - No stratifiers
    # Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
    dat_keys_perf <- unique(datGiven[!is.na(clusVar_Perf), get(clusVar_Perf)]) %>% as.data.table() %>% tidyr::drop_na(); colnames(dat_keys_perf) <- clusVar_Perf
    dat_keys_def <- unique(datGiven[get(excDef)==F & !is.na(get(clusVar_Def)), get(clusVar_Def)]) %>% as.data.table() %>% tidyr::drop_na(); colnames(dat_keys_def) <- clusVar_Def
    # Use simple random sampling with the stratifiers to select the performance- and default spell IDs that ought to be in the subsampled dataset
    dat_smp_keys_perf <- dat_keys_perf %>% slice_sample(prop=smp_perc) %>% as.data.table()
    dat_smp_keys_def <- dat_keys_def %>% slice_sample(prop=smp_perc) %>% as.data.table()
  } else { # - Stratifiers
    # Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
    dat_keys_perf <- unique(datGiven[!is.na(clusVar_Perf), mget(c(clusVar_Perf, stratifiers_Perf))]) %>% as.data.table() %>% tidyr::drop_na()
    dat_keys_def <- unique(datGiven[get(excDef)==F & !is.na(get(clusVar_Def)), mget(c(clusVar_Def, stratifiers_Def))]) %>% as.data.table() %>% tidyr::drop_na()
    # Use simple random sampling with the stratifiers to select the performance- and default spell IDs that ought to be in the subsampled dataset
    dat_smp_keys_perf <- dat_keys_perf %>% group_by(across(all_of(stratifiers_Perf))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
    dat_smp_keys_def <- dat_keys_def %>% group_by(across(all_of(stratifiers_Def))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
  }
  
  # - Create two subsampled datasets from the sampled performance- and default spell IDs
  datGiven_smp_perf <- datGiven %>% subset(get(clusVar_Perf) %in% dat_smp_keys_perf[, get(clusVar_Perf)])
  datGiven_smp_def <- datGiven %>% subset(get(clusVar_Def) %in% dat_smp_keys_def[, get(clusVar_Def)])
  # - Fuse the subsampled datasets
  datGiven_smp <- funion(datGiven_smp_perf, datGiven_smp_def, all = F) # Using a union to concatenate the two training datasets (duplicate rows are removed)
  # - Arranging the subsample and removing defExc (this will be recreated for the resampling scheme)
  datGiven_smp <- datGiven_smp %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
  # - Creating a variable for identifying FALSE performance- and default spells
  datGiven_smp <- False_Perf_Def(datGiven_smp, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                 Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                 PerfSpell_Max_Date = "timeVar_Perf_Max", DefSpell_Max_Date = "timeVar_Def_Max")
  # - Identifying single observation default spells as they are not taken into account directly through the resampling scheme
  datGiven_smp[, (excDef) := F] # Creating a variable for identifying these observations
  datGiven_smp[get(defCounter)==1 & get(timeVar)==timeVar_Def_Max, (excDef) := T] # Identifying FALSE default spells
  
  
  # --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold
  if (!(all(is.na(stratifiers_Perf)) | all(is.na(stratifiers_Def)))){
    # - Obtaining the stratum that are below the minimum
    selectionVar_smp_Perf <- c(clusVar_Perf, timeVar, stratifiers_Perf)
    datStrata_smp_min_Perf <- datGiven_smp[get(perfCounter)==1, mget(selectionVar_smp_Perf)][, list(Freq = .N), by=stratifiers_Perf][Freq<minStrata_size_Perf,]
    # - Conditionally applying the exclusions
    if (sum(datStrata_smp_min_Perf[,Freq]) > 0){
      # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
      datGiven_smp_old_n_Perf <- datGiven_smp[!is.na(get(clusVar_Perf)),.N]; datGiven_smp_prior_Perf <- (datGiven_smp[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table())[[2]] # NOTE: Change this indexing number manually for the target variable
      # Initiating a vector which will contain the exclusion IDs
      dat_keys_exc_Perf <- NA
      # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
      for (i in 1:datStrata_smp_min_Perf[,.N]){
        class_type_Perf <- sapply(datStrata_smp_min_Perf[,1:length(stratifiers_Perf)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
        
        excCond_Perf <- datStrata_smp_min_Perf[i,1:length(stratifiers_Perf)] # Getting the values of the ith minimum strata
        excCond_Perf <- data.table(Stratifier = colnames(excCond_Perf), # Building a dataset
                                   Value = unname(t(excCond_Perf)),
                                   Class = class_type_Perf)
        excCond_Perf[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
        excCond_Perf[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
        excCond2_Perf <- parse(text = paste0(paste0(excCond_Perf$Condition, collapse = ""), perfCounter, "==1")) # Compiling the ith condition
        
        dat_keys_exc_Perf <- c(dat_keys_exc_Perf, as.vector(datGiven_smp[eval(excCond2_Perf), get(clusVar_Perf)]))
      }
      dat_keys_exc_Perf <- dat_keys_exc_Perf[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
      
      # Applying the exclusions to the subsampled dataset
      datGiven_smp <- copy(datGiven_smp[!(get(clusVar_Perf) %in% dat_keys_exc_Perf),])
      
    }
  
    # - Cleanup
    suppressWarnings(rm(selectionVar_smp_Perf, datGiven_smp_old_n_Perf, datGiven_smp_prior_Perf, class_type_Perf, excCond_Perf, excCond2_Perf))
    
    
    # --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold for DEFAULT SPELLS
    # - Obtaining the stratum that are below the minimum
    selectionVar_smp_Def <- c(clusVar_Def, timeVar, stratifiers_Def)
    datStrata_smp_min_Def <- datGiven_smp[get(defCounter)==1, mget(selectionVar_smp_Def)][, list(Freq = .N), by=stratifiers_Def][Freq<minStrata_size_Def,]
    
    # - Conditionally applying the exclusions
    if (sum(datStrata_smp_min_Def[,Freq]) > 0){
      # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
      datGiven_smp_old_n_Def <- datGiven_smp[!is.na(get(clusVar_Def)),.N]; datGiven_smp_prior_Def <- (datGiven_smp[get(defCounter)==1, get(resolDef)] %>% table() %>% prop.table())[[2]] # NOTE: Change this indexing number manually for the target variable
      # Initiating a vector which will contain the exclusion IDs
      dat_keys_exc_Def <- NA
      # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
      for (i in 1:datStrata_smp_min_Def[,.N]){
        class_type_Def <- sapply(datStrata_smp_min_Def[,1:length(stratifiers_Def)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
        
        excCond_Def <- datStrata_smp_min_Def[i,1:length(stratifiers_Def)] # Getting the values of the ith minimum strata
        excCond_Def <- data.table(Stratifier = colnames(excCond_Def), # Building a dataset
                                  Value = unname(t(excCond_Def)),
                                  Class = class_type_Def)
        excCond_Def[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
        excCond_Def[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
        excCond2_Def <- parse(text = paste0(paste0(excCond_Def$Condition, collapse = ""), defCounter, "==1")) # Compiling the ith condition
        
        dat_keys_exc_Def <- c(dat_keys_exc_Def, as.vector(datGiven_smp[eval(excCond2_Def), get(clusVar_Def)]))
      }
      dat_keys_exc_Def <- dat_keys_exc_Def[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
      
      # Applying the exclusions to the subsampled dataset
      datGiven_smp <- copy(datGiven_smp[!(get(clusVar_Def) %in% dat_keys_exc_Def),])
    }
    # - Clean up
    suppressWarnings(rm(selectionVar_smp_Def, datGiven_smp_old_n_Def, datGiven_smp_prior_Def, class_type_Def, excCond_Def, excCond2_Def))
  }
  
  # --- Implementing the resampling scheme
  # - Setting the seed
  set.seed(1)
  # - Conditional loop for stratifiers
  if ((all(is.na(stratifiers_Perf)) | all(is.na(stratifiers_Def)))){ # No stratifiers
    # - Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
    dat_keys_smp_perf <- unique(datGiven_smp[!is.na(get(clusVar_Perf)), get(clusVar_Perf)]) %>% as.data.table() %>% tidyr::drop_na(); colnames(dat_keys_smp_perf) <- clusVar_Perf
    dat_keys_smp_def <- unique(datGiven_smp[!is.na(get(clusVar_Def)) & get(excDef)==F, get(clusVar_Def)]) %>% as.data.table() %>% tidyr::drop_na(); colnames(dat_keys_smp_def) <- clusVar_Def
    # - Creating the training- and validation datasets
    # Use simple random sampling from the unique performance- and default spell keys
    dat_train_keys_perf <- dat_keys_smp_perf %>% slice_sample(prop=smp_frac) %>% as.data.table()
    dat_train_keys_def <- dat_keys_smp_def %>% slice_sample(prop=smp_frac) %>% as.data.table()
  } else { # Stratifiers
    # - Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
    dat_keys_smp_perf <- unique(datGiven_smp[!is.na(get(clusVar_Perf)), mget(c(clusVar_Perf, stratifiers_Perf))])
    dat_keys_smp_def <- unique(datGiven_smp[!is.na(get(clusVar_Def)) & get(excDef)==F, mget(c(clusVar_Def, stratifiers_Def))])
    # - Creating the training- and validation datasets
    # Use simple random sampling from the unique performance- and default spell keys
    dat_train_keys_perf <- dat_keys_smp_perf %>% group_by(across(all_of(stratifiers_Perf))) %>% slice_sample(prop=smp_frac) %>% as.data.table()
    dat_train_keys_def <- dat_keys_smp_def %>% group_by(across(all_of(stratifiers_Def))) %>% slice_sample(prop=smp_frac) %>% as.data.table()
  } # Create two training datasets from the sampled performance- and default spell IDs
  datGiven_train_perf <- datGiven_smp %>% subset(get(clusVar_Perf) %in% dat_train_keys_perf[, get(clusVar_Perf)])
  datGiven_train_def <- datGiven_smp %>% subset(get(clusVar_Def) %in% dat_train_keys_def[, get(clusVar_Def)])
  # Create two validation datasets from the sampled performance- and default spell IDs
  datGiven_valid_perf <- datGiven_smp %>% subset(!is.na(get(clusVar_Perf)) & !(get(clusVar_Perf) %in% datGiven_train_perf[, get(clusVar_Perf)]))
  datGiven_valid_def <- datGiven_smp %>% subset(!is.na(get(clusVar_Def)) & !(get(clusVar_Def) %in% datGiven_train_def[, get(clusVar_Def)]))
  # Fuse the two training- and validation datasets
  datGiven_train <- funion(datGiven_train_perf, datGiven_train_def, all = F) # Using a union to concatenate the two training datasets (duplicate rows are removed)
  datGiven_valid <- funion(datGiven_valid_perf, datGiven_valid_def, all = F) # Using a union to concatenate the two validation datasets (duplicate rows are removed)
  # - Arranging the subsample and removing defExc (this will be recreated for the resampling scheme)
  datGiven_train <- datGiven_train %>% subset(select=-which(names(datGiven_train) %in% c(excDef, "PerfSpell_F", "DefSpell_F"))) %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
  datGiven_valid <- datGiven_valid %>% subset(select=-which(names(datGiven_valid) %in% c(excDef, "PerfSpell_F", "DefSpell_F"))) %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
  # - Creating a variable for identifying FALSE performance- and default spells
  datGiven_train <- False_Perf_Def(datGiven_train, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                    Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                    PerfSpell_Max_Date = "timeVar_Perf_Max", DefSpell_Max_Date = "timeVar_Def_Max")
  datGiven_valid <- False_Perf_Def(datGiven_valid, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                    Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                    PerfSpell_Max_Date = "timeVar_Perf_Max", DefSpell_Max_Date = "timeVar_Def_Max")
  # - [CHECK] Reconstitution of the cross-validation scheme
  # check.5_a <- datGiven_smp[!is.na(get(clusVar_Perf)), .N] == datGiven_train[!is.na(get(clusVar_Perf)), .N] + datGiven_valid[!is.na(get(clusVar_Perf)), .N] # Performance spells
  # check.5_b <- datGiven_smp[!is.na(get(clusVar_Def)), .N] == datGiven_train[!is.na(get(clusVar_Def)), .N] + datGiven_valid[!is.na(get(clusVar_Def)), .N] # Default spells
  # cat(ifelse((check.5_a & check.5_b), "SAFE: Resampling scheme reconstitues subsample dataset", "WARNING: Not all observations accounted for in resampling scheme"))
  
  # --- Create the output table
  if (!(all(is.na(stratifiers_Perf)) | all(is.na(stratifiers_Def)))){ # Stratifiers
    datTemp <- data.table("SubSample_Size"=datGiven_smp[get(counter)==1,.N], "SubSample_FullSize"=datGiven_smp[,.N],
                          "Exclusions_Perf"=sum(datStrata_smp_min_Perf$Freq), "Exclusions_Def"=sum(datStrata_smp_min_Def$Freq))
  } else { # No stratifiers
    datTemp <- data.table("SubSample_Size"=datGiven_smp[get(counter)==1,.N], "SubSample_FullSize"=datGiven_smp[,.N])
  }
  
  # --- Calculate error measure 1: Difference in prior probabilities between population and training (as subsampled + resampled)
  # Calculate prior probabilities within each relevant dataset, e.g., proportion of defaults across all time
  # - Population
  if (is.na(prior_prop)) {
    prior_prop <- datGiven[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
    prior_prop <- prior_prop[names(prior_prop)[names(prior_prop) == resolPerf_targetVar]][[1]] # Subsetting for only the target variable
  }
  
  # - Subsampled & resampled training set
  prior_prop_train <- datGiven_train[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table() # Computing the prior probabilities of the performance spell resolution outcomes
  prior_prop_train <- prior_prop_train[names(prior_prop_train)[names(prior_prop_train) == resolPerf_targetVar]][[1]] # Subsetting for only the target variable
  
  # - Compare population with training set using chosen error measure
  err_priorProb_AE <- abs(prior_prop - prior_prop_train) # absolute error
  err_priorProb_SqrdErr <- (prior_prop - prior_prop_train)^2 # squared error
  
  # - Append output table (preliminary)
  datTemp <- data.table(datTemp, "SampleFrac"=smp_frac, "Stratifiers"=paste(stratifiers_Perf, collapse="; "),
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
test <- subSmp_strat_adv(smp_size=65000, smp_frac=smp_frac, stratifiers_Perf=stratifiers_Perf, stratifiers_Def=stratifiers_Def, targetVar=targetVar, clusVar=clusVar,
                         clusVar_Perf=clusVar_Perf, clusVar_Def=clusVar_Def, CurStatus=CurStatus, Counter=counter, perfCounter=perfCounter, defCounter=defCounter, timeVar=timeVar,
                         seed=1, prior_prop=NA, eventRate_prop=NA, minStrata_size_Perf=0, minStrata_size_Def=0, datGiven=datCredit, EventRate=T)
proc.time() - ptm  #IGNORE: for computation time calculation
### RESULTS:  Resampling Tool : SubSample_Size = 437900; Exclusions = 0; Err_EventRate_PopTrain_MAE = 0.000526; Err_EventRate_PopValid_MAE = 0.000718; Err_EventRate_TrainValid_MAE = 0.001026
### RESULTS:  Optimal Sampling: SubSample_Size = 4378819 ; Exclusions_Perf = 0; Exclusions_Def = 0; Err_EventRate_PopTrain_MAE = 0.0005263027; Err_EventRate_PopValid_MAE = 0.0007184041; Err_EventRate_TrainValid_MAE = 0.001026239
### CONCLUCION: SAFE, the functions in this script reproduces the results of the advanced resampling tool


# --- Cleanup
rm(test, ptm); gc()


# --- Main Loop (outer function call)
cl.port <- makeCluster(cpu.threads-1)
registerDoParallel(cl.port)

cat(paste0("1 (", Sys.time(),"). Iterating across subsample sizes ..."),
    file="subsampleLoop_Adv.txt", append=F)

ptm <- proc.time() #IGNORE: for computation time calculation

# - Multithreaded looping procedure using the foreach-package
datResults <- foreach(it=1:(length(seed_v)*length(smp_size_v)), .combine='rbind', .verbose=F, .inorder=T, 
                      .packages=c("dplyr","data.table","tidyselect","lubridate", "scales")) %dopar%
  {
    # - Testing 
    # it <- 1
    
    # - Set indices
    iSeed <- (it-1) %% length(seed_v) + 1 # modulo operation
    iSize <- (it-1) %/% length(seed_v) + 1 # integer-valued division
    
    # - Iterate 
    temp <- subSmp_strat_adv(smp_size=smp_size_v[iSize], smp_frac=smp_frac, stratifiers_Perf=stratifiers_Perf, stratifiers_Def=stratifiers_Def, targetVar=targetVar, clusVar=clusVar,
                             clusVar_Perf=clusVar_Perf, clusVar_Def=clusVar_Def, CurStatus=CurStatus, Counter=counter, perfCounter=perfCounter, defCounter=defCounter, timeVar=timeVar,
                             seed=seed_v[iSeed], prior_prop=NA, eventRate_prop=NA, minStrata_size_Perf=0, minStrata_size_Def=0, datGiven=datCredit, EventRate=T)
    
    # - Reporting
    if (iSeed == length(seed_v)) {
      cat(paste0("\n2 (", Sys.time(),"). Subsample size: ", comma(smp_size_v[iSize]), " tested ",length(seed_v), " times."),
          file="subsampleLoop_Adv.txt", append=T) 
    }
    
    return(temp)
  }  

t <- proc.time() - ptm  #IGNORE: for computation time calculation
cat(paste0("\n3 (", Sys.time(),"). ForEach-loop done. Elapsed time: ", round(t[3]/60), " minutes."),
    file="subsampleLoop_Adv.txt", append=T)

# - Save to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0("C:/Users/R5532132/OneDrive - FRG/", "subSampleSizes_Adv2"), datResults); gc()
stopCluster(cl.port)




# ------ 3. Graphing
# --- Load in Dataset
if (!exists('datResults')) unpack.ffdf(paste0(genObjPath,"subSampleSizes_Adv"), tempPath)
# - Confidence interval parameter
# confLevel <- 0.95

# datResults1 <- copy(datResults); rm(datResults); gc()
# if (!exists('datResults')) unpack.ffdf(paste0(genObjPath,"subSampleSizes_Adv2"), tempPath)
# datResults2 <- copy(datResults); rm(datResults); gc()
# datResults <- rbind(datResults1, datResults2); rm(datResults1, datResults2); gc()

# --- Populating the results dataset for each iteration
datResults[, SubSample_Size := as.vector(sapply(X=smp_size_v, function(x) rep(x, length(seed_v))))]

# --- Aggregate to subsample size level
datGraph <- datResults[, list(PriorProb_MAE = mean(Err_PriorProb_AE , na.rm=T), PriorProb_MAE_SD = sd(Err_PriorProb_AE , na.rm=T),
                              PriorProb_RMSE = sqrt(sum(Err_PriorProb_SqrdErr, na.rm=T)/.N), PriorProb_RMSE_SE = sd(Err_PriorProb_SqrdErr, na.rm=T),
                              EventRate_PopTrain_MAE_Mean = mean(Err_EventRate_PopTrain_MAE, na.rm=T), EventRate_PopTrain_MAE_SD = sd(Err_EventRate_PopTrain_MAE, na.rm=T),
                              EventRate_TrainValid_MAE_Mean = mean(Err_EventRate_TrainValid_MAE, na.rm=T), EventRate_TrainValid_MAE_SD = sd(Err_EventRate_TrainValid_MAE, na.rm=T),
                              SubSample_Size_Mean = round(mean(SubSample_Size, na.rm=T)), SubSample_FullSize_Mean = round(mean(SubSample_FullSize, na.rm=T)), N=.N),
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
datAnnotate <- datGraph[, list(`"Size "*italic(s)`=comma(SubSample_FullSize_Mean), 
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
datEventRate_PopTrain_1st <- datGraph[order(SubSample_Size), list(Gradient = diff(EventRate_PopTrain_MAE_Mean) / diff(SubSample_Size_Mean))]$Gradient
datEventRate_TrainValid_1st <- datGraph[order(SubSample_Size), list(Gradient = diff(EventRate_TrainValid_MAE_Mean) / diff(SubSample_Size_Mean))]$Gradient
# plot(diff(datEventRate_PopTrain_1st), type="b", main="2nd derivative") # 2nd derivative neither smooth nor monotonic
# plot(diff(datEventRate_TrainValid_1st), type="b", main="2nd derivative") # 2nd derivative is smooth but not monotonic
# - Find index of stationary points x such that f''(x) <= epislon
statPoint_PopTrain <- which(diff(datEventRate_PopTrain_1st) < 10^(-10.25))[2]
statPoint_TrainValid <- which(diff(datEventRate_TrainValid_1st) < 10^(-10.5))[1] # 2nd derivative is neither smooth nor monotonic
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
(g1 <- ggplot(datGraph3, aes(x=SubSample_FullSize_Mean, y=Value, group=Set)) + theme_minimal() + 
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
               nudge_y=0.0002, nudge_x=1050000, size=2) + 
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
