# ===================================== DATA FUSION =====================================
# Subsample the full credit dataset before fusing it with the base macroeconomic
# variables from the macroeconomic data and then fusing it to the wider "input space" data.
# Finally, apply feature engineering to create some new variables for the input space.
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Esmerelda Oberholzer, Roelinde Bester, Marcel Muller
# ---------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2d(i).Data_Fusion.R

# -- Inputs:
#   - datInput.raw | raw input space imported in script 1
#   - datCredit_real | prepared credit data from script 2d(i)
#   - various parameters set in the setup script 0
#   - datMV | prepared feature engineered macroeconomic data from script 2a
#
# -- Outputs:
#   - datCredit_finalc | enriched credit dataset, fused with macroeconomic data and
#                        various input fields, with feature engineering applied
# ---------------------------------------------------------------------------------------
# NOTE: This script predominantly comes from another project (Kasmeer).
# =======================================================================================


# ------ 1. Preliminaries
# --- Load in Dataset
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

# --- Some feature engineering
# - Max date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Max_Date := max(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date:= NA]

# - Max date of each default spell (used as a stratifier and for identifying FALSE default spells)
datCredit_real[!is.na(DefSpell_Key), DefSpell_Max_Date := max(Date, na.rm=T), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), DefSpell_Max_Date:= NA]

# - Min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Min_Date := min(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Min_Date:= NA]

# - Min date of each default spell (used as a stratifier)
datCredit_real[!is.na(DefSpell_Key), DefSpell_Min_Date := min(Date, na.rm=T), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), DefSpell_Min_Date:= NA]

# - Creating a curing indicator for default spells (enables graphing of curing events in default spells)
datCredit_real[, Cured_Ind := ifelse(!is.na(DefSpell_Key) & DefSpellResol_Type_Hist=="Cured" & Date==DefSpell_Max_Date,1,0)] # | Reduce conditions to only the second and third (check feasibility)

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  TRUE ~ "Settled & Other"))
# Checking the proportions of the newly created variable
datCredit_real$PerfSpellResol_Type_Hist2 %>% table() %>% prop.table()

# --- Identifying single observation default spells as they are not taken into account directly through the subsampling/ resampling scheme
datCredit_real[, DefSpell_Exc := F] # Creating a variable for identifying these observations
datCredit_real[DefSpell_Counter==1 & Date==DefSpell_Max_Date, DefSpell_Exc := T] # Identifying FALSE default spells
# - Checking how many single-observation default spells exist and assessing their impact
(check.1 <- datCredit_real[DefSpell_Exc==T, .N] / datCredit_real[DefSpell_Counter==1,.N])
cat(sprintf("%.4f", check.1), "% (", datCredit_real[DefSpell_Exc==T, .N], " of ", datCredit_real[DefSpell_Counter==1,.N], ")", "of default spells are to be excluded and are indirectly taken in to account via the subsampling/ resampling scheme.")

# --- Clean up
rm(check.1); gc()




# ------ 2. Subsampling
# --- Preliminaries
# - Field names
stratifiers_Perf <- c("PerfSpell_Max_Date", "PerfSpellResol_Type_Hist2") #First variable should be of type "date"
stratifiers_Def <- c("DefSpell_Max_Date", "DefSpellResol_Type_Hist")
targetVar <- c("DefaultStatus1_lead_12_max") # Field name of the main target (i.e., the 12-month default rate)
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
perfMax_timeVar <- "PerfSpell_Max_Date" # Field name of last observed time of the associated performance spell
defMax_timeVar <- "DefSpell_Max_Date" # Field name of last observed time of the associated default spell
# - Other parameters
smp_size <- 65000 # fixed size of downsampled set in terms of the number of unique performance spells
cat(smp_size, " is ", sprintf("%.4f", smp_size/length(unique(datCredit_real[,get(clusVar_Perf)]))*100), "% of all performance spells.")
smp_frac <- 0.7 # sampling fraction for resampling scheme
minStrata_size_Perf <- 20 # Minimum size for performance spells
minStrata_size_Def <- 5 # Minimum size for default spells
# - Implied sampling fraction for downsampling step
smp_perc <- smp_size/length(unique(datCredit_real[,get(clusVar_Perf)]))

# --- Subsampling
# - Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
dat_keys_perf <- unique(datCredit_real[!is.na(clusVar_Perf), mget(c(clusVar_Perf, stratifiers_Perf))])
dat_keys_def <- unique(datCredit_real[get(excDef)==F & !is.na(get(clusVar_Def)), mget(c(clusVar_Def, stratifiers_Def))])

# - Setting the seed
set.seed(1)
# - Use simple random sampling with the stratifiers to select the performance- and default spell IDs that ought to be in the subsampled dataset
dat_smp_keys_perf <- dat_keys_perf %>% group_by(across(all_of(stratifiers_Perf))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
dat_smp_keys_def <- dat_keys_def %>% group_by(across(all_of(stratifiers_Def))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
# - Create two subsampled datasets from the sampled performance- and default spell IDs
datCredit_smp_perf <- datCredit_real %>% subset(get(clusVar_Perf) %in% dat_smp_keys_perf[, get(clusVar_Perf)])
datCredit_smp_def <- datCredit_real %>% subset(get(clusVar_Def) %in% dat_smp_keys_def[, get(clusVar_Def)])
# - Fuse the subsampled datasets
datCredit_smp <- funion(datCredit_smp_perf, datCredit_smp_def, all = F) # Using a union to concatenate the two training datasets (duplicate rows are removed)
# - Arranging the subsample and removing defExc (this will be recreated for the resampling scheme)
datCredit_smp <- datCredit_smp %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
# - Creating a variable for identifying FALSE performance- and default spells
datCredit_smp <- False_Perf_Def(datCredit_smp, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                PerfSpell_Max_Date = perfMax_timeVar, DefSpell_Max_Date = defMax_timeVar)

# - [CHECK] Proportion of default spells in the full dataset and the subsequent training- and validation datasets
check.2_a <- datCredit_real[DefSpell_Counter==1,.N]/(datCredit_real[PerfSpell_Counter==1, .N] + datCredit_real[DefSpell_Counter==1, .N])
check.2_b <- datCredit_smp[DefSpell_Counter==1, .N]/(datCredit_smp[PerfSpell_Counter==1, .N] + datCredit_smp[DefSpell_Counter==1 & !is.na(DefSpell_Key), .N])
cat(sprintf("%.4f", check.2_a*100), "% of all spells in the full dataset are default spells. \n",
    sprintf("%.4f", check.2_b*100), "% of all spells in the subsampled dataset are performance spells. \n")

# - Identifying single observation default spells as they are not taken into account directly through the resampling scheme
datCredit_smp[, (excDef) := F] # Creating a variable for identifying these observations
datCredit_smp[get(defCounter)==1 & get(timeVar)==get(defMax_timeVar), (excDef) := T] # Identifying FALSE default spells
# Checking how many single observations default spells exist and assessing their impact
(check.3 <- datCredit_smp[get(excDef)==T, .N] / datCredit_smp[get(defCounter)==1,.N])
cat(sprintf("%.4f", check.3), "% (", datCredit_smp[get(excDef)==T, .N], " of ", datCredit_smp[get(defCounter)==1,.N], ")", "of default spells are to be excluded and are indirectly taken into account via the subsampling/ resampling scheme.")

# - Clean up
rm(dat_smp_keys_perf, dat_smp_keys_def,  datCredit_smp_perf, datCredit_smp_def, check.2_a, check.2_b, check.3)


# --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold for PERFORMANCE SPELLS
# - Obtaining the stratum that are below the minimum
selectionVar_smp_Perf <- c(clusVar_Perf, timeVar, stratifiers_Perf)
datStrata_smp_min_Perf <- datCredit_smp[get(perfCounter)==1, mget(selectionVar_smp_Perf)][, list(Freq = .N), by=stratifiers_Perf][Freq<minStrata_size_Perf,]
cat(sum(datStrata_smp_min_Perf[,Freq]), "accounts of ", datCredit_smp[get(perfCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Perf[,Freq])/datCredit_smp[get(perfCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Perf)

# - Conditionally applying the exclusions
if (sum(datStrata_smp_min_Perf[,Freq]) > 0){
  # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
  datCredit_smp_old_n_Perf <- datCredit_smp[!is.na(get(clusVar_Perf)),.N]; datCredit_smp_prior_Perf <- (datCredit_smp[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table())[[2]] # NOTE: Change this indexing number manually for the target variable
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
    
    dat_keys_exc_Perf <- c(dat_keys_exc_Perf, as.vector(datCredit_smp[eval(excCond2_Perf), get(clusVar_Perf)]))
  }
  dat_keys_exc_Perf <- dat_keys_exc_Perf[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
  
  # Applying the exclusions to the subsampled dataset
  datCredit_smp <- copy(datCredit_smp[!(get(clusVar_Perf) %in% dat_keys_exc_Perf),])
  
  cat(datCredit_smp_old_n_Perf-datCredit_smp[!is.na(get(clusVar_Perf)),.N], " observations removed (", sprintf("%.4f", (datCredit_smp_old_n_Perf-datCredit_smp[!is.na(get(clusVar_Perf)),.N])/datCredit_smp_old_n_Perf*100), "% ) \n",
      "Prior probability = ", sprintf("%.4f", datCredit_smp_prior_Perf*100), "% comapred to ", sprintf("%.4f", (datCredit_smp[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table())[[2]]*100), "%")
}
# - Obtaining the stratum that are below the minimum
datStrata_smp_min_Perf <- datCredit_smp[get(perfCounter)==1, mget(selectionVar_smp_Perf)][, list(Freq = .N), by=stratifiers_Perf][Freq<minStrata_size_Perf,]
cat(sum(datStrata_smp_min_Perf[,Freq]), "accounts of ", datCredit_smp[get(perfCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Perf[,Freq])/datCredit_smp[get(perfCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Perf)

# - Cleanup
suppressWarnings(rm(selectionVar_smp_Perf, datStrata_smp_min_Perf, datCredit_smp_old_n_Perf, datCredit_smp_prior_Perf, class_type_Perf, excCond_Perf, excCond2_Perf))


# --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold for DEFAULT SPELLS
# - Obtaining the stratum that are below the minimum
selectionVar_smp_Def <- c(clusVar_Def, timeVar, stratifiers_Def)
datStrata_smp_min_Def <- datCredit_smp[get(defCounter)==1, mget(selectionVar_smp_Def)][, list(Freq = .N), by=stratifiers_Def][Freq<minStrata_size_Def,]
cat(sum(datStrata_smp_min_Def[,Freq]), "accounts of ", datCredit_smp[get(defCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Def[,Freq])/datCredit_smp[get(defCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Def)

# - Conditionally applying the exclusions
if (sum(datStrata_smp_min_Def[,Freq]) > 0){
  # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
  datCredit_smp_old_n_Def <- datCredit_smp[!is.na(get(clusVar_Def)),.N]; datCredit_smp_prior_Def <- (datCredit_smp[get(defCounter)==1, get(resolDef)] %>% table() %>% prop.table())[[2]] # NOTE: Change this indexing number manually for the target variable
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
    
    dat_keys_exc_Def <- c(dat_keys_exc_Def, as.vector(datCredit_smp[eval(excCond2_Def), get(clusVar_Def)]))
  }
  dat_keys_exc_Def <- dat_keys_exc_Def[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
  
  # Applying the exclusions to the subsampled dataset
  datCredit_smp <- copy(datCredit_smp[!(get(clusVar_Def) %in% dat_keys_exc_Def),])
  
  cat(datCredit_smp_old_n_Def-datCredit_smp[!is.na(get(clusVar_Def)),.N], " observations removed (", sprintf("%.4f", (datCredit_smp_old_n_Def-datCredit_smp[!is.na(get(clusVar_Def)),.N])/datCredit_smp_old_n_Def*100), "% ) \n",
      "Prior probability = ", sprintf("%.4f", datCredit_smp_prior_Def*100), "% comapred to ", sprintf("%.4f", (datCredit_smp[get(defCounter)==1, get(resolDef)] %>% table() %>% prop.table())[[2]]*100), "%")
}
# - Obtaining the stratum that are below the minimum
datStrata_smp_min_Def <- datCredit_smp[get(defCounter)==1, mget(selectionVar_smp_Def)][, list(Freq = .N), by=stratifiers_Def][Freq<minStrata_size_Def,]
cat(sum(datStrata_smp_min_Def[,Freq]), "accounts of ", datCredit_smp[get(defCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Def[,Freq])/datCredit_smp[!is.na(get(clusVar_Def)) & get(defCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Def)

# - Clean up
suppressWarnings(rm(selectionVar_smp_Def, datStrata_smp_min_Def, datCredit_smp_old_n_Def, datCredit_smp_prior_Def, class_type_Def, excCond_Def, excCond2_Def,
                    dat_keys_smp_perf, dat_keys_smp_defm, datCredit_real))




# ------- 3. Fusing credit with macroeconomic information
# - Loan in datasets
if (!exists('datMV')) unpack.ffdf(paste0(genPath,"datMV"), tempPath)

# - Find intersection between fields in the credit dataset and the macroeconomic dataset
(overlap_flds <- intersect(colnames(datCredit_smp), colnames(datMV))) # no overlapping fields except Date

# - Merge on Date by performing a left-join | Only subsetting the fundamental/ base macroeconomic variables; the lags therof are fused in script 4a.Data_Fusion2.
datCredit_smp <- merge(datCredit_smp, datMV, by="Date", all.x=T); gc()

# - Create Interest Rate margin using the repo rate + 3.5% (Prime Rate's definition in South Africa)
datCredit_smp <- datCredit_smp %>% mutate(InterestRate_Margin = round(InterestRate_Nom - (M_Repo_Rate+0.035), digits=4)) %>%
  relocate(InterestRate_Margin, .after=InterestRate_Nom)

# - Create AgeToTerm variable
datCredit_smp <- datCredit_smp %>% mutate(AgeToTerm = Age_Adj/Term)

# - Create BalanceToTerm variable
datCredit_smp <- datCredit_smp %>% mutate(BalanceToTerm = Balance/Term)

# - Validate merging success by checking for missingness (should be zero)
list_merge_variables <- list(colnames(datMV))
results_missingness <- list()
for (i in 1:length(list_merge_variables)){
  output <- sum(is.na(datCredit_smp$list_merge_variables[i]))
  results_missingness[[i]] <- output
}
length(which(results_missingness > 0)) == 0 # confirmed, no missing values

# - Clean-up
rm(datMV, list_merge_variables, results_missingness); gc()




# ------- 4. Fusing credit dataset with additional input fields
# - Confirm if the input space data is loaded into memory
if (!exists('datInput.raw')) unpack.ffdf(paste0(genPath,"creditdata_input1"), tempPath)

# [SANITY CHECK] Prevalence of overlapping fields in the input space and the main credit dataset
# Find intersection between fields in input space and those perhaps already in the main credit dataset
overlap_flds <- intersect(colnames(datCredit_smp), colnames(datInput.raw))
check.fuse1 <- length(overlap_flds) == 0 # FALSE; duplicate columns exists.
cat(check.fuse1 %?% 'SAFE: No overlapping fields in the input space and the main credit dataset' %:%
      'WARNING: Overlapping field(s) detected in the input space and the main credit dataset.')
# Conditional reporting
if (check.fuse1 == 0) {cat('NOTE: The following fields overlap: ', overlap_flds,"\n",sep="\t")}

# - Remove any additional variables that are not going to be used
suppressWarnings( datInput.raw[, `:=`(slc_status_final_pred7 = NULL, slc_status_final = NULL, 
                                      slc_curing_ind = NULL, datex = NULL)])

# - Ensure variables are not present in dataset before fusion (useful during debugging)
suppressWarnings( datCredit_smp[, `:=`(slc_pmnt_method = NULL, slc_past_due_amt = NULL, slc_days_excess = NULL,
                                       slc_status_final_pred7 = NULL, slc_status_final = NULL, slc_curing_ind = NULL,
                                       slc_acct_pre_lim_perc = NULL, slc_acct_roll_ever_24 = NULL,
                                       slc_acct_arr_dir_3 = NULL, slc_acct_prepaid_perc_dir_12 = NULL, 
                                       ccm_ute_lvl_40_cnt_24m = NULL, ccm_worst_arrears_6m = NULL, ccm_worst_arrears_24m = NULL)])

# - Format the date in the correct format for merging
datInput.raw[, date := as.Date(date, format="%Y-%m-%d")]

# - Rename the datasets for merging
colnames(datInput.raw)[colnames(datInput.raw) %in% c("date", "acct_no")] <- c("Date", "LoanID")

# - Check the data grain
data_grain_check <- datInput.raw[, list(Freq = .N), by=list(LoanID, Date)][Freq>1,]
sum(is.na(data_grain_check$LoanID))
# the data grain is broken in the cases where a Loan_ID does not exist - we are not interested in these accounts in any case

# - Merge on LoanID and Date by performing a left-join
datCredit_smp <- merge(datCredit_smp, datInput.raw, by=c("Date", "LoanID"), all.x=T); gc()

# - Check the data grain
NROW(data_grain_check_merge <- datCredit_smp[, list(Freq = .N), by=list(LoanID, Date)][Freq>1,])==0
# success, the data grain check is passed

# - Save intermediary snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_smp_a"), datCredit_smp)

# - Clean-up
rm(datInput.raw, data_grain_check, data_grain_check_merge); gc()




# ------- 5. Feature engineering for modelling purposes (SLC input space)
# - Confirm if the main dataset is loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_smp_a"), tempPath)

# --- Missing value diagnostics & treatments
# - Missing value indicators for the input space variables
# NOTE: There are a lot of missing values for these variables because of system changes etc.
datCredit_smp[, value_ind_slc_past_due_amt := ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 0, 1)]
datCredit_smp[, value_ind_slc_days_excess := ifelse(is.na(slc_days_excess) | slc_days_excess == "", 0, 1)]
datCredit_smp[, value_ind_slc_acct_pre_lim_perc := ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 0, 1)]
datCredit_smp[, value_ind_slc_acct_prepaid_perc_dir_12 := ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 0, 1)]

# - Check the missingness of the variables
# If they are more than 50% missing - remove
table(datCredit_smp$value_ind_slc_past_due_amt) %>% prop.table()             # missingness: 11.58% - keep the variable (numeric)
table(datCredit_smp$value_ind_slc_days_excess) %>% prop.table()              # missingness: 74.37% - discard the variable
table(datCredit_smp$value_ind_slc_acct_pre_lim_perc) %>% prop.table()        # missingness: 11.58 - keep the variable (numeric)
table(datCredit_smp$value_ind_slc_acct_prepaid_perc_dir_12) %>% prop.table() # missingness: 11.68% - keep the variable (numeric)

# - Remove the variables that have missingness > 50%
suppressWarnings( datCredit_smp[, `:=`(value_ind_slc_days_excess = NULL, slc_days_excess = NULL, 
                                        value_ind_ccm_ute_lvl_40_cnt_24m = NULL, ccm_ute_lvl_40_cnt_24m = NULL,
                                        value_ind_ccm_worst_arrears_6m = NULL, ccm_worst_arrears_6m = NULL,
                                        value_ind_ccm_worst_arrears_24m = NULL, ccm_worst_arrears_24m = NULL)]); gc()


# --- Computing basic statistics of each variable
# - Categorical variables
varList_SLC_Cat <- c('slc_pmnt_method','slc_acct_arr_dir_3')
varSLC_Info_Cat <- describe(subset(datCredit_smp, select = varList_SLC_Cat))

# - Numeric variables
varList_SLC_Num <- c('slc_past_due_amt','slc_acct_pre_lim_perc','slc_acct_prepaid_perc_dir_12','slc_acct_roll_ever_24')
varSLC_Info_Num <- describe(subset(datCredit_smp, select = varList_SLC_Num))

# - [SANITY CHECK] - Check for ovelaps between the categorical and numeric variables
check.fuse3 <- intersect(varList_SLC_Cat,varList_SLC_Num)
cat(if(identical(check.fuse3, character(0))) {'SAFE: No overlapping fields in the categorical and numerical variables'} else {
                                             'WARNING: Overlapping field(s) detected in the input space and the main credit dataset.'})

# --- Analysis and Treatments of Categorical variables
# - [slc_pmnt_method] - Payment method
varSLC_Info_Cat$slc_pmnt_method # [slc_pmnt_method] has 7 levels and 506737 observations have missing values.
# [TREATMENT] Binning "Unknown" values and missing values into one level
datCredit_smp[, slc_pmnt_method := 
                 ifelse(is.na(slc_pmnt_method) | slc_pmnt_method == "" | slc_pmnt_method == "Unknown",
                        "MISSING_DATA", slc_pmnt_method)]
(varSLC_Info_Cat$slc_pmnt_method <- describe(datCredit_smp$slc_pmnt_method)) # [slc_pmnt_method] has 7 levels and 0 observations have missing values.
# [TREATMENT] Apply factor transformation
datCredit_smp[,slc_pmnt_method:=factor(slc_pmnt_method)]
### RESULTS: ~ Missing values imputed and facorisation applied to the levels of the variable.
###            No further treatment necessary

# - [slc_acct_arr_dir_3] - Account-level arrears direction vs three months ago
varSLC_Info_Cat$slc_acct_arr_dir_3 # [slc_acct_arr_dir_3] has 4 levels and 506737 observations have missing values.
# [TREATMENT] Binning "N/A" values and missing values into one level
datCredit_smp[, slc_acct_arr_dir_3 := 
                 ifelse(is.na(slc_acct_arr_dir_3) | slc_acct_arr_dir_3 == "" | slc_acct_arr_dir_3 == "N/A", 
                        "MISSING_DATA", slc_acct_arr_dir_3)]
(varSLC_Info_Cat$slc_acct_arr_dir_3 <- describe(datCredit_smp$slc_acct_arr_dir_3)) # [slc_acct_arr_dir_3] has 4 levels and 0 observations have missing values.
# [TREATMENT] Apply factor transformation
datCredit_smp[,slc_acct_arr_dir_3:=factor(slc_acct_arr_dir_3)]
### RESULTS: ~ Missing values imputed and facorisation applied to the levels of the variable.
###            No further treatment necessary



# --- Analysis and Treatments of Numeric variables
# - [slc_past_due_amt]
varSLC_Info_Num$slc_past_due_amt # [slc_past_due_amt] has scale [0;2545589.79] and 506737 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 2718; note a VERY large single outlier in the right hand-side tail
# [TREATMENT] Median imputation for missing values (due to non-normality, i.e. high skewness) - Median imputation proceeds from the distribution of values which constitutes performace spells only (all associated values from default spells are excluded)
datCredit_smp[, slc_past_due_amt_imputed := 
                 ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 
                        median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_past_due_amt) | slc_past_due_amt == ""), slc_past_due_amt], na.rm=TRUE), slc_past_due_amt)]
(varSLC_Info_Num$slc_past_due_amt_imputed <- describe(datCredit_smp$slc_past_due_amt_imputed)) # [slc_past_due_amt_imputed] has scale [0;2545589.79] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 2403
hist(datCredit_smp$slc_past_due_amt_imputed, breaks=500); skewness(datCredit_smp$slc_past_due_amt_imputed, na.rm = T); datCredit_smp[slc_past_due_amt_imputed==0,.N]/datCredit_smp[,.N] # [slc_past_due_amt_imputed] is skewed to the right; Skewness = 27.67188; 91.32243% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values) - Winsorisation only applied to observation in default spell
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_past_due_amt) | slc_past_due_amt == ""), slc_past_due_amt], probs = 0.999)
datCredit_smp[, slc_past_due_amt_imputed_wins := ifelse(slc_past_due_amt_imputed < wins_quant, slc_past_due_amt_imputed, wins_quant)]
(varSLC_Info_Num$slc_past_due_amt_imputed_wins <- describe(datCredit_smp$slc_past_due_amt_imputed_wins)) # [slc_past_due_amt_imputed_wins] has scale [0;32027.06] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1012
hist(datCredit_smp$slc_past_due_amt_imputed_wins, breaks=500); skewness(datCredit_smp$slc_past_due_amt_imputed_wins, na.rm = T); datCredit_smp[slc_past_due_amt_imputed_wins==0,.N]/datCredit_smp[,.N] # [slc_past_due_amt_imputed_wins] is skewed to the right; Skewness = 5.521292; 91.32243% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with no shifting due to high prevalence of zero-values)
datCredit_smp[, slc_past_due_amt_imputed_wins_mm_scaled := scaler(slc_past_due_amt_imputed_wins,shift = F)]
(varSLC_Info_Num$slc_past_due_amt_imputed_wins_mm_scaled <- describe(datCredit_smp$slc_past_due_amt_imputed_wins_mm_scaled)) 
hist(datCredit_smp$slc_past_due_amt_imputed_wins_mm_scaled, breaks=500)
### RESULTS: ~ Median imputation used to treat missing values.
###            Imputed variable was scaled using min-max scaling with no shifting.
###            [slc_past_due_amt_imputed_wins_mm_scaled] has scale [0;1] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.03161
###            No further treatment necessary


# - [slc_acct_pre_lim_perc] - Prepaid/available funds to limit
varSLC_Info_Num$slc_acct_pre_lim_perc # [slc_acct_pre_lim_perc] has scale [0;1] and 506737 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.09307
# [TREATMENT] Median imputation for missing values (due to non-normality, i.e. high skewness) - Median imputation proceeds from the distribution of values which constitutes performace spells only (all associated values from default spells are excluded)
datCredit_smp[, slc_acct_pre_lim_perc_imputed := 
                 ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 
                        median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == ""), slc_acct_pre_lim_perc], na.rm=TRUE), slc_acct_pre_lim_perc)]
(varSLC_Info_Num$slc_acct_pre_lim_perc_imputed <- describe(datCredit_smp$slc_acct_pre_lim_perc_imputed))
hist(datCredit_smp$slc_acct_pre_lim_perc_imputed, breaks=500); skewness(datCredit_smp$slc_acct_pre_lim_perc_imputed, na.rm = T); datCredit_smp[slc_acct_pre_lim_perc_imputed==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 3.099905; 69.4518% of variables have zero values.
### RESULTS: ~ Median imputation used to treat missing values.
###            No outliers present, the variable has a reasonable scale and hence no scaling is applied.
###            [slc_acct_pre_lim_perc_imputed] has scale [0;1] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.07596
###            No further treatment necessary

# - [slc_acct_roll_ever_24]
varSLC_Info_Num$slc_acct_roll_ever_24 # [slc_acct_roll_ever_24] has scale [0;4] and 507101 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.4893
# [TREATMENT] Median imputation for missing values -  Median imputation proceeds from the distribution of values which constitutes performace spells only (all associated values from default spells are excluded)
datCredit_smp[, slc_acct_roll_ever_24_imputed := 
                 ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "" | slc_acct_roll_ever_24 == "Unknown", 
                        median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "") ,slc_acct_roll_ever_24], na.rm=TRUE), slc_acct_roll_ever_24)]
(varSLC_Info_Cat$slc_acct_roll_ever_24_imputed <- describe(datCredit_smp$slc_acct_roll_ever_24_imputed))
hist(datCredit_smp[,slc_acct_roll_ever_24_imputed], breaks=500); skewness(datCredit_smp$slc_acct_roll_ever_24_imputed, na.rm = T); datCredit_smp[slc_acct_roll_ever_24_imputed==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 2.501518; 82.32882% of variables have zero values.
### RESULTS: ~ Median imputation used to treat missing values.
###            No outliers present, since the variable represents counts and only has 4 levels.
###            For the same reason there, no scaling is applied.
###            [slc_acct_roll_ever_24] has scale [0;4] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.4326
###            No further treatment necessary

# - [slc_acct_prepaid_perc_dir_12] The prepaid/available funds of an account compared to 12 months ago (as a %)
varSLC_Info_Num$slc_acct_prepaid_perc_dir_12 # [slc_acct_pre_lim_perc] has scale [0;3609540717775.48] and 506737 observations have missing values; with 0 at 50% quantile and 0.07214 at 75% quantile and mean of 1216054
# [TREATMENT] Median imputation for missing values - Median imputation proceeds from the distribution of values which constitutes performace spells only (all associated values from default spells are excluded)
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed := 
                 ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 
                        median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == ""), slc_acct_prepaid_perc_dir_12], na.rm=TRUE), slc_acct_prepaid_perc_dir_12)]
(varSLC_Info_Num$slc_acct_prepaid_perc_dir_12_imputed <- describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed)) # [slc_acct_prepaid_perc_dir_12_imputed] has scale [0;3609540717775.48] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1075259
hist(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed, breaks=500); skewness(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed, na.rm = T); datCredit_smp[slc_acct_prepaid_perc_dir_12_imputed==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 2066.594; 76.1361% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values) - Winsorisation only applied to observations in default spell
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == ""), slc_acct_prepaid_perc_dir_12], probs = 0.99)
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed_wins := ifelse(slc_acct_prepaid_perc_dir_12_imputed < wins_quant, slc_acct_prepaid_perc_dir_12_imputed, wins_quant)] # Note that this is applied to the top 1%, this is since the scale is still very large when applying Winsorisation to the top 0.01%/
(varSLC_Info_Num$slc_acct_prepaid_perc_dir_12_imputed_wins <- describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_wins)) # [slc_acct_prepaid_perc_dir_12_imputed_wins] has scale [0;87.56] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1.455
hist(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_wins, breaks=500); skewness(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_wins, na.rm = T); datCredit_smp[slc_acct_prepaid_perc_dir_12_imputed_wins==0,.N]/datCredit_smp[,.N] # [slc_acct_prepaid_perc_dir_12_imputed_wins] is skewed to the right; Skewness = 8.708705; 76.1361% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with no shifting due to high prevalence of zero-values)
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled := scaler(slc_acct_prepaid_perc_dir_12_imputed_wins,shift = F)]
(varSLC_Info_Num$slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled <- describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled)) 
hist(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled, breaks=500)
### RESULTS: ~ Median imputation used to treat missing values.
###            Imputed variable was scaled using min-max scaling with no shifting.
###            [slc_acct_prepaid_perc_dir_12_imputed_wins_mm_scaled] has scale [~0;1] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.01662.
###            No further treatment necessary

# - Saving description objects on the disk for fast disk-based retrieval
pack.ffdf(paste0(genObjPath,"SLC_Var_Cat_Descript"), varSLC_Info_Cat)
pack.ffdf(paste0(genObjPath,"SLC_Var_Num_Descript"), varSLC_Info_Num)

# - Clean-up
rm(varList_SLC_Cat, varSLC_Info_Cat, varList_SLC_Num, varSLC_Info_Num); gc()




# ------- 6. Feature engineering for modelling purposes (Credit data input space)
# --- Computing basic statistics of each variable
# - Categorical variables
varList_Credit_Cat <- c('g0_Delinq','DefaultStatus1', 'DefaultStatus1_lead_12', 'DefaultStatus1_lead_12_max',
                        'PerfSpell_LeftTrunc', 'PerfSpellResol_Type_Hist','HasLeftTruncPerfSpell',
                        'DefSpell_LeftTrunc', 'DefSpellResol_Type_Hist','HasLeftTruncDefSpell',
                        'Event_Type','WOff_Ind',
                        'EarlySettle_Ind','FurtherLoan_Ind','Redraw_Ind','LN_TPE','Repaid_Ind',
                        'g0_Delinq_Shift', 'PerfSpell_g0_Delinq_Shift', 'DefSpell_g0_Delinq_Shift')
varCredit_Info_Cat <- describe(subset(datCredit_smp, select = varList_Credit_Cat))

# - Numeric variables
varList_Credit_Num <- c('Age_Adj', 'PerfSpell_Num','Term','InterestRate_Nom','InterestRate_Margin','Principal','Instalment','Receipt_Inf','Arrears',
                        'Balance','TimeInPerfSpell','PerfSpell_TimeEnd','PerfSpell_Age', 'PerfSpell_Counter',
                        'TimeInDefSpell', 'DefSpellResol_TimeEnd', 'DefSpell_Age', 'DefSpell_Counter',
                        'Event_Time','FurtherLoan_Amt','Redrawn_Amt', 'BalanceToTerm', 'AgeToTerm',
                        'g0_Delinq_Num', 'PerfSpell_g0_Delinq_Num', 'DefSpell_g0_Delinq_Num',
                        'g0_Delinq_SD', 'PerfSpell_g0_Delinq_SD', 'DefSpell_g0_Delinq_SD')

# - [SANITY CHECK] - Check for ovelaps between the categorical and numeric variables
check.fuse4 <- intersect(varList_Credit_Cat,varList_Credit_Num) # No overlapping variables.
cat(if(identical(check.fuse4, character(0))) {'SAFE: No overlapping fields in the categorical and numerical variables'} else {
  'WARNING: Overlapping field(s) detected in the input space and the main credit dataset.'})

# --- Analysis and Treatments of Categorical variables
# - [g0_Delinq]
varCredit_Info_Cat$g0_Delinq # [g0_Delinq] has 4 levels and no missing values. 89.7% of observations have no delinquency and 4% of observations have delinquency status of 3.
datCredit_smp$g0_Delinq %>% table() %>% prop.table()
barplot(table(datCredit_smp$g0_Delinq)) # Most accounts have delinquency status of zero.
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [DefaultStatus1]
varCredit_Info_Cat$DefaultStatus1 # [DefaultStatus1] has 2 levels and no missing values. 94.74% of observations in default
datCredit_smp$DefaultStatus1 %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [DefaultStatus1_lead_12]
varCredit_Info_Cat$DefaultStatus1_lead_12 # [DefaultStatus1_lead_12] has 2 levels and 52673 missing values. 93.27% of observations is in default 12-months from their reporting date
datCredit_smp$DefaultStatus1_lead_12 %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [DefaultStatus1_lead_12_max]
varCredit_Info_Cat$DefaultStatus1_lead_12_max # [DefaultStatus1_lead_12_max] has 2 levels and 52673 missing values. 91.80% of observations in default in the next the 12-months from their reporting date
datCredit_smp$DefaultStatus1_lead_12_max %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [g0_Delinq_Shift]
varCredit_Info_Cat$g0_Delinq_Shift # [g0_Delinq_Shift] has 2 levels and 52673 missing values. 91.80% of observations in default in the next the 12-months from their reporting date
datCredit_smp$g0_Delinq_Shift %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [PerfSpell_g0_Delinq_Shift]
varCredit_Info_Cat$PerfSpell_g0_Delinq_Shift # [PerfSpell_g0_Delinq_Shift] has 2 levels and 52673 missing values. 91.80% of observations in default in the next the 12-months from their reporting date
datCredit_smp$PerfSpell_g0_Delinq_Shift %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [DefSpell_g0_Delinq_Shift]
varCredit_Info_Cat$DefSpell_g0_Delinq_Shift # [DefSpell_g0_Delinq_Shift] has 2 levels and 52673 missing values. 91.80% of observations in default in the next the 12-months from their reporting date
datCredit_smp$PerfSpell_g0_Delinq_Shift %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [PerfSpell_LeftTrunc] - Indicates if the account is left-truncated.
varCredit_Info_Cat$PerfSpell_LeftTrunc # [PerfSpell_LeftTrunc] has 2 levels and 207334 missing values. 45.14% of performance spells have left-truncations.
datCredit_smp$PerfSpell_LeftTrunc %>% table() %>% prop.table()
### RESULTS: ~ No further treatment necessary

# - [HasLeftTruncPerfSpell] - Indicates if the performance spell is left-truncated.
varCredit_Info_Cat$HasLeftTruncPerfSpell # [HasLeftTruncPerfSpell] has 2 levels and no missing values. 51.76% of performance spells have left-truncations.
datCredit_smp$HasLeftTruncPerfSpell %>% table() %>% prop.table()
# [TREATMENT] Apply factor transformation
datCredit_smp[,HasLeftTruncPerfSpell:=factor(HasLeftTruncPerfSpell)]
### RESULTS:~ [PerfSpell_LeftTrunc] and [HasLeftTruncPerfSpell] are statistically the same, but they are used at different levels (see script 2c).
###         ~ No further treatment necessary

# - [PerfSpellResol_Type_Hist] - Categorical variable showing how a performance spell has ended (defaulted, settled, written-off, or anything else (paid-up))
varCredit_Info_Cat$PerfSpellResol_Type_Hist # [PerfSpellResol_Type_Hist] has 5 levels and no missing values. 12.52% of performance spells defaulted; 41.56% of performance spells have been right-censored.
datCredit_smp$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
barplot(table(datCredit_smp$PerfSpellResol_Type_Hist))
# [TREATMENT] Apply factor transformation
datCredit_smp[,PerfSpellResol_Type_Hist:=factor(PerfSpellResol_Type_Hist)]
### RESULTS: ~ No further treatment necessary

# - [DefSpell_LeftTrunc] - Indicates if the account is left-truncated.
varCredit_Info_Cat$DefSpell_LeftTrunc # [DefSpell_LeftTrunc] has 2 levels and 4146416 missing values. 3.89% of default spells have left-truncations.
datCredit_smp$DefSpell_LeftTrunc %>% table() %>% prop.table()
### RESULTS: ~ No further treatment necessary

# - [HasLeftTruncDefSpell] - Indicates if the default spell is left-truncated.
varCredit_Info_Cat$HasLeftTruncDefSpell # [HasLeftTruncPerfSpell] has 2 levels and no missing values. 4.19%% of default spells have left-truncations.
datCredit_smp$HasLeftTruncDefSpell %>% table() %>% prop.table()
# [TREATMENT] Apply factor transformation
datCredit_smp[,HasLeftTruncDefSpell:=factor(HasLeftTruncDefSpell)]
### RESULTS:~ [DefSpell_LeftTrunc] and [HasLeftTruncDefSpell] seem statistically the same, but they are used at different levels.
###         ~ No further treatment necessary

# - [DefSpellResol_Type_Hist]
varCredit_Info_Cat$DefSpellResol_Type_Hist # [DefSpellResol_Type_Hist] has 3 levels and 4146416 missing values. 61.43% of default spells have cured; 59.53% of default spells are written-off
datCredit_smp$DefSpellResol_Type_Hist %>% table() %>% prop.table()
barplot(table(datCredit_smp$DefSpellResol_Type_Hist))
# [TREATMENT] Apply factor transformation
datCredit_smp[,DefSpellResol_Type_Hist:=factor(DefSpellResol_Type_Hist)]
### RESULTS: ~ No further treatment necessary

# - [Event_Type]
varCredit_Info_Cat$Event_Type # [Event_Type] has 4 levels and no missing values. 45% of loan accounts were active (censored); 47.5% were settled; 4.2% were written-off.
datCredit_smp$Event_Type %>% table() %>% prop.table()
barplot(table(datCredit_smp$Event_Type))
# [TREATMENT] Apply factor transformation
datCredit_smp[,Event_Type:=factor(Event_Type)]
### RESULTS: ~ No further treatment necessary

# - [WOff_Ind]
varCredit_Info_Cat$WOff_Ind # [WOff_Ind] has 2 levels and no missing values. 99.94% of observations did not experience write-off.
datCredit_smp$WOff_Ind %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [EarlySettle_Ind]
varCredit_Info_Cat$EarlySettle_Ind # [EarlySettle_Ind] has 2 levels and no missing values. 99.26% of observations did not experience early settlement.
datCredit_smp$EarlySettle_Ind %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [FurtherLoan_Ind]
varCredit_Info_Cat$FurtherLoan_Ind # [FurtherLoan_Ind] has 2 levels and no missing values. 99.88% of observations did not "re-mortgage".
datCredit_smp$FurtherLoan_Ind %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [Redraw_Ind]
varCredit_Info_Cat$Redraw_Ind # [Redraw_Ind] has 3 levels and no missing values. 93.15% of observations are in level 1, 3.95% are in level 2, and 2.90% of observation are in level 3.
datCredit_smp$Redraw_Ind %>% table() %>% prop.table()
barplot(table(datCredit_smp$Redraw_Ind))
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [Repaid_Ind]
varCredit_Info_Cat$Repaid_Ind # [Repaid_Ind] has 2 levels and no missing values. 99.97456% of observations experienced no repayment.
datCredit_smp$Repaid_Ind %>% table() %>% prop.table()
### RESULTS: ~ Do not apply transformation due to the special/ important nature of this field.
###            No further treatment necessary

# - [LN_TPE]
varCredit_Info_Cat$LN_TPE # [LN_TPE] has 2 levels and no missing values. 89.1% of observations are in level 'CHL' and 10.9% are in level 'WHL'
datCredit_smp$LN_TPE %>% table() %>% prop.table()
# [TREATMENT] Apply factor transformation
datCredit_smp[,LN_TPE:=factor(LN_TPE)]
### RESULTS: ~ No further treatment necessary


# - Saving description objects on the disk for fast disk-based retrieval
pack.ffdf(paste0(genObjPath,"Credit_Var_Cat_Descript"), varCredit_Info_Cat)



# --- Analysis and Treatments of Numeric variables (i) Delinquency level feature engineering
# - Indicator for when a shift in the state of g0_Delinq occurs (target event) --- Performance spell level
datCredit_real[, PerfSpell_g0_Delinq_Shift := ifelse(lag(g0_Delinq, n=1)==g0_Delinq,0,1), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_g0_Delinq_Shift), PerfSpell_g0_Delinq_Shift := 0] # All first observations have g0_Delinq_Shift = NA; set these values to zero.
datCredit_real[is.na(PerfSpell_Key), PerfSpell_g0_Delinq_Shift := NA] # Ensure that this performance spell level variable is undefined for default spells.
datCredit_real$PerfSpell_g0_Delinq_Shift %>% table() %>% prop.table() # [PerfSpell_g0_Delinq_Shift] has 2 levels and no missing values. 96.64% of observations have a shift in delinquency and 3.36% of observations have no shifts in delinquency.

# - Indicator for when a shift in the state of g0_Delinq occurs (target event) --- Default spell level
datCredit_real[, DefSpell_g0_Delinq_Shift := ifelse(lag(g0_Delinq, n=1)==g0_Delinq,0,1), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_g0_Delinq_Shift), DefSpell_g0_Delinq_Shift := 0] # All first observations have g0_Delinq_Shift = NA; set these values to zero.
datCredit_real[is.na(DefSpell_Key), DefSpell_g0_Delinq_Shift := NA] # Ensure that this default spell level variable is undefined for performance spells.
datCredit_real$DefSpell_g0_Delinq_Shift %>% table() %>% prop.table() # [DefSpell_g0_Delinq_Shift] has 2 levels and no missing values. 85.82% of observations have a shift in delinquency and  14.18% of observations have no shifts in delinquency.

# - State number --- Performance spell level
datCredit_real[, PerfSpell_g0_Delinq_Num := cumsum(PerfSpell_g0_Delinq_Shift), by=list(PerfSpell_Key)] # Assign state numbers over each performance spell
datCredit_real[, PerfSpell_g0_Delinq_Num := PerfSpell_g0_Delinq_Num + 1] # Small adjustment to ensure that spell numbers can't be zero.
datCredit_real[is.na(PerfSpell_g0_Delinq_Shift), PerfSpell_g0_Delinq_Num := NA] # Ensure that this performance spell level variable is undefined for default spells.
hist(datCredit_real$PerfSpell_g0_Delinq_Num, breaks='FD')

# - State number --- Default spell level
datCredit_real[, DefSpell_g0_Delinq_Num := cumsum(DefSpell_g0_Delinq_Shift), by=list(DefSpell_Key)] # Assign state numbers over each default spell
datCredit_real[, DefSpell_g0_Delinq_Num := DefSpell_g0_Delinq_Num + 1] # Small adjustment to ensure that spell numbers can't be zero.
datCredit_real[is.na(DefSpell_g0_Delinq_Shift), DefSpell_g0_Delinq_Num := NA] # Ensure that this default spell level variable is undefined for performance spells.
hist(datCredit_real$DefSpell_g0_Delinq_Num, breaks='FD')

# - State volatility --- Performance spell level
datCredit_real[, PerfSpell_g0_Delinq_SD := sd(g0_Delinq), by=PerfSpell_Key]
datCredit_real[is.na(PerfSpell_g0_Delinq_Shift), PerfSpell_g0_Delinq_SD := NA] # Ensure that this performance spell level variable is undefined for default spells.
hist(datCredit_real$PerfSpell_g0_Delinq_SD, breaks='FD')

# - State volatility --- Default spell level
datCredit_real[, DefSpell_g0_Delinq_SD := sd(g0_Delinq), by=DefSpell_Key]
datCredit_real[is.na(DefSpell_g0_Delinq_Shift), DefSpell_g0_Delinq_SD := NA] # Ensure that this default spell level variable is undefined for performance spells.
hist(datCredit_real$DefSpell_g0_Delinq_SD, breaks='FD')


# --- Analysis and Treatments of Numeric variables (ii) Analysis
# - Compute the descriptive statistics for the numeric variables of the credit input space
varCredit_Info_Num <- describe(subset(datCredit_smp, select = varList_Credit_Num))

# - [Age_Adj]
varCredit_Info_Num$Age_Adj # [Age_Adj] has scale [1;1305], with 68 at 50% quantile and 127 at 75% quantile and mean of 85.23.
hist(datCredit_smp$Age_Adj, breaks = 'FD') # Distribution of [Age_Adj] is skewed to the right.
### RESULTS: ~ No further treatment necessary => This variable relates to time and as such should not be transformed

# - [PerfSpell_Num]
varCredit_Info_Num$PerfSpell_Num # [PerfSpell_Num] has scale [1;8], with a mean of 1.067
hist(datCredit_smp$PerfSpell_Num, breaks = 'FD') # Distribution of [PerfSpell_Num] is skewed to the right.
### RESULTS: ~ No further treatment necessary => This variable relates to "a measure of time" and as such should not be transformed

# - [g0_Delinq_Num]
varCredit_Info_Num$g0_Delinq_Num # [g0_Delinq_Num] has scale [1;8], with a mean of 1.067
hist(datCredit_smp$g0_Delinq_Num, breaks = 'FD') # Distribution of [PerfSpell_Num] is skewed to the right.
### RESULTS: ~ No further treatment necessary => This variable relates to "a measure of time" and as such should not be transformed

# - [PerfSpell_g0_Delinq_Num]
varCredit_Info_Num$PerfSpell_g0_Delinq_Num # [PerfSpell_g0_Delinq_Num] has scale [1;8], with a mean of 1.067
hist(datCredit_smp$PerfSpell_g0_Delinq_Num, breaks = 'FD') # Distribution of [PerfSpell_Num] is skewed to the right.
### RESULTS: ~ No further treatment necessary => This variable relates to "a measure of time" and as such should not be transformed

# - [DefSpell_g0_Delinq_Num]
varCredit_Info_Num$DefSpell_g0_Delinq_Num # [DefSpell_g0_Delinq_Num] has scale [1;8], with a mean of 1.067
hist(datCredit_smp$DefSpell_g0_Delinq_Num, breaks = 'FD') # Distribution of [PerfSpell_Num] is skewed to the right.
### RESULTS: ~ No further treatment necessary => This variable relates to "a measure of time" and as such should not be transformed

# - [Term]
varCredit_Info_Num$Term # [Term] has scale [1;970], with 240 at 50% quantile and 240 at 75% quantile and mean of 238.2.
hist(datCredit_smp$Term, breaks = 'FD') # Vast majority of loans centered around 240, with some other loans with a lower term and very few with a larger term.
### RESULTS:~ No further treatment necessary => This variable relates to time and as such should not be transformed

# - [InterestRate_Nom]
varCredit_Info_Num$InterestRate_Nom # [InterestRate_Nom] has scale [0;0.2365], with 0.09 at 50% quantile and 0.1050 at 75% quantile and mean of 0.09307.
hist(datCredit_smp$InterestRate_Nom, breaks = 'FD') # Note some outliers at 0%, which may be a special case of debt restructuring to help defaulted borroweres
### RESULTS:~ No further treatment necessary

# - [InterestRate_Margin]
varCredit_Info_Num$InterestRate_Margin # [InterestRate_Margin] has scale [-0.1550;0.1080], with -0.0080 at 50% quantile and 0 at 75% quantile and mean of -0.007391.
hist(datCredit_smp$InterestRate_Margin, breaks = 'FD') # Distribution is centered around zero.
### RESULTS: ~ No further treatment necessary

# - [Principal]
varCredit_Info_Num$Principal # [Principal] has scale [0.01;208477000], with 510000 at 50% quantile and 850000 at 75% quantile and mean of 646260.
hist(datCredit_smp$Principal, breaks='FD'); skewness(datCredit_smp$Principal, na.rm = T); datCredit_smp[Principal==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 15.13721; 0% of observations have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Principal) | Principal == ""), Principal], probs = 0.999) # 99.9% quantile = 4 680 480
datCredit_smp[, Principal_wins := ifelse(Principal < wins_quant, Principal, wins_quant)] # Note that this is applied to the top 1%, this is since the scale is still very large when applying Winsorisation to the top 0.01%/
(varCredit_Info_Num$Principal_wins <- describe(datCredit_smp$Principal_wins)) # [Principal_wins] has scale [0.01;4680480.00] and 0 observations have missing values; with 510000 at 50% quantile and 850000 at 75% quantile and mean of 644358.
hist(datCredit_smp$Principal_wins, breaks='FD'); skewness(datCredit_smp$Principal_wins, na.rm = T); datCredit_smp[Principal_wins==0,.N]/datCredit_smp[,.N] # [Principal_wins] is skewed to the right; Skewness = 2.026844; 0% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with shifting)
datCredit_smp[, Principal_wins_mm_scaled := scaler(Principal_wins,shift = T)]
(varCredit_Info_Num$Principal_wins_mm_scaled <- describe(datCredit_smp$Principal_wins_mm_scaled)) 
hist(datCredit_smp$Principal_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile.
###            Winsorised variable scaled using min-max scaling with no shifting.
###            [Principal_wins_mms_scaled] has scale [0;1] and 0 observations have missing values; with 0.10896 at 50% quantile and 0.18161 at 75% quantile and mean of 0.1377.
###            No further treatment necessary

# - [Instalment]
varCredit_Info_Num$Instalment # [Instalment] has scale [0;16417873.80], with 4719.1 at 50% quantile and 7921.4 at 75% quantile and mean of 5947.
hist(datCredit_smp$Instalment, breaks ='FD'); skewness(datCredit_smp$Instalment, na.rm = T); datCredit_smp[Instalment==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 490.0671; 0.31% of observations have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Instalment) | Instalment == ""), Instalment], probs = 0.999) # 99.9% quantile = 43 296.14
datCredit_smp[, Instalment_wins := ifelse(Instalment < wins_quant, Instalment, wins_quant)]
(varCredit_Info_Num$Instalment_wins <- describe(datCredit_smp$Instalment_wins)) # [Instalment_wins] has scale [0;43296.14] and 0 observations have missing values; with 4719.1 at 50% quantile and 7921.4 at 75% quantile and mean of 5870.
hist(datCredit_smp$Instalment_wins, breaks='FD'); skewness(datCredit_smp$Instalment_wins, na.rm = T); datCredit_smp[Instalment_wins==0,.N]/datCredit_smp[,.N] # [Instalment_wins] is skewed to the right; Skewness = 2.005984; 0.03121732% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with no shifting)
datCredit_smp[, Instalment_wins_mm_scaled := scaler(Instalment_wins,shift = F)]
(varCredit_Info_Num$Instalment_wins_mm_scaled <- describe(datCredit_smp$Instalment_wins_mm_scaled)) 
hist(datCredit_smp$Instalment_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile.
###            Winsorised variable scaled using min-max scaling with no shifting.
###            [Instalment_mm_scaled] has scale [0;1] and 0 observations have missing values; with 0.10900 at 50% quantile and 0.18296 at 75% quantile and mean of 0.1356.
###            No further treatment necessary

# - [Receipt_Inf] Not an important variable; keep in script (caution in using variable)
varCredit_Info_Num$Receipt_Inf # [Receipt_Inf] has scale [0;60785555.01], with 4371 at 50% quantile and 8284 at 75% quantile and mean of 14262.
hist(datCredit_smp$Receipt_Inf, breaks='FD'); skewness(datCredit_smp$Receipt_Inf, na.rm = T); datCredit_smp[Receipt_Inf==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 64.9649; 13.22% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Receipt_Inf) | Receipt_Inf == ""), Receipt_Inf], probs = 0.99) # 99% quantile = 153 344.9
datCredit_smp[, Receipt_Inf_wins := ifelse(Receipt_Inf < wins_quant, Receipt_Inf, wins_quant)]
(varCredit_Info_Num$Receipt_Inf_wins <- describe(datCredit_smp$Receipt_Inf_wins)) # [Receipt_Inf_wins] has scale [0;153344.92] and 0 observations have missing values; with 4371 at 50% quantile and 8284 at 75% quantile and mean of 8204.
hist(datCredit_smp$Receipt_Inf_wins, breaks='FD'); skewness(datCredit_smp$Receipt_Inf_wins, na.rm = T); datCredit_smp[Receipt_Inf_wins==0,.N]/datCredit_smp[,.N] # [Receipt_Inf_wins] is skewed to the right; Skewness = 6.341937; 13.22% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with no shifting)
datCredit_smp[, Receipt_Inf_wins_mm_scaled := scaler(Receipt_Inf_wins,shift = T)]
(varCredit_Info_Num$Receipt_Inf_wins_mm_scaled <- describe(datCredit_smp$Receipt_Inf_wins_mm_scaled)) 
hist(datCredit_smp$Receipt_Inf_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 1% quantile.
###            Winsorised variable was scaled using min-max scaling with no shifting.
###            [Receipt_Inf_wins_mm_scaled] has scale [0;1] and 0 observations have missing values; with 0.028503 at 50% quantile and 0.054022 at 75% quantile and mean of 0.0535.
###            No further treatment necessary

# - [Arrears]
varCredit_Info_Num$Arrears # [Arrears] has scale [-1.96;2545589.79], with 0 at 50% quantile and 0 at 75% quantile and mean of 2495.
hist(datCredit_smp$Arrears, breaks='FD'); skewness(datCredit_smp$Arrears, na.rm = T); datCredit_smp[Arrears==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 49.40433; 90.028965% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Arrears) | Arrears == ""), Arrears], probs = 0.999) # 99.9% quantile = 31 765.79
datCredit_smp[, Arrears_wins := ifelse(Arrears < wins_quant, Arrears, wins_quant)]
(varCredit_Info_Num$Arrears_wins <- describe(datCredit_smp$Arrears_wins)) # [Arrears_wins] has scale [-1.96;31765.79] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1071.
hist(datCredit_smp$Arrears_wins, breaks='FD'); skewness(datCredit_smp$Arrears_wins, na.rm = T); datCredit_smp[Arrears_wins==0,.N]/datCredit_smp[,.N] # [Arrears_wins] is skewed to the right; Skewness = 5.340883; 90.30% of variables have zero values.
# [TREATMENT] Scale the data using standardisation (the variable has negative values and therefore min-max is inappropriate; also no shifting is used)
datCredit_smp[, Arrears_wins_norm_scaled := scaler.norm(Arrears_wins,shift = F)]
(varCredit_Info_Num$Arrears_wins_norm_scaled <- describe(datCredit_smp$Arrears_wins_norm_scaled)) 
hist(datCredit_smp$Arrears_wins_norm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.01% quantile.
###            Variable was scaled using min-max scaling with no shifting.
###            [Arrears_wins_norm_scaled] has scale [--0.000401275260;6.503483176980] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.2192.
###            No further treatment necessary

# - [Balance]
varCredit_Info_Num$Balance # [Balance] has scale [-146.71;41687282.40], with 365231.5 at 50% quantile and 692208.4 at 75% quantile and mean of 489701.
hist(datCredit_smp$Balance, breaks = 'FD'); skewness(datCredit_smp$Balance, na.rm = T); datCredit_smp[Balance==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 4.310693; 3.79% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Balance) | Balance == ""), Balance], probs = 0.999) # 99.9% quantile = 33 844 818
datCredit_smp[, Balance_wins := ifelse(Balance < wins_quant, Balance, wins_quant)]
(varCredit_Info_Num$Balance_wins <- describe(datCredit_smp$Balance_wins)) # [Balance_wins] has scale [-146.71;3844818.17] and 0 observations have missing values; with 365231.5 at 50% quantile and 692208.4 at 75% quantile and mean of 488207.
hist(datCredit_smp$Balance_wins, breaks='FD'); skewness(datCredit_smp$Balance_wins, na.rm = T); datCredit_smp[Balance_wins==0,.N]/datCredit_smp[,.N] # [Balance_wins] is skewed to the right; Skewness = 1.941194; 3.79% of variables have zero values.
# [TREATMENT] Scale the data using standardisation (the variable has negative values and therefore min-max is inappropriate; also no shifting is used)
datCredit_smp[,Balance_wins_norm_scaled := scaler.norm(Balance_wins,shift = F)]
(varCredit_Info_Num$Balance_wins_norm_scaled <- describe(datCredit_smp$Balance_wins_norm_scaled)) 
hist(datCredit_smp$Balance_wins_norm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.01% quantile.
###            Variable was scaled using min-max scaling with no shifting.
###            [Balance_wins_norm_scaled] has scale [-0.00029373479839;7.69788622964179] and 0 observations have missing values; with 0.73124668 at 50% quantile and 1.38590199 at 75% quantile and mean of 0.9775.
###            No further treatment necessary

# - [TimeInPerfSpell]
varCredit_Info_Num$TimeInPerfSpell # [TimeInPerfSpell] has scale [1;1305], with 62 at 50% quantile and 118 at 75% quantile and mean of 79.55.
hist(datCredit_smp$TimeInPerfSpell, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [PerfSpell_TimeEnd]
varCredit_Info_Num$PerfSpell_TimeEnd # [PerfSpell_TimeEnd] has scale [1;1305], with 136 at 50% quantile and 194 at 75% quantile and mean of 139.2.
hist(datCredit_smp$PerfSpell_TimeEnd, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [PerfSpell_Age]
varCredit_Info_Num$PerfSpell_Age # [PerfSpell_Age] has scale [1;1305], with 129 at 50% quantile and 191 at 75% quantile and mean of 134.5.
hist(datCredit_smp$PerfSpell_Age, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [Event_Time]   
varCredit_Info_Num$Event_Time # [Event_Time] has scale [1;1305, with 142 at 50% quantile and 197 at 75% quantile and mean of 143.8.
hist(datCredit_smp$Event_Time, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [TimeInDefSpell]
varCredit_Info_Num$TimeInDefSpell # [TimeInDefSpell] has scale [1;571], with 11 at 50% quantile and 24 at 75% quantile and mean of 21.04.
hist(datCredit_smp$TimeInDefSpell, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [DefSpellResol_TimeEnd]
varCredit_Info_Num$DefSpellResol_TimeEnd # [DefSpellResol_TimeEnd] has scale [1;1202], with 98 at 50% quantile and 178 at 75% quantile and mean of 121.4.
hist(datCredit_smp$DefSpellResol_TimeEnd, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [DefSpell_Age]
varCredit_Info_Num$DefSpell_Age # [DefSpell_Age] has scale [1;402], with 26 at 50% quantile and 47 at 75% quantile and mean of 37.64.
hist(datCredit_smp$DefSpell_Age, breaks = 'FD')
### RESULTS: ~ DO NOT APPLY ANY TRANSFORMATIONS DUE TO THE SPECIAL NATURE OF THIS VARIABLE!

# - [FurtherLoan_Amt]
varCredit_Info_Num$FurtherLoan_Amt # [FurtherLoan_Amt] has scale [0;2081864.654560], with 0 at 50% quantile and 0 at 75% quantile and mean of 196.2.
hist(datCredit_smp$FurtherLoan_Amt, breaks='FD'); skewness(datCredit_smp$FurtherLoan_Amt, na.rm = T); datCredit_smp[FurtherLoan_Amt==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 61.26655; 96.05% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with no shifting)
datCredit_smp[, FurtherLoan_Amt_mm_scaled := scaler(FurtherLoan_Amt,shift = F)]
(varCredit_Info_Num$FurtherLoan_Amt_mm_scaled <- describe(datCredit_smp$FurtherLoan_Amt_mm_scaled)) 
hist(datCredit_smp$FurtherLoan_Amt_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected, but this is expected since loans are rarely "extended".
###            Winsorised variable was scaled using min-max scaling with no shifting.
###            [FurtherLoan_Amt_mm_scaled] has scale [0;1] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.00009425.
###            No further treatment necessary

# - [Redrawn_Amt]        
varCredit_Info_Num$Redrawn_Amt # [Redrawn_Amt] has scale [0;9417751.24377500], with 0 at 50% quantile and 0 at 75% quantile and mean of 2180.
hist(datCredit_smp$Redrawn_Amt, breaks='FD'); skewness(datCredit_smp$Redrawn_Amt, na.rm = T); datCredit_smp[Redrawn_Amt==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 61.26655; 96.05% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.01% quantile to remove the influence of extreme outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Redrawn_Amt) | Redrawn_Amt == ""), Redrawn_Amt], probs = 0.9999) # 99.99% quantile = 1 075 488
datCredit_smp[, Redrawn_Amt_wins := ifelse(Redrawn_Amt < wins_quant, Redrawn_Amt, wins_quant)]
(varCredit_Info_Num$Redrawn_Amt_wins <- describe(datCredit_smp$Redrawn_Amt_wins)) # [Redrawn_Amt_wins] has scale [0;1075488.46220894] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 2120.
hist(datCredit_smp$Redrawn_Amt_wins, breaks='FD'); skewness(datCredit_smp$Redrawn_Amt_wins, na.rm = T); datCredit_smp[Redrawn_Amt_wins==0,.N]/datCredit_smp[,.N] # [Redrawn_Amt_wins] is skewed to the right; Skewness = 23.3661; 96.05% of values have zero values.
# [TREATMENT] Scale the data using min-max scaling (with shifting)
datCredit_smp[, Redrawn_Amt_wins_mm_scaled := scaler(Redrawn_Amt_wins,shift = F)]
(varCredit_Info_Num$Redrawn_Amt_wins_mm_scaled <- describe(datCredit_smp$Redrawn_Amt_wins_mm_scaled)) 
hist(datCredit_smp$Redrawn_Amt_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.01% quantile.
###            Winsorised variable was scaled using min-max scaling with no shifting.
###            [Redrawn_Amt_wins_mm_scaled] has scale [0;1] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.001971.
###            No further treatment necessary

# - [AgeToTerm]
varCredit_Info_Num$AgeToTerm # [AgeToTerm] has scale [0.002673797;22.400000000], with 0.28750 at 50% quantile and 0.53750 at 75% quantile and mean of 0.3647.
hist(datCredit_smp$AgeToTerm, breaks='FD'); skewness(datCredit_smp$AgeToTerm, na.rm = T); datCredit_smp[AgeToTerm==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 8.137409; 0% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 1% quantile to remove the influence of extreme outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(AgeToTerm) | AgeToTerm == ""), AgeToTerm], probs = 0.999) # 99.9% quantile = 2.483333
datCredit_smp[, AgeToTerm_wins := ifelse(AgeToTerm < wins_quant, AgeToTerm, wins_quant)]
(varCredit_Info_Num$AgeToTerm_wins <- describe(datCredit_smp$AgeToTerm_wins)) # [AgeToTerm_wins] has scale [0.002673797;2.483333333] and 0 observations have missing values; with 0.28750 at 50% quantile and 0.53750 at 75% quantile and mean of 0.3628.
hist(datCredit_smp$AgeToTerm_wins, breaks='FD'); skewness(datCredit_smp$AgeToTerm_wins, na.rm = T); datCredit_smp[AgeToTerm_wins==0,.N]/datCredit_smp[,.N] # [AgeToTerm_wins] is skewed to the right; Skewness = 1.489611; 0% of values have zero values.
# [TREATMENT] Scale the data using min-max scaling (with shifting)
datCredit_smp[, AgeToTerm_wins_mm_scaled := scaler(AgeToTerm_wins,shift = F)]
(varCredit_Info_Num$AgeToTerm_wins_mm_scaled <- describe(datCredit_smp$AgeToTerm_wins_mm_scaled)) 
hist(datCredit_smp$AgeToTerm_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile.
###            Winsorised variable was scaled using min-max scaling with no shifting.
###            [AgeToTerm_wins_mm_scaled] has scale [0;1] and 0 observations have missing values; with 0.11590 at 50% quantile and 0.21668 at 75% quantile and mean of 0.1462.
###            No further treatment necessary

# - [BalanceToTerm]
varCredit_Info_Num$BalanceToTerm # [BalanceToTerm] has scale [-0.61129166667;173697.00999999998], with 1540.4410 at 50% quantile and 2910.6708 at 75% quantile and mean of 2070.
hist(datCredit_smp$BalanceToTerm, breaks='FD'); skewness(datCredit_smp$BalanceToTerm, na.rm = T); datCredit_smp[BalanceToTerm==0,.N]/datCredit_smp[,.N]# Distribution is skewed to the right; Skewness = 5.934091; 3.552881% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of extreme outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(BalanceToTerm) | BalanceToTerm == ""), BalanceToTerm], probs = 0.999) # 99.9% quantile = 17 659.6 
datCredit_smp[, BalanceToTerm_wins := ifelse(BalanceToTerm < wins_quant, BalanceToTerm, wins_quant)]
(varCredit_Info_Num$BalanceToTerm_wins <- describe(datCredit_smp$BalanceToTerm_wins)) # [BalanceToTerm_wins] has scale [-0.61129166667;17659.6] and 0 observations have missing values; with 1540.4410 at 50% quantile and 2910.6708 at 75% quantile and mean of 2062.
hist(datCredit_smp$BalanceToTerm_wins, breaks='FD'); skewness(datCredit_smp$BalanceToTerm_wins, na.rm = T); datCredit_smp[BalanceToTerm_wins==0,.N]/datCredit_smp[,.N] # [BalanceToTerm_wins] is skewed to the right; Skewness = 2.087682; 3.48% of values have zero values.
# [TREATMENT] Scale the data using standardisation (the variable has negative values and therefore min-max is inappropriate; also no shifting is used)
datCredit_smp[,BalanceToTerm_wins_norm_scaled := scaler.norm(BalanceToTerm_wins,shift = F)]
(varCredit_Info_Num$BalanceToTerm_wins_norm_scaled <- describe(datCredit_smp$BalanceToTerm_wins_norm_scaled)) 
hist(datCredit_smp$BalanceToTerm_wins_norm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile.
###            Winsorised variable was scaled using min-max scaling with no shifting.
###            [BalanceToTerm_wins_mm_scaled] has scale [-0.00028717652441;8.29623799867224] and 0 observations have missing values; with 0.7236782 at 50% quantile and 1.3673936 at 75% quantile and mean of 0.9685.
###            No further treatment necessary


# --- Save object to disk
pack.ffdf(paste0(genObjPath,"Credit_Var_Num_Descript"), varCredit_Info_Num)


# --- Combine objects into a singe object and store to disk -> This object then describes all covariates
if (!exists('varSLC_Info_Cat')) unpack.ffdf(paste0(genObjPath,"SLC_Var_Cat_Descript"), tempPath)
if (!exists('varSLC_Info_Num')) unpack.ffdf(paste0(genObjPath,"SLC_Var_Num_Descript"), tempPath)
if (!exists('varCredit_Info_Cat')) unpack.ffdf(paste0(genObjPath,"Credit_Var_Cat_Descript"), tempPath)
if (!exists('varCredit_Info_Num')) unpack.ffdf(paste0(genObjPath,"Credit_Var_Num_Descript"), tempPath)

Covariate_Info <- c(varSLC_Info_Cat,varSLC_Info_Num,varCredit_Info_Cat,varCredit_Info_Num)

# --- Save object to disk
pack.ffdf(paste0(genObjPath,"All_Var_Descript"), Covariate_Info)

# --- Save data checkpoint
pack.ffdf(paste0(genPath,"creditdata_final_smp_b"), datCredit_smp)

# --- Clean up
suppressWarnings(rm(varSLC_Info_Cat, varSLC_Info_Num, varCredit_Info_Cat, varCredit_Info_Num, check.fuse1, check.fuse3, check.fuse4, lookup_IDs,
                    Covariate_Info, lookup, lookup2))
gc()




# ------ 7. Implementing a resampling scheme with spell-level clustering
# --- Implement resampling scheme using given main sampling fraction
# - Setting the seed
set.seed(1)
# - Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
dat_keys_smp_perf <- unique(datCredit_smp[!is.na(get(clusVar_Perf)), mget(c(clusVar_Perf, stratifiers_Perf))])
dat_keys_smp_def <- unique(datCredit_smp[!is.na(get(clusVar_Def)) & get(excDef)==F, mget(c(clusVar_Def, stratifiers_Def))])
# - Creating the training- and validation datasets
# Use simple random sampling from the unique performance- and default spell keys
dat_train_keys_perf <- dat_keys_smp_perf %>% group_by(across(all_of(stratifiers_Perf))) %>% slice_sample(prop=smp_frac) %>% as.data.table()
dat_train_keys_def <- dat_keys_smp_def %>% group_by(across(all_of(stratifiers_Def))) %>% slice_sample(prop=smp_frac) %>% as.data.table()
# Create two training datasets from the sampled performance- and default spell IDs
datCredit_train_perf <- datCredit_smp %>% subset(get(clusVar_Perf) %in% dat_train_keys_perf[, get(clusVar_Perf)])
datCredit_train_def <- datCredit_smp %>% subset(get(clusVar_Def) %in% dat_train_keys_def[, get(clusVar_Def)])
# Create two validation datasets from the sampled performance- and default spell IDs
datCredit_valid_perf <- datCredit_smp %>% subset(!is.na(get(clusVar_Perf)) & !(get(clusVar_Perf) %in% datCredit_train_perf[, get(clusVar_Perf)]))
datCredit_valid_def <- datCredit_smp %>% subset(!is.na(get(clusVar_Def)) & !(get(clusVar_Def) %in% datCredit_train_def[, get(clusVar_Def)]))
# Fuse the two training- and validation datasets
datCredit_train <- funion(datCredit_train_perf, datCredit_train_def, all = F) # Using a union to concatenate the two training datasets (duplicate rows are removed)
datCredit_valid <- funion(datCredit_valid_perf, datCredit_valid_def, all = F) # Using a union to concatenate the two validation datasets (duplicate rows are removed)
# - Arranging the subsample and removing defExc (this will be recreated for the resampling scheme)
datCredit_train <- datCredit_train %>% subset(select=-which(names(datCredit_train) %in% c(excDef, "PerfSpell_F", "DefSpell_F"))) %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
datCredit_valid <- datCredit_valid %>% subset(select=-which(names(datCredit_valid) %in% c(excDef, "PerfSpell_F", "DefSpell_F"))) %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
# - Creating a variable for identifying FALSE performance- and default spells
datCredit_train <- False_Perf_Def(datCredit_train, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                  Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                  PerfSpell_Max_Date = perfMax_timeVar, DefSpell_Max_Date = defMax_timeVar)
datCredit_valid <- False_Perf_Def(datCredit_valid, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                  Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                  PerfSpell_Max_Date = perfMax_timeVar, DefSpell_Max_Date = defMax_timeVar)

# - [CHECK] Proportion of default spells in the full dataset and the subsequent training- and validation datasets
check.4_a <- datCredit_smp[get(defCounter)==1,.N]/(datCredit_smp[get(perfCounter)==1, .N] + datCredit_smp[get(defCounter)==1, .N])
check.4_b <- datCredit_train[get(defCounter)==1, .N]/(datCredit_train[get(perfCounter)==1, .N] + datCredit_train[get(defCounter)==1, .N])
check.4_c <- datCredit_valid[get(defCounter)==1, .N]/(datCredit_valid[get(perfCounter)==1, .N] + datCredit_valid[get(defCounter)==1, .N])
cat(sprintf("%.4f", check.4_a*100), "% of all spells in the full dataset are default spells. \n",
    sprintf("%.4f", check.4_b*100), "% of all spells in the training dataset are default spells. \n",
    sprintf("%.4f", check.4_c*100), "% of all spells in the validation dataset are default spells. \n")
# - [CHECK]
check.5_a <- datCredit_smp[!is.na(get(clusVar_Perf)), .N] == datCredit_train[!is.na(get(clusVar_Perf)), .N] + datCredit_valid[!is.na(get(clusVar_Perf)), .N]
check.5_b <- datCredit_smp[!is.na(get(clusVar_Def)), .N] == datCredit_train[!is.na(get(clusVar_Def)), .N] + datCredit_valid[!is.na(get(clusVar_Def)), .N]
cat(ifelse((check.5_a & check.5_b), "SAFE: Resampling scheme reconstitues subsample dataset", "WARNING: Not all observations accounted for in resampling scheme"))


# --- Saving the cross-validation scheme
# - Training dataset
pack.ffdf(paste0(genPath,"creditdata_train"), datCredit_train)
# - Validation dataset
pack.ffdf(paste0(genPath,"creditdata_valid"), datCredit_valid)


# --- Clean up
suppressWarnings(rm(dat_keys_smp_perf, dat_keys_smp_perf,  dat_train_keys_perf, dat_train_keys_def, datCredit_train_perf, datCredit_train_def,  datCredit_valid_perf, datCredit_valid_def,
                    check.4_a, check.4_b, check.4_c, check.5_a, check.5_b))





















































# ------ 2. Default spell exclusions
# --- Load in Dataset
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final4d"), tempPath)

# --- Applying the exclusions
datCredit_smp <- subset(datCredit_smp, !is.na(PerfSpell_Key))

# --- [SANITY CHECK]
# - Creating checks
sum_defaults <- sum(is.na(datCredit_smp$PerfSpell_Num) > 0)
check.fuse1 <- sum_defaults == 0

# - Reporting
cat(check.fuse1 %?% 'SAFE: All default spells removed, fusion successfull.' %:%
      'WARNING: Default spells detected, fusion not successfull.')

# Conditional reporting
if (check.fuse1 == 0) {
  cat('NOTE: ', check.fuse1, 'observations from default spells detected',"\n",sep="\t")
}

# --- Saving endpoint
pack.ffdf(paste0(genPath,"creditdata_analytics"), datCredit_smp)

