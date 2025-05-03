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




# ------- 5. Feature engineering for modelling purposes
# --- 5.1 Confirm if the main dataset is loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_smp_a"), tempPath)


# --- 5.2 Computing basic statistics of each variable
# - Categorical variables
varList_Cat <- c("DefaultStatus1", "DefaultStatus1_lead_12", "DefaultStatus1_lead_12_max",
                 "g0_Delinq_Shift", "PerfSpellResol_Type_Hist",
                 "Event_Type","WOff_Ind", "EarlySettle_Ind", "FurtherLoan_Ind", "Redraw_Ind", "LN_TPE", "Repaid_Ind",
                 "slc_pmnt_method", "slc_acct_arr_dir_3")
varInfo_Cat <- describe(subset(datCredit_smp, select = varList_Cat))

# - Numeric variables
varList_Num <- c("g0_Delinq", "g0_Delinq_Num", "g0_Delinq_SD", "g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12",
                 "Age_Adj", "PerfSpell_Num", "Term", "InterestRate_Nom", "InterestRate_Margin", "Principal", "Instalment", "Arrears",
                 "Balance", "TimeInPerfSpell", "PerfSpell_TimeEnd", "PerfSpell_Age", "PerfSpell_Counter",
                 "Redrawn_Amt", "BalanceToTerm", "AgeToTerm",
                 "slc_past_due_amt", "slc_acct_pre_lim_perc", "slc_acct_prepaid_perc_dir_12", "slc_acct_roll_ever_24")
varInfo_Num <- describe(subset(datCredit_smp, select = varList_Num))

# - [SANITY CHECK] - Check for ovelaps between the categorical and numeric variables
check.fuse3 <- intersect(varList_Cat,varList_Num)
cat(if(identical(check.fuse3, character(0))) {'SAFE: No overlapping fields in the categorical and numerical variables'} else {
  'WARNING: Overlapping field(s) detected in the input space and the main credit dataset.'})


# --- 5.3 Missing value diagnostics & treatments
# - Diagnostics of missing values in the additional engineered "SLC" input space | If missingness > 50% missing remove variable
# Categorical variables
table(is.na(datCredit_smp$slc_pmnt_method)) %>% prop.table()              # missingness: 11.58% - keep variable
table(is.na(datCredit_smp$slc_acct_arr_dir_3)) %>% prop.table()           # missingness: 11.58% - keep variable
# Numerical variables
table(is.na(datCredit_smp$slc_past_due_amt)) %>% prop.table()             # missingness: 11.58% - keep variable
table(is.na(datCredit_smp$slc_days_excess)) %>% prop.table()              # missingness: 74.37% - discard variable
table(is.na(datCredit_smp$slc_acct_pre_lim_perc)) %>% prop.table()        # missingness: 11.58 - keep variable
table(is.na(datCredit_smp$slc_acct_prepaid_perc_dir_12)) %>% prop.table() # missingness: 11.68% - keep variable

# - Categorical variables
# [slc_pmnt_method] - Missing value indicators
varInfo_Cat$slc_pmnt_method
### RESULTS: [slc_pmnt_method] has 7 levels and 506737 observations have missing values. Bin the missing values with the already existent "unknown" bin.
# [TREATMENT] Binning "Unknown" values and missing values into one level
datCredit_smp[, slc_pmnt_method := 
                ifelse(is.na(slc_pmnt_method) | slc_pmnt_method == "" | slc_pmnt_method == "Unknown",
                       "MISSING_DATA", slc_pmnt_method)]
(varInfo_Cat$slc_pmnt_method <- describe(datCredit_smp$slc_pmnt_method))
### RESULTS: [slc_pmnt_method] has 7 levels and 0 observations have missing values.
# [TREATMENT] Apply factor transformation
datCredit_smp[,slc_pmnt_method:=factor(slc_pmnt_method)]
### RESULTS: Missing values imputed and facorisation applied to the levels of the variable.

# [slc_acct_arr_dir_3] - Missing value indicators
varInfo_Cat$slc_acct_arr_dir_3
### RESULTS: [slc_acct_arr_dir_3] has 4 levels and xx observations have missing values.
### [TREATMENT] Binning "N/A" values and missing values into one level
datCredit_smp[, slc_acct_arr_dir_3 := 
                ifelse(is.na(slc_acct_arr_dir_3) | slc_acct_arr_dir_3 == "" | slc_acct_arr_dir_3 == "N/A",
                       "MISSING_DATA", slc_acct_arr_dir_3)]
(varInfo_Cat$slc_acct_arr_dir_3 <- describe(datCredit_smp$slc_acct_arr_dir_3))
### RESULTS: [slc_acct_arr_dir_3] has 4 levels and 506737 missing values. Bin the missing values with the already existent "N/A" bin.
# [TREATMENT] Apply factor transformation
datCredit_smp[,slc_acct_arr_dir_3:=factor(slc_acct_arr_dir_3)]
### RESULTS: Missing values imputed and facorisation applied to the levels of the variable.

# - Numerical variables
# [slc_past_due_amt] - Missing value indicators
datCredit_smp[, value_ind_slc_past_due_amt := ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 0, 1)]
(varInfo_Cat$value_ind_slc_past_due_amt <- describe(datCredit_smp$value_ind_slc_past_due_amt))
### RESULTS: No missing values, variable created successfully.

# [slc_past_due_amt] - Missing value imputation
varInfo_Num$slc_past_due_amt
### RESULTS:    [slc_past_due_amt] has scale [0;2545590] and 506737 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 2718; note a VERY large single outlier in the right hand-side tail
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values (due to non-normality, i.e. high skewness)
datCredit_smp[, slc_past_due_amt_imputed_med := 
                ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_past_due_amt) | slc_past_due_amt == ""), slc_past_due_amt], na.rm=TRUE), slc_past_due_amt)]
(varInfo_Num$slc_past_due_amt_imputed_med <- describe(datCredit_smp$slc_past_due_amt_imputed_med))
### RESULTS: [slc_past_due_amt_imputed_med] has scale [0;2545589.79] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 2403
hist(datCredit_smp$slc_past_due_amt_imputed_med, breaks=500); skewness(datCredit_smp$slc_past_due_amt_imputed, na.rm = T); datCredit_smp[slc_past_due_amt_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS:    [slc_past_due_amt_imputed_med] is skewed to the right; Skewness = 27.67188; 91.32243% of variables have zero values.
### CONCLUISON: Use winsorisation on the upper 0.1% quantile to remove influence of outliers.
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_past_due_amt) | slc_past_due_amt == ""), slc_past_due_amt], probs = 0.999))
datCredit_smp[, slc_past_due_amt_imputed_med_wins := ifelse(slc_past_due_amt_imputed_med < wins_quant, slc_past_due_amt_imputed_med, wins_quant)]
(varInfo_Num$slc_past_due_amt_imputed_med_wins <- describe(datCredit_smp$slc_past_due_amt_imputed_med_wins))
### RESULTS: [slc_past_due_amt_imputed_wins] has scale [0;32027.06] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1012
### SUMMARY: Median imputation used to treat missing values.
###          Winsorisation applied to 99.9% quantile (value = 32027.06)
###          No further treatment necessary.

# [slc_acct_pre_lim_perc] - Missing value indicator
datCredit_smp[, value_ind_slc_acct_pre_lim_perc := ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 0, 1)]
(varInfo_Cat$value_ind_slc_acct_pre_lim_perc <- describe(datCredit_smp$value_ind_slc_acct_pre_lim_perc))
### RESTULS: Binning successful, no missing values in [value_ind_slc_acct_pre_lim_perc]

# [slc_acct_pre_lim_perc] - Missing value imputation
varInfo_Num$slc_acct_pre_lim_perc
### RESULTS:    [slc_acct_pre_lim_perc] has scale [0;1] and 506737 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.09247
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values (due to non-normality, i.e. high skewness)
datCredit_smp[, slc_acct_pre_lim_perc_imputed_med := 
                ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == ""), slc_acct_pre_lim_perc], na.rm=TRUE), slc_acct_pre_lim_perc)]
(varInfo_Num$slc_acct_pre_lim_perc_imputed_med <- describe(datCredit_smp$slc_acct_pre_lim_perc_imputed_med))
### RESTULS: [slc_acct_pre_lim_perc_imputed_med] has scale [0;0.711808] and 0 observations have missing values; with 0 at 50% quantile and 0 .005343at 75% quantile and mean of 0.08176.
hist(datCredit_smp$slc_acct_pre_lim_perc_imputed_med, breaks=500); skewness(datCredit_smp$slc_acct_pre_lim_perc_imputed_med, na.rm = T); datCredit_smp[slc_acct_pre_lim_perc_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS: Distribution is skewed to the right; Skewness = 3.099905; 69.4518% of variables have zero values.
### SUMMARY: Median imputation used to treat missing values.
###          No outliers present, the variable has a reasonable scale and hence no scaling is applied.
###          [slc_acct_pre_lim_perc_imputed] has scale [0;1] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.07596
###          No further treatment necessary

# [slc_acct_prepaid_perc_dir_12] - Missing value indicator
datCredit_smp[, value_ind_slc_acct_prepaid_perc_dir_12 := ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 0, 1)]
(varInfo_Cat$value_ind_slc_acct_prepaid_perc_dir_12 <- describe(datCredit_smp$value_ind_slc_acct_prepaid_perc_dir_12))
### RESTULS: Binning successful, no missing values in [value_ind_slc_acct_prepaid_perc_dir_12]

# [slc_acct_prepaid_perc_dir_12] - Missing value imputation
varSLC_Info_Num$slc_acct_prepaid_perc_dir_12
### RESULTS:    [slc_acct_prepaid_perc_dir_12] has scale [0;3609540717775.48] and 506737 observations have missing values; with 0 at 50% quantile and 0.07214 at 75% quantile and mean of 1216054
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performace spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed_med := 
                ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == ""), slc_acct_prepaid_perc_dir_12], na.rm=TRUE), slc_acct_prepaid_perc_dir_12)]
(varInfo_Num$slc_acct_prepaid_perc_dir_12_imputed_med <- describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med))
### RESULTS: [slc_acct_prepaid_perc_dir_12_imputed] has scale [0;3609540717775.48] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1075259
hist(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med, breaks=500); skewness(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med, na.rm = T); datCredit_smp[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 2066.594; 76.1361% of variables have zero values.
### CONCLUSION: Apply Winsorisation to the upper 1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values) - Winsorisation only applied to observations in default spell
###             NApply Winsorisation to the top 1%, this is since the scale is still very large when applying Winsorisation to the top 0.01%
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == ""), slc_acct_prepaid_perc_dir_12], probs = 0.99))
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed_med_wins := ifelse(slc_acct_prepaid_perc_dir_12_imputed_med < wins_quant, slc_acct_prepaid_perc_dir_12_imputed_med, wins_quant)]
(varInfo_Num$slc_acct_prepaid_perc_dir_12_imputed_wins <- describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med_wins))
### RESULTS: [slc_acct_prepaid_perc_dir_12_imputed_wins] has scale [0;87.56] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 1.455.
hist(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med_wins, breaks=500); skewness(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med_wins, na.rm = T); datCredit_smp[slc_acct_prepaid_perc_dir_12_imputed_med_wins==0,.N]/datCredit_smp[,.N]
### RESULTS: [slc_acct_prepaid_perc_dir_12_imputed_med_wins] is skewed to the right; Skewness = 8.708705; 76.1361% of variables have zero values.
### SUMMARY: Median imputation used to treat missing values.
###          Winsorisation applied to the 99% quantile (value = 87.556633).
###          No further treatment necessary.

# [slc_acct_roll_ever_24] - Missing value indicator
datCredit_smp[, value_ind_slc_acct_roll_ever_24 := ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "", 0, 1)]
(varInfo_Cat$value_ind_slc_acct_roll_ever_24 <- describe(datCredit_smp$value_ind_slc_acct_roll_ever_24))
### RESTULS: Binning successful, no missing values in [value_ind_slc_acct_roll_ever_24]

# [slc_acct_roll_ever_24] - Missing value imputation
varInfo_Num$slc_acct_roll_ever_24
### RESULTS:    [slc_acct_roll_ever_24] has scale [0;4] and 507101 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.4893
### CONCLUSION: Median imputation proceeds from the distribution of values which constitutes performance spells only (all associated values from default spells are excluded)
# [TREATMENT] Median imputation for missing values
datCredit_smp[, slc_acct_roll_ever_24_imputed_med := 
                ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "" | slc_acct_roll_ever_24 == "Unknown", 
                       median(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "") ,slc_acct_roll_ever_24], na.rm=TRUE), slc_acct_roll_ever_24)]
(varInfo_Cat$slc_acct_roll_ever_24_imputed_med <- describe(datCredit_smp$slc_acct_roll_ever_24_imputed_med))
### RESULTS: [slc_acct_roll_ever_24_imputed_med] has scale [0;3609540717775.48] and 0 observations have missing values; with 0 at 50% quantile and 0 at 75% quantile and mean of 0.4326
hist(datCredit_smp[,slc_acct_roll_ever_24_imputed_med], breaks=500); skewness(datCredit_smp$slc_acct_roll_ever_24_imputed_med, na.rm = T); datCredit_smp[slc_acct_roll_ever_24_imputed_med==0,.N]/datCredit_smp[,.N]
### RESULTS: Distribution is skewed to the right; Skewness = 2.501518; 82.32882% of variables have zero values.
### SUMMARY: Median imputation used to treat missing values.
###          No outliers present, since the variable represents counts and only has 4 levels.
###          No further treatment necessary

# - [InterestRate_Margin] (incorporating risk-based pricing info)
varInfo_Num$InterestRate_Margin; hist(datCredit_smp$InterestRate_Margin, breaks="FD")
datCredit_smp[is.na(InterestRate_Margin), .N] / datCredit_smp[,.N] * 100
### RESULTS:    Highly right-skewed distribution (as expected), with mean of -0.007 vs median of -0.008, 
###             bounded by [-0.02, 0.01] for 5%-95% percentiles; some negative outliers distort shape of distribution
### CONCLUSION: Use median imputation, given 0.53% missingness degree
datCredit_smp[, InterestRate_Margin_imputed_mean := 
                ifelse(is.na(InterestRate_Margin) | InterestRate_Margin == "", 
                       median(InterestRate_Margin, na.rm=TRUE), InterestRate_Margin)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(InterestRate_Margin_imputed_mean), .N] == 0) %?% 
       'SAFE: Treatment successful for [InterestRate_Margin_imputed_mean].\n' %:% 
       'ERROR: Treatment failed for [InterestRate_Margin_imputed_mean] \n' )
(varInfo_Num$InterestRate_Margin_imputed_mean <- describe(datCredit_smp$InterestRate_Margin_imputed_mean)); hist(datCredit_smp$InterestRate_Margin_imputed_mean, breaks="FD")
### RESULTS: Highly right-skewed distribution (as expected), with mean of -0.007 vs median of -0.008, bounded by [-0.02, 0.01] for 5%-95% percentiles.

# - Remove the variables that have missingness > 50%
suppressWarnings( datCredit_smp[, `:=`(value_ind_slc_days_excess = NULL, slc_days_excess = NULL, 
                                        value_ind_ccm_ute_lvl_40_cnt_24m = NULL, ccm_ute_lvl_40_cnt_24m = NULL,
                                        value_ind_ccm_worst_arrears_6m = NULL, ccm_worst_arrears_6m = NULL,
                                        value_ind_ccm_worst_arrears_24m = NULL, ccm_worst_arrears_24m = NULL)]); gc()




# --- 6. Feature Engineering: Binning and factorisation
# - Condense the payment group
datCredit_smp[, pmnt_method_grp := 
                case_when(slc_pmnt_method == "Debit Order FNB account" | slc_pmnt_method == "Debit Order other bank" ~ "Debit Order",
                          slc_pmnt_method == "Salary" | slc_pmnt_method == "Suspense" ~ "Salary/Suspense",
                          TRUE ~ slc_pmnt_method)]
# [SANITY CHECK] Check new feature for illogical values
cat( ( datCredit_smp[is.na(pmnt_method_grp), .N] == 0) %?% 
       'SAFE: New feature [pmnt_method_grp] has logical values.\n' %:% 
       'WARNING: New feature [pmnt_method_grp] has illogical values \n' )
(varInfo_Cat$pmnt_method_grp <- describe(datCredit_smp$pmnt_method_grp))
### RESULTS: Bins grouped logically such that each bin now has sufficient observations

# - [PerfSpellResol_Type_Hist]
varInfo_Cat$PerfSpellResol_Type_Hist
datCredit_smp$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
barplot(table(datCredit_smp$PerfSpellResol_Type_Hist))
### RESULTS [PerfSpellResol_Type_Hist] has 5 levels and no missing values. 12.52% of performance spells defaulted; 41.56% of performance spells have been right-censored.
# [TREATMENT] Apply factor transformation
datCredit_smp[,PerfSpellResol_Type_Hist:=factor(PerfSpellResol_Type_Hist)]
### RESULTS: No further treatment necessary

# - [LN_TPE]
varInfo_Cat$LN_TPE
datCredit_smp$LN_TPE %>% table() %>% prop.table()
### RESULTS: [LN_TPE] has 2 levels and no missing values. 89.1% of observations are in level 'CHL' and 10.9% are in level 'WHL'
# [TREATMENT] Apply factor transformation
datCredit_smp[,LN_TPE:=factor(LN_TPE)]
### RESULTS: ~ No further treatment necessary

# - Factorised [g0_Delinq] variable
datCredit_smp[,g0_Delinq_fac := as.factor(g0_Delinq)]
(varInfo_Cat$g0_Delinq_fac <- describe(datCredit_smp$g0_Delinq_fac))

# - Bin [InterestRate_Margin_imputed] | Binning the variable into three equally sized bins
datCredit_smp[, InterestRate_Margin_imputed_bin := factor(ntile(InterestRate_Margin_imputed_mean, n=3))]
(varInfo_Num$InterestRate_Margin_Imputed_Bin <- describe(datCredit_smp$InterestRate_Margin_imputed_bin))


# --- Saving description objects on the disk for fast disk-based retrieval
pack.ffdf(paste0(genObjPath,"Credit_Var_Cat_Descript"), varInfo_Cat)




# --- 7. Analysis and Treatments of Numeric variables
# - [Principal]
varInfo_Num$Principal
### RESULTS: [Principal] has scale [0.01;208477000], with 510000 at 50% quantile and 850000 at 75% quantile and mean of 646260.
hist(datCredit_smp$Principal, breaks='FD'); skewness(datCredit_smp$Principal, na.rm = T); datCredit_smp[Principal==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 15.13721; 0% of observations have zero values.
### CONCLUSION: Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Principal) | Principal == ""), Principal], probs = 0.999)) # 99.9% quantile = 4680480
datCredit_smp[, Principal_wins := ifelse(Principal < wins_quant, Principal, wins_quant)] # Note that this is applied to the top 1%, this is since the scale is still very large when applying Winsorisation to the top 0.01%/
(varInfo_Num$Principal_wins <- describe(datCredit_smp$Principal_wins))
### RESULTS: [Principal_wins] has scale [0.01;4680480.00] and 0 observations have missing values; with 510000 at 50% quantile and 850000 at 75% quantile and mean of 644358.
hist(datCredit_smp$Principal_wins, breaks='FD'); skewness(datCredit_smp$Principal_wins, na.rm = T); datCredit_smp[Principal_wins==0,.N]/datCredit_smp[,.N]
### RESULTS: [Principal_wins] is skewed to the right; Skewness = 2.026844; 0% of variables have zero values.
### SUMMARY: Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile (value = 4680480).
###          No further treatment necessary.

# - [Instalment]
varInfo_Num$Instalment
### RESULTS: [Instalment] has scale [0;16417873.80], with 4719.1 at 50% quantile and 7921.4 at 75% quantile and mean of 5947.
hist(datCredit_smp$Instalment, breaks ='FD'); skewness(datCredit_smp$Instalment, na.rm = T); datCredit_smp[Instalment==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 490.0671; 0.31% of observations have zero values.
### CONCLUSION: Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Instalment) | Instalment == ""), Instalment], probs = 0.999)) # 99.9% quantile = 43296.14
datCredit_smp[, Instalment_wins := ifelse(Instalment < wins_quant, Instalment, wins_quant)]
(varCredit_Info_Num$Instalment_wins <- describe(datCredit_smp$Instalment_wins))
### RESULTS: [Instalment_wins] has scale [0;43296.14] and 0 observations have missing values; with 4719.1 at 50% quantile and 7921.4 at 75% quantile and mean of 5870.
hist(datCredit_smp$Instalment_wins, breaks='FD'); skewness(datCredit_smp$Instalment_wins, na.rm = T); datCredit_smp[Instalment_wins==0,.N]/datCredit_smp[,.N]
### RESULTS: [Instalment_wins] is skewed to the right; Skewness = 2.005984; 0.03121732% of variables have zero values .
### SUMMARY: Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile (value = 43296.14).
###          No further treatment necessary.

# - [Balance]
varCredit_Info_Num$Balance
### RESULTS: [Balance] has scale [-146.71;41687282.40], with 365231.5 at 50% quantile and 692208.4 at 75% quantile and mean of 489701.
hist(datCredit_smp$Balance, breaks = 'FD'); skewness(datCredit_smp$Balance, na.rm = T); datCredit_smp[Balance==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 4.310693; 3.79% of variables have zero values.
### CONCLUSION: Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Balance) | Balance == ""), Balance], probs = 0.999)) # 99.9% quantile = 3844818
datCredit_smp[, Balance_wins := ifelse(Balance < wins_quant, Balance, wins_quant)]
(varInfo_Num$Balance_wins <- describe(datCredit_smp$Balance_wins))
### RESULTS: [Balance_wins] has scale [-146.71;3844818.17] and 0 observations have missing values; with 365231.5 at 50% quantile and 692208.4 at 75% quantile and mean of 488207.
hist(datCredit_smp$Balance_wins, breaks='FD'); skewness(datCredit_smp$Balance_wins, na.rm = T); datCredit_smp[Balance_wins==0,.N]/datCredit_smp[,.N]
### RESULTS: [Balance_wins] is skewed to the right; Skewness = 1.941194; 3.79% of variables have zero values.
### SUMMARY: Extreme outliers detected and Winsorisation applied to deal with the upper 0.01% quantile (value = 3844818 ).
###          Variable was scaled using min-max scaling with no shifting.
###          No further treatment necessary




# --- 8. Feature Engineering: ratio-type varaibles (Period-level)
# - [AgeToTerm]
varInfo_Num$AgeToTerm
### RESULTS: [AgeToTerm] has scale [0.002673797;22.400000000], with 0.28750 at 50% quantile and 0.53750 at 75% quantile and mean of 0.3647.
hist(datCredit_smp$AgeToTerm, breaks='FD'); skewness(datCredit_smp$AgeToTerm, na.rm = T); datCredit_smp[AgeToTerm==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 8.137409; 0% of variables have zero values.
### CONCLUSION: Apply Winsorisation to the upper 1% quantile to remove the influence of extreme outliers (which for this variable does not seem to be influential values)
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(AgeToTerm) | AgeToTerm == ""), AgeToTerm], probs = 0.999)) # 99.9% quantile = 2.483333
datCredit_smp[, AgeToTerm_wins := ifelse(AgeToTerm < wins_quant, AgeToTerm, wins_quant)]
(varInfo_Num$AgeToTerm_wins <- describe(datCredit_smp$AgeToTerm_wins))
### RESULTS: [AgeToTerm_wins] has scale [0.002673797;2.483333333] and 0 observations have missing values; with 0.28750 at 50% quantile and 0.53750 at 75% quantile and mean of 0.3628.
hist(datCredit_smp$AgeToTerm_wins, breaks='FD'); skewness(datCredit_smp$AgeToTerm_wins, na.rm = T); datCredit_smp[AgeToTerm_wins==0,.N]/datCredit_smp[,.N]
### RESULTS: [AgeToTerm_wins] is skewed to the right; Skewness = 1.489611; 0% of values have zero values.
### SUMMARY: Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile (value = 2.483333).
###          Winsorised variable was scaled using min-max scaling with no shifting.
###          No further treatment necessary

# - [BalanceToTerm]
varInfo_Num$BalanceToTerm
## RESULTS: [BalanceToTerm] has scale [-0.61129166667;173697.00999999998], with 1540.4410 at 50% quantile and 2910.6708 at 75% quantile and mean of 2070.
hist(datCredit_smp$BalanceToTerm, breaks='FD'); skewness(datCredit_smp$BalanceToTerm, na.rm = T); datCredit_smp[BalanceToTerm==0,.N]/datCredit_smp[,.N]
### RESULTS:    Distribution is skewed to the right; Skewness = 4.784716; 3.789549% of variables have zero values.
### CONCLUSION: Apply Winsorisation to the upper 0.1% quantile to remove the influence of extreme outliers (which for this variable does not seem to be influential values)
# [TREATMENT] Apply Winsorisation
(wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(BalanceToTerm) | BalanceToTerm == ""), BalanceToTerm], probs = 0.999)) # 99.9% quantile = 17659.6 
datCredit_smp[, BalanceToTerm_wins := ifelse(BalanceToTerm < wins_quant, BalanceToTerm, wins_quant)]
(varInfo_Num$BalanceToTerm_wins <- describe(datCredit_smp$BalanceToTerm_wins))
### RESULTS: [BalanceToTerm_wins] has scale [-0.61129166667;17659.6] and 0 observations have missing values; with 1540.4410 at 50% quantile and 2910.6708 at 75% quantile and mean of 2062.
hist(datCredit_smp$BalanceToTerm_wins, breaks='FD'); skewness(datCredit_smp$BalanceToTerm_wins, na.rm = T); datCredit_smp[BalanceToTerm_wins==0,.N]/datCredit_smp[,.N]
### RESULTS: [BalanceToTerm_wins] is skewed to the right; Skewness = 2.087682; 3.48% of values have zero values.
### SUMMARY: Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile (value = 17659.6).
###          Winsorised variable was scaled using min-max scaling with no shifting.
###          No further treatment necessary


# --- Save object to disk
pack.ffdf(paste0(genObjPath,"Credit_Var_Num_Descript"), varInfo_Num)


# --- Combine objects into a singe object and store to disk -> This object then describes all covariates
if (!exists('varInfo_Cat')) unpack.ffdf(paste0(genObjPath,"Credit_Var_Cat_Descript"), tempPath)
if (!exists('varInfo_Num')) unpack.ffdf(paste0(genObjPath,"Credit_Var_Num_Descript"), tempPath)

Covariate_Info <- c(varInfo_Cat,varInfo_Num)

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
# - Load in fused subsampled dataset (if not already in memory)
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final_smp_b"), tempPath)
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



