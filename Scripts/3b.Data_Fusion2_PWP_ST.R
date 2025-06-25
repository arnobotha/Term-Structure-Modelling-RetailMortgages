# ================================ DATA FUSION - PWPST ====================================
# Subsample the full credit dataset before fusing it to the wider "input space" data.
# Finally, analyse and fix basic deficiencies within the input space, followed by
# some basic feature engineering towards cultivating a comprehensive input space, suitable
# for Cox regression modelling (Prentice-Williams-Peterson [PWP] recurrent event subtype).
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB), Marcel Muller (MM), Roland Breedt (RB), 
#                   Bernard Scheepers (BS)
#
# DESCRIPTION:
# This script performs the following high-level tasks:
#   1) Subsample main dataset into a more manageable but still representative dataset
#   2) Fuse subsampled set with input space
#   3) Screen input space against missingness and apply appropriate treatments
#   4) Perform feature engineering (transforms, ratio-type)
# ---------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R

# -- Inputs:
#   - datInput.raw | raw input space imported in script 1
#   - datCredit_real | prepared credit data from script 3b
#   - various parameters set in the setup script 0
#   - datMV | prepared feature engineered macroeconomic data from script 2a
#
# -- Outputs:
#   - datCredit_final_c_PWPST | enriched credit dataset, fused with macroeconomic data and
#                        various input fields, with feature engineering applied
#   - datCredit_smp_PWPST   | subsampled dataset
# ---------------------------------------------------------------------------------------
# NOTE: This script predominantly comes from another project (Kasmeer).
# =======================================================================================





# ------ 1. Preliminaries

ptm <- proc.time() # for runtime calculations (ignore)

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)
doDescribe <- F # whether to run describe as part of production run or not

# - Field Names and other parameters
clusVar <- "LoanID" # for clustered random sampling on subject-level
clusVar_Spell <- "PerfSpell_Key" # for clustered random sampling on spell-level within a subject's history
timeVar <- "Date" # calendar time variable over which main event is observed & predicted, used in strata analysis
counter <- "Counter" # subject-level row index variable, used in strata analysis
# Optional field named (stratification)
stratifiers <- NA # First variable should be of type "date"
### NOTE: Assign stratifiers to NA if no stratifiers are desired. Good candidates include: "Event_Type", "LN_TPE", and "HasDefaulted_Ever"

# - Facet specification field names (for graphing purposes of the resolution rates)
resolType <- "PerfSpellResol_Type_Hist2" # Performance spell resolution types, used in strata analysis
resolType_Val <- "Defaulted" # Reference value in the performance spell resolution type field [resolType], used in strata analysis
### NOTE: Set to NA if not interested in creating additional facets for the performance spells using stopping time

# - Subsampling parameters
smp_size <- 90000 # fixed size of downsampled set in terms of the number of unique Loans
# Implied sampling fraction for downsampling step
smp_perc <- smp_size/length(unique(datCredit_real[!is.na(get(clusVar_Spell)),get(clusVar)]))
cat("NOTE: A fixed sample size of", comma(smp_size), "is", percent(smp_perc, accuracy=0.1), "of all loans.\n")
### RESULTS: 90k constitutes 14% of all loans

# - Resampling, stratification, and other general parameters
smp_frac <- 0.7 # sampling fraction for resampling scheme
minStrata_size <- 0 # Minimum strata size specified for subsample, used in strata analysis
timeDef_TFD <- F # Special logic during resampling for "TFD". 







# ------ 2. Clustered subsampling scheme with n-way stratified random sampling

# --- 1. Apply specified subsampling scheme onto data

# - Set seed for sampling 
set.seed(1, kind="Mersenne-Twister")

# - Training Key population
if (all(is.na(stratifiers))){ # No stratifiers
  # Get unique subject IDs or keys from the full dataset
  datKeys <- unique(datCredit_real[!is.na(get(clusVar_Spell)), mget(c(clusVar))])
  # Use simple random sampling to select at random some keys from which the training set will be populated 
  datKeys_sampled <- datKeys %>% slice_sample(prop=smp_perc) %>% as.data.table()
} else { # Stratifiers
  # Get unique subject IDs or keys from the full dataset
  datKeys <- unique(datCredit_real[!is.na(get(clusVar_Spell)), mget(c(clusVar, stratifiers))])
  # Use stratified random sampling to select at random some keys from which the training set will be populated 
  datKeys_sampled <- datKeys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
}

# - Obtain the associated loan records in creating the subsampled dataset
datCredit_smp <- copy(datCredit_real[get(clusVar) %in% datKeys_sampled[, get(clusVar),]])

# - Save intermediary snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1a"), datCredit_smp)
pack.ffdf(paste0(genPath,"creditdata_final_PWPST_keys"), datKeys_sampled)

# House keeping
rm(datCredit_real)


# --- 2. Strata analysis and subsequent exclusions to ensure adherence to specified minimum strata size
if (all(!is.na(stratifiers))){ # - Conditional loop for strata
  selectionVar_smp <- c(clusVar, timeVar, stratifiers)
  # - Test for exclusions given violations in the minimum strata size, as provided
  datStrata_smp_min <- datCredit_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
  cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datCredit_smp[get(counter)==1,.N], "(", 
      sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datCredit_smp[get(counter)==1,.N]*100), 
      "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
  
  # - Conditionally applying the exclusions
  if (sum(datStrata_smp_min[,Freq]) > 0){
    # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
    datCredit_smp_old_n <- datCredit_smp[,.N]
    # Computing the prior probabilities of the performance spell resolution outcomes
    datCredit_smp_prior <- datCredit_smp[get(timeVar)==timeVar_Perf_Min, get(resolType)] %>% table() %>% prop.table()
    datCredit_smp_prior <- datCredit_smp_prior[names(datCredit_smp_prior)[names(datCredit_smp_prior) == resolType_Val]][[1]]
    # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
    for (i in 1:datStrata_smp_min[,.N]){
      # Getting the type of class of each stratifier (used for building the ith condition)
      class_type <- sapply(datStrata_smp_min[,1:length(stratifiers)], function(x) {class(x[[1]])}) 
      
      excCond <- datStrata_smp_min[i,1:length(stratifiers)] # Getting the values of the ith minimum strata
      excCond <- data.table(Stratifier = colnames(excCond), # Building a dataset
                            Value = unname(t(excCond)), # Ensure that the column name is Value instead of Value.V1
                            Class = class_type)
      excCond[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), 
                                   paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
      excCond[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
      # Compiling the ith condition
      excCond2 <- parse(text = paste0(paste0(excCond$Condition, collapse = ""), counter,"==1"))
      # Add the excluded subject key to our list
      if (i==1) dat_keys_exc <- as.vector(datCredit_smp[eval(excCond2), get(clusVar)]) # set new list
      else dat_keys_exc <- c(dat_keys_exc, as.vector(datCredit_smp[eval(excCond2), get(clusVar)]) ) # append to existing list
    }
    
    # Applying the exclusions to the subsampled dataset
    datCredit_smp <- copy(datCredit_smp[!(get(clusVar) %in% dat_keys_exc),])
    
    cat(datCredit_smp_old_n-datCredit_smp[,.N], " observations removed (", 
        sprintf("%.4f", (datCredit_smp_old_n-datCredit_smp[,.N])/datCredit_smp_old_n*100), "% ) \n",
        "Prior probability = ", sprintf("%.4f", datCredit_smp_prior*100), "% comapred to ", 
        sprintf("%.4f", (datCredit_smp[get(timeVar)==timeVar_Perf_Min, get(resolType)] %>% table() %>% prop.table())[[2]]*100), "%")
  }
  # [SANITY CHECK] Are there still violations in minimum strata sizes?
  datStrata_smp_min <- datCredit_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
  if (NROW(datStrata_smp_min) > 0) cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datCredit_smp[get(counter)==1,.N], "(", 
                                       sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datCredit_smp[get(counter)==1,.N]*100), 
                                       "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
  # - Cleanup
  suppressWarnings( rm(datStrata_smp_min, datStrata_smp_min, datCredit_smp_old_n, datCredit_smp_prior, 
                       dat_keys_exc, class_type, excCond, excCond2))
}




# ------- 3. Fuse the input space with the subsampled prepared dataset

# - Confirm that required data objects are loaded into memory
if (!exists('datInput.raw')) unpack.ffdf(paste0(genPath,"creditdata_input1"), tempPath)
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1a"), tempPath)

# [SANITY CHECK] Prevalence of overlapping fields in the input space and the main credit dataset
# Find intersection between fields in input space and those perhaps already in the main credit dataset
overlap_flds <- intersect(colnames(datCredit_smp), colnames(datInput.raw))
check.fuse1 <- length(overlap_flds) == 0 # FALSE; duplicate columns exists.
cat(check.fuse1 %?% 'SAFE: No overlapping fields in the input space and the main credit dataset' %:%
      'WARNING: Overlapping field(s) detected in the input space and the main credit dataset.')
# Conditional reporting
if (check.fuse1 == 0) {cat('NOTE: The following fields overlap: ', overlap_flds,"\n",sep="\t")}
### RESULTS: slc_past_due_amt overlap

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
sum(is.na(data_grain_check$LoanID)); gc()
# the data grain is broken in the cases where a Loan_ID does not exist - we are not interested in these accounts in any case
# - Merge on LoanID and Date by performing a left-join
datCredit_smp <- merge(datCredit_smp, datInput.raw, by=c("Date", "LoanID"), all.x=T); gc()
# - Check the data grain
NROW(data_grain_check_merge <- datCredit_smp[, list(Freq = .N), by=list(LoanID, Date)][Freq>1,])==0
# success, the data grain check is passed

# - Save intermediary snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1b"), datCredit_smp)
# - Clean-up
rm(datInput.raw, data_grain_check, data_grain_check_merge); gc()




# ------- 4. Feature engineering for modelling purposes

# - Confirm that required data objects are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1a"), tempPath)


# --- 1. Missing value indicators for the input space variables
# NOTE: There are a lot of missing values for these variables because of system changes, data migrations, etc.
datCredit_smp[, value_ind_slc_pmnt_method := ifelse(is.na(slc_pmnt_method) | slc_pmnt_method == "", 1, 0)]
datCredit_smp[, value_ind_slc_days_excess := ifelse(is.na(slc_days_excess) | slc_days_excess == "", 1, 0)]
datCredit_smp[, value_ind_slc_acct_pre_lim_perc := ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 1, 0)]
datCredit_smp[, value_ind_slc_acct_roll_ever_24 := ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "", 1, 0)]
datCredit_smp[, value_ind_slc_acct_arr_dir_3 := ifelse(is.na(slc_acct_arr_dir_3) | slc_acct_arr_dir_3 == "", 1, 0)]
datCredit_smp[, value_ind_slc_acct_prepaid_perc_dir_12 := ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 1, 0)]
datCredit_smp[, value_ind_slc_past_due_amt := ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 1, 0)]

# - Check the missingness of the variables
# If they are more than 50% missing - remove
table(datCredit_smp$value_ind_slc_pmnt_method) %>% prop.table()              # missingness: 11.8% - keep the variable (categorical)
table(datCredit_smp$value_ind_slc_days_excess) %>% prop.table()              # missingness: 74.9% - discard the variable
table(datCredit_smp$value_ind_slc_acct_pre_lim_perc) %>% prop.table()        # missingness: 11.8% - keep the variable (numeric) 
table(datCredit_smp$value_ind_slc_acct_roll_ever_24) %>% prop.table()        # missingness: 11.8% - keep the variable (numeric + delinquency theme)     
table(datCredit_smp$value_ind_slc_acct_arr_dir_3) %>% prop.table()           # missingness: 11.8% - keep the variable (categorical + delinquency theme)        
table(datCredit_smp$value_ind_slc_acct_prepaid_perc_dir_12) %>% prop.table() # missingness: 11.8% - keep the variable (numeric)

# - Remove the variables that have missingness > 50%
suppressWarnings( datCredit_smp[, `:=`(value_ind_slc_days_excess = NULL, slc_days_excess = NULL)]); gc()



# --- 2. Missing value treatment (categorical variables)
# Treatment: create "missing"-bin for all N/A values

# - Payment method
if (doDescribe) describe(datCredit_smp$slc_pmnt_method)
# Merge with existing "Unknown" bin or empty values
datCredit_smp[, slc_pmnt_method := 
                ifelse(is.na(slc_pmnt_method) | slc_pmnt_method == "" | slc_pmnt_method == "Unknown",
                       "MISSING_DATA", slc_pmnt_method)]
# [SANITY CHECK] Confirm treatment success
cat( (sum(datCredit_smp$slc_pmnt_method == "" | is.na(datCredit_smp$slc_pmnt_method) | 
            datCredit_smp$slc_pmnt_method == "Unknown") == 0) %?% 
       'SAFE: Treatment successful for [slc_pmnt_method].\n' %:% 'ERROR: Treatment failed for [slc_pmnt_method] \n' )
### RESULTS: Treatment for missingness was successful


# - Account-level arrears direction vs three months ago
if (doDescribe) describe(datCredit_smp$slc_acct_arr_dir_3)
# Merge with existing "N/A" bin or empty values
datCredit_smp[, slc_acct_arr_dir_3 := 
                ifelse(is.na(slc_acct_arr_dir_3) | slc_acct_arr_dir_3 == "" | slc_acct_arr_dir_3 == "N/A", 
                       "MISSING_DATA", slc_acct_arr_dir_3)]
# [SANITY CHECK] Confirm treatment success
cat( ( sum(datCredit_smp$slc_acct_arr_dir_3 == "" | is.na(datCredit_smp$slc_acct_arr_dir_3) |
             datCredit_smp$slc_acct_arr_dir_3 == "N/A") == 0) %?% 
       'SAFE: Treatment successful for [slc_acct_arr_dir_3].\n' %:% 'ERROR: Treatment failed for [slc_acct_arr_dir_3] \n' )
### RESULTS: Treatment for missingness was successful

# - Save to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath, "creditdata_final_PWPST_smp1b"), datCredit_smp); gc()



# --- 3. Missing value treatment (numeric variables)
# Analyse whether to use mean or median value imputation

# - Confirm that required data objects are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1b"), tempPath)

# - Prepaid/available funds to limit
if (doDescribe) describe(datCredit_smp$slc_acct_pre_lim_perc); hist(datCredit_smp$slc_acct_pre_lim_perc, breaks='FD')
datCredit_smp[is.na(slc_acct_pre_lim_perc), .N] / datCredit_smp[,.N] * 100
### RESULTS: Highly right-skewed distribution, with mean of ~0.098 vs median of 0,
# bounded by [0, 0.82] for 5%-95% percentiles; no outliers, other than at 0 and 1
# Use median imputation, given 11.8% missingness degree, trading off the minor distributional distortion as a result
datCredit_smp[, slc_acct_pre_lim_perc_imputed_med := 
                ifelse(is.na(slc_acct_pre_lim_perc) | slc_acct_pre_lim_perc == "", 
                       median(slc_acct_pre_lim_perc, na.rm=TRUE), slc_acct_pre_lim_perc)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(slc_acct_pre_lim_perc_imputed_med), .N ] == 0) %?% 
       'SAFE: Treatment successful for [slc_acct_pre_lim_perc_imputed_med].\n' %:% 
       'ERROR: Treatment failed for [slc_acct_pre_lim_perc_imputed_med] \n' )
if (doDescribe) describe(datCredit_smp$slc_acct_pre_lim_perc_imputed_med); 
hist(datCredit_smp$slc_acct_pre_lim_perc_imputed_med, breaks='FD')
### RESULTS: Imputation successful, with mean of 0.086 vs median of ,
# bounded by [0, 0.74] for 5%-95% percentiles; no outliers, other than at 0 and 1

# - Number of times an account was in arrears over last 24 months
if (doDescribe) describe(datCredit_smp$slc_acct_roll_ever_24); hist(datCredit_smp$slc_acct_roll_ever_24, breaks='FD')
datCredit_smp[is.na(slc_acct_roll_ever_24), .N] / datCredit_smp[,.N] * 100
### RESULTS: Highly right-skewed distribution with mean of 0.3032, though discrete values with 84% of data having 0-value.
# Use mean imputation, given 11.8% missingness degree, trading off the minor distributional distortion as a result
datCredit_smp[, slc_acct_roll_ever_24_imputed_mean := 
                ifelse(is.na(slc_acct_roll_ever_24) | slc_acct_roll_ever_24 == "", 
                       mean(slc_acct_roll_ever_24, na.rm=TRUE), slc_acct_roll_ever_24)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(slc_acct_roll_ever_24_imputed_mean), .N ] == 0) %?% 
       'SAFE: Treatment successful for [slc_acct_roll_ever_24_imputed_mean].\n' %:% 
       'ERROR: Treatment failed for [slc_acct_roll_ever_24_imputed_mean] \n' )
if (doDescribe) describe(datCredit_smp$slc_acct_roll_ever_24_imputed_mean); 
hist(datCredit_smp$slc_acct_roll_ever_24_imputed_mean, breaks='FD')
### RESULTS: Imputation successful, categorical variable now has 6 distinct classes, with majority having 0-value, while
# the imputed cases (value of 0.30) being the second most prevalent.

# - Percentage-valued direction of prepaid/available funds - current compared to 12 months ago
if (doDescribe) describe(datCredit_smp$slc_acct_prepaid_perc_dir_12); 
hist(datCredit_smp[slc_acct_prepaid_perc_dir_12<=5, slc_acct_prepaid_perc_dir_12])
datCredit_smp[is.na(slc_acct_prepaid_perc_dir_12), .N] / datCredit_smp[,.N] * 100
### RESULTS: Highly right-skewed distribution, with mean of ~13m vs median of 0, 
# bounded by [0, 3.62] for 5%-95% percentiles; some very large outliers
# Use median imputation, given 11.8% missingness degree, trading off the minor distributional distortion as a result
datCredit_smp[, slc_acct_prepaid_perc_dir_12_imputed_med := 
                ifelse(is.na(slc_acct_prepaid_perc_dir_12) | slc_acct_prepaid_perc_dir_12 == "", 
                       median(slc_acct_prepaid_perc_dir_12, na.rm=TRUE), slc_acct_prepaid_perc_dir_12)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(slc_acct_prepaid_perc_dir_12_imputed_med), .N] == 0) %?% 
       'SAFE: Treatment successful for [slc_acct_prepaid_perc_dir_12_imputed_med].\n' %:% 
       'ERROR: Treatment failed for [slc_acct_prepaid_perc_dir_12_imputed_med] \n' )
if (doDescribe) describe(datCredit_smp$slc_acct_prepaid_perc_dir_12_imputed_med); 
hist(datCredit_smp[slc_acct_prepaid_perc_dir_12_imputed_med<=5, slc_acct_prepaid_perc_dir_12_imputed_med])
### RESULTS: Imputation successful, with mean of ~12m vs median of 0,
# bounded by [0, 3.07] for 5%-95% percentiles; extreme outliers


# - Amount by which the account is overdue at the associated reporting date
if (doDescribe) describe(datCredit_smp$slc_past_due_amt); hist(datCredit_smp$slc_past_due_amt, breaks='FD')
datCredit_smp[is.na(slc_past_due_amt), .N] / datCredit_smp[,.N] * 100
### RESULTS: Highly right-skewed distribution, with mean of 329.4 vs median of 0, 
# bounded by [0, 171.4] for 5%-95% percentiles; some very large outliers
### MM: Scope for eventual extreme value treatment if those outliers are correct; or use winsorized mean
# Use median imputation, given 11.8% missingness degree, trading off the minor distributional distortion as a result
datCredit_smp[, slc_past_due_amt_imputed_med := 
                ifelse(is.na(slc_past_due_amt) | slc_past_due_amt == "", 
                       median(slc_past_due_amt, na.rm=TRUE), slc_past_due_amt)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(slc_past_due_amt_imputed_med), .N] == 0) %?% 
       'SAFE: Treatment successful for [slc_past_due_amt_imputed_med].\n' %:% 
       'ERROR: Treatment failed for [slc_past_due_amt_imputed_med] \n' )
if (doDescribe) describe(datCredit_smp$slc_past_due_amt_imputed_med); 
hist(datCredit_smp$slc_past_due_amt_imputed_med[datCredit_smp$slc_past_due_amt_imputed_med>0], breaks='FD')
### RESULTS: Imputation successful, with mean of 290.6 vs median of 0,
# bounded by [0, 0] for 5%-95% percentiles; extreme outliers

# - InterestRate_Margin (incorporating risk-based pricing info)
if (doDescribe) describe(datCredit_smp$InterestRate_Margin); hist(datCredit_smp$InterestRate_Margin, breaks="FD")
datCredit_smp[is.na(InterestRate_Margin), .N] / datCredit_smp[,.N] * 100
### RESULTS: Highly right-skewed distribution (as expected), with mean of -0.007 vs median of -0.008, 
# bounded by [-0.02, 0.01] for 5%-95% percentiles; some negative outliers distort shape of distribution
# Use median imputation for any possible missingness
datCredit_smp[, InterestRate_Margin_imputed_mean := 
                ifelse(is.na(InterestRate_Margin) | InterestRate_Margin == "", 
                       median(InterestRate_Margin, na.rm=TRUE), InterestRate_Margin)]
# [SANITY CHECK] Confirm treatment success
cat( ( datCredit_smp[is.na(InterestRate_Margin_imputed_mean), .N] == 0) %?% 
       'SAFE: Treatment successful for [InterestRate_Margin_imputed_mean].\n' %:% 
       'ERROR: Treatment failed for [InterestRate_Margin_imputed_mean] \n' )
if (doDescribe) describe(datCredit_smp$InterestRate_Margin_imputed_mean);
hist(datCredit_smp$InterestRate_Margin_imputed_mean, breaks="FD")
### RESULTS: Imputation successful, with mean of -0.007 vs median of -0.008,
# bounded by [-0.02, 0.01] for 5%-95% percentiles; some negative outliers distort shape of distribution

# - Save to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath, "creditdata_final_PWPST_smp1c"), datCredit_smp); gc()



# --- 4. Feature Engineering: ratio-type variables (Period-level)

# - Confirm that required data objects are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1c"), tempPath)

# - Loan age to loan term
datCredit_smp[, AgeToTerm := Age_Adj/Term] # where the loan is in its lifetime
# [SANITY CHECK] Check new feature for illogical values
cat( ( datCredit_smp[is.na(AgeToTerm), .N] == 0) %?% 
       'SAFE: New feature [AgeToTerm] has logical values.\n' %:% 
       'WARNING: New feature [AgeToTerm] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$AgteToTerm); hist(datCredit_smp[AgeToTerm<2, AgeToTerm], breaks='FD')
### RESULTS: Highly right-skewed distribution as expected, with mean of 0.37 vs median of 0.29,
# bounded by [0.03, 0.9] for 5%-95% percentiles; some large outliers (max: 224)

# - Balance to loan term | how much is still outstanding compared to Principal/Limit
datCredit_smp[, BalanceToPrincipal := Balance/Principal]
# [SANITY CHECK] Check new feature for illogical values
cat( ( datCredit_smp[is.na(BalanceToPrincipal), .N] == 0) %?% 
       'SAFE: New feature [BalanceToPrincipal] has logical values.\n' %:% 
       'WARNING: New feature [BalanceToPrincipal] has illogical values \n' )
# distributional analysis
if (doDescribe) describe(datCredit_smp$BalanceToPrincipal); hist(datCredit_smp$BalanceToPrincipal, breaks='FD')
### RESULTS: Highly left-skewed distribution, with mean of 0.7 vs median of 0.85,
# bounded by [~0, 1] for 5%-95% percentiles; no outliers



# --- 5. Feature Engineering: Binning and factorisation
# - Condense the payment group
datCredit_smp[, pmnt_method_grp := 
                case_when(slc_pmnt_method == "Debit Order FNB account" | slc_pmnt_method == "Debit Order other bank" ~ "Debit Order",
                          slc_pmnt_method == "Salary" | slc_pmnt_method == "Suspense" ~ "Salary/Suspense",
                          TRUE ~ slc_pmnt_method)]
# [SANITY CHECK] Check new feature for illogical values
cat((datCredit_smp[is.na(pmnt_method_grp), .N] == 0) %?% 
      'SAFE: New feature [pmnt_method_grp] has logical values.\n' %:% 
      'WARNING: New feature [pmnt_method_grp] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$pmnt_method_grp)
### RESULTS: Bins grouped logically such that each bin now has sufficient observations, with proportions:
# Debit Order: 68%; MISSING_DATA: 13%; Salary/Suspense: 6%; Statement: 12%
### CONCLUSION: Given the greater utility of this newly-binned variable, rather use this variable than [scl_pmnt_method]

# - Factorised [g0_Delinq] variable
datCredit_smp[,g0_Delinq_fac := as.factor(g0_Delinq)]
if (doDescribe) describe(datCredit_smp$g0_Delinq_fac)
### RESULTS: Proportion of observations in each g0 bucket:
# 0: 94%; 1: 5%; 2: 0.6%; 3: 0.4%
# [SANITY CHECK] Check new feature for illogical values
cat((anyNA(datCredit_smp$g0_Delinq_fac)) %?% 'WARNING: New feature [g0_Delinq_fac] has missing values. \n' %:%
      'SAFE: New feature [g0_Delinq_fac] has no missing values. \n')
### RESULTS: [g0_Delinq_fac] created without any missingness

# - Bin [InterestRate_Margin_imputed] | Binning the variable into three equally sized bins
datCredit_smp[, InterestRate_Margin_imputed_bin := factor(ntile(InterestRate_Margin_imputed_mean, n=3))]
if (doDescribe) describe(datCredit_smp$InterestRate_Margin_imputed_bin)
### RESULTS: Binning was successful into three equal sized bins
# [SANITY CHECK] Check new feature for illogical values
cat((anyNA(datCredit_smp$InterestRate_Margin_imputed_bin)) %?% 'WARNING: New feature [InterestRate_Margin_imputed_bin] has missing values. \n' %:%
      'SAFE: New feature [InterestRate_Margin_imputed_bin] has no missing values. \n')
### RESULTS: [InterestRate_Margin_imputed_bin] created without any missingness

# - Bin [PerfSpell_Num] based on previous analysis (script 4a(i)) towards grouping later spells together
# [SANITY CHECK] Check new feature for illogical values
datCredit_smp[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]
cat( (all.equal(datCredit_smp[PerfSpell_Num < 4,PerfSpell_Num], datCredit_smp[PerfSpell_Grp < 4,PerfSpell_Grp]) & 
        datCredit_smp[!is.na(PerfSpell_Key) & is.na(PerfSpell_Grp), .N] == 0) %?%
       'SAFE: New feature [PerfSpell_Grp] has no missing values and is binned as intended. \n' %:%
       'WARNING: New feature [PerfSpell_Grp] either has missing valeus or its binning failed. \n ' )

# - Create indicator variable for specific arrears direction
datCredit_smp[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]

# - Create binned version of performance spell number
datCredit_smp[, PerfSpell_Num_binned := ifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]


# --- 6. Feature Engineering: Inflating time-sensitive monetary variables to the latest date
# - Confirm that required data objects are loaded into memory
if (!exists('datMV')) unpack.ffdf(paste0(genPath,"datMV"), tempPath)

# - Getting a range of inflation factors for each date in the sampling window
date_range <- ceiling_date(unique(datCredit_smp$Date), unit="month")-days(1)
datInflation <- data.table(Date=date_range)
datInflation[,Inf_Factor:=adjInflation_MV(datMacro=datMV, time="Date", Inflation_Growth="M_Inflation_Growth", g_start=Date, g_stop = date_range[length(date_range)]), by=Date]
datCredit_smp <- merge(datCredit_smp, datInflation, all.x=T, by="Date")
# [SANITY CHECK] Are inflation factors missing?
cat((anyNA(datCredit_smp$Inf_Factor)) %?% paste0('WARNING: Inflation factor(s) is(are) missing for ', unique(datCredit_smp[is.na(Inf_Factor),Date]), '. \n') %:%
      'SAFE: Inflation factors created successfully. \n')
### RESULTS: [Inf_Factor] variables  missing for 2022-12-31.

# - Deflate the relevant variables using the pre-calculated inflation factors
datCredit_smp[, Principal_Real := Principal*Inf_Factor]
datCredit_smp[, Balance_Real := Balance*Inf_Factor]
datCredit_smp[, Instalment_Real := Instalment*Inf_Factor]

# [SANITY CHECK] Missingness in any new variables?
cat( (all(anyNA(datCredit_smp$Principal_Real), anyNA(datCredit_smp$Balance_Real), anyNA(datCredit_smp$Instalment_Real)))
     %?% paste0('WARNING: Some values of [Principal_Real], [Balance_Real], and/or [Instalment_Real] not created successfully. \n') %:%
       'SAFE: Variables inflated successfully. \n')
if (doDescribe) describe(datCredit_smp$Principal_Real); hist(datCredit_smp$Principal_Real[datCredit_smp$Principal_Real<5000000], breaks="FD")
if (doDescribe) describe(datCredit_smp$Balance_Real); hist(datCredit_smp$Balance_Real[datCredit_smp$Balance_Real< 5000000], breaks="FD")
if (doDescribe) describe(datCredit_smp$Instalment_Real); hist(datCredit_smp$Instalment_Real[datCredit_smp$Instalment_Real<25000], breaks="FD")
### RESULTS:  Some values of [Principal_Real], [Balance_Real], and/or [Instalment_Real] not created successfully, refer to line 464.
# [Principal_Real]# Highly right-skewed distribution, with mean of 8.95m vs median of 7.9m
#                 bounded by [127k, 2.25m] for 5%-95% percentiles; severe outliers to the right: 95.6m
# [Balance_Real]: Highly right-skewed distribution, with mean of 6.89m vs median of 5.93m
#                 bounded by [71, 1.5m] for 5%-95% percentiles; severe outliers to the right: 94m
# [Instalment_Real]: Highly right-skewed distribution, with mean of 8.4k vs median of 7.3k
#                 bounded by [497, 21k] for 5%-95% percentiles; severe outliers to the right: 20m; left: 0

# - Save to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath, "creditdata_final_PWPST_smp1d"), datCredit_smp); gc()

# - Clean up
rm(datMV, date_range, datInflation)



# --- 8. Featuring Engineering: Portfolio-level information

# - Confirm that required data objects are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1d"), tempPath)

# - Pre default delinquency rate
#Note: Creating an aggregated dataset with which to fuse to the full dataset
dat_g0_Delinq_Aggr <- data.table(datCredit_smp[DefaultStatus1==0, list(sum(g0_Delinq>0, na.rm=T)/.N), by=list(Date)])
colnames(dat_g0_Delinq_Aggr) <- c("Date", "g0_Delinq_Any_Aggr_Prop")
# Applying various lags
lags <- c(1,2,3,4,5,6,9,12) # Lags
ColNames <- colnames(dat_g0_Delinq_Aggr)[-1] # Names of the columns
for (i in seq_along(lags)){ # Looping over the specified lags and applying each to each of the specified columns
  for (j in seq_along(ColNames)){
    # Impute NA's with the non lagged value
    dat_g0_Delinq_Aggr[, (paste0(ColNames[j],"_Lag_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"),get(ColNames[j]))]
  }
}
# [Sanity Check] Check for any missing values before merging the dat_g0_Delinq_Aggr dataset to datCredit_smp
cat((anyNA(dat_g0_Delinq_Aggr)) %?% 'WARNING: One of the new [g0_Delinq_Any_Aggr_Prop] features has missing values. \n' %:%
      'SAFE: New [g0_Delinq_Any_Aggr_Prop] features created sucessfully without any missing values. \n')
### RESULTS: [g0_Delinq_Any_Aggr_Prop] variables created successfully without any missingness
# Fusing the aggregated variable with its various lags to the full dataset
datCredit_smp <- merge(datCredit_smp, dat_g0_Delinq_Aggr, by="Date", all.x=T)
# [SANITY CHECK] Check new feature for illogical values
cat( ( sum(datCredit_smp[DefaultStatus1==0, sum(g0_Delinq_Any_Aggr_Prop + sum(g0_Delinq==0)/.N, na.rm=T), by=Date][,2]) ==
         sum(datCredit_smp[DefaultStatus1==0,.N,by=Date][,2]) & (sum(is.na(datCredit_smp$g0_Delinq_Any_Aggr_Prop))==0)) %?% 
       'SAFE: New feature [g0_Delinq_Any_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [g0_Delinq_Any_Aggr_Prop] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$g0_Delinq_Any_Aggr_Prop); plot(unique(datCredit_smp$g0_Delinq_Any_Aggr_Prop), type="b")
### RESULTS: Variable has a logical trend, with mean of 0.059 vs median of 0.049, 
# bounded by [0.038, 0.12] for 5%-95% percentiles; no large outliers
# [SANITY CHECK] Check new feature for missingness after fusion
cat((anyNA(datCredit_smp$g0_Delinq_Any_Aggr_Prop)) %?% 'WARNING: New feature [g0_Delinq_Any_Aggr_Prop] has missing values. \n' %:%
      'SAFE: New feature [g0_Delinq_Any_Aggr_Prop] has no missing values. \n')
### RESULTS: [g0_Delinq_Any_Aggr_Prop] created without any missingness


# - Average pre-default delinquency level
datCredit_smp[,g0_Delinq_Ave:=mean(ifelse(DefaultStatus1==0,g0_Delinq,0), na.rm=T), by=Date]
# [SANITY CHECK] Check new feature for illogical values
cat( (sum(datCredit_smp[, sum(is.na(g0_Delinq_Ave)), by=Date][,2])==0) %?% 
       'SAFE: New feature [g0_Delinq_Ave] has logical values.\n' %:% 
       'WARNING: New feature [g0_Delinq_Ave] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$g0_Delinq_Ave); hist(datCredit_smp$g0_Delinq_Ave, breaks="FD")
### RESULTS: Follows a logical trend, with mean of 0.06397 vs median of 0.05269,
# bounded by [0.04357, 0.12998] for 5%-95% percentiles; no outliers


# - Create a lagged version of the aggregated default rate created in script 2f
dat_DefaultRate <- datCredit_smp[!duplicated(Date),list(Date,DefaultStatus1_Aggr_Prop)]
# Applying various lags
lags <- c(1,2,3,4,5,6,9,12) # Lags
ColNames <- colnames(dat_DefaultRate)[-1] # Names of the columns
for (i in seq_along(lags)){ # Looping over the specified lags and applying each to each of the specified columns
  for (j in seq_along(ColNames)){
    dat_DefaultRate[, (paste0(ColNames[j],"_Lag_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"),get(ColNames[j]))] # Impute NA's with the non lagged value
  }
}
# [Sanity Check] Check for any missing values before merging the dat_DefaultRate dataset to datCredit_smp
cat((anyNA(dat_DefaultRate)) %?% 'WARNING: One of the new DefaultRate features has missing values. \n' %:%
      'SAFE: New DefaultRate features created sucessfully without any missing values. \n')
### RESULTS: No missingness, continue with merge
# Fusing the various lagged versions of the DefaultRate to the full dataset
datCredit_smp <- merge(datCredit_smp, dat_DefaultRate[,-"DefaultStatus1_Aggr_Prop"], by="Date", all.x=T)
# [Sanity Check] Check if merge was successful by checking for missingness in the 12-month lagged version of the default rate
cat((anyNA(datCredit_smp$DefaultStatus1_Aggr_Prop_Lag_12)) %?% 'WARNING: Merge unsuccessful, NA values present. \n' %:%
      'SAFE: Merge successful, no NA values present. \n')
if (doDescribe) describe(datCredit_smp$DefaultStatus1_Aggr_Prop_Lag_12); hist(datCredit_smp$DefaultStatus1_Aggr_Prop_Lag_12, breaks="FD")
### RESULTS: No missingness after merge, merge successful. Original DefaultRate_12 has mean of 0.05 and median of 0.049; 
# bounded by [0.028, 0.084] for 5%-95% percentiles; no outliers


# - Ratio type variables (portfolio-level) during performance spells
# (Total) Arrears to (Total) Balance; (Total) Instalments to (Total) Balance
# NOTE: These portfolio-level aggregated variables are engineered to capture/ aggregate information only for accounts that are in a performance spell
# The resulting aggregated dataset can be fused to the full dataset
dat_Aggr <- data.table(datCredit_smp[DefaultStatus1==0, list(sum(Arrears, na.rm=T)/sum(Balance, na.rm=T)), by=list(Date)], # [ArrearsToBalance_Aggr]
                       datCredit_smp[DefaultStatus1==0, list(sum(Instalment, na.rm=T)/sum(Balance)), by=list(Date)][,2]) # [InstalmentToBalance_Aggr]
colnames(dat_Aggr) <- c("Date", "ArrearsToBalance_Aggr_Prop", "InstalmentToBalance_Aggr_Prop")
# Fusing the aggregated dataset to the full dataset
datCredit_smp <- merge(datCredit_smp, dat_Aggr, by="Date", all.x=T)
# [SANITY CHECK] Check new feature for illogical values
cat( (sum(datCredit_smp[, sum(is.na(ArrearsToBalance_Aggr_Prop)), by=Date][,2])==0) %?% 
       'SAFE: New feature [ArrearsToBalance_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [ArrearsToBalance_Aggr_Prop] has illogical values \n' )
cat( (sum(datCredit_smp[, sum(is.na(InstalmentToBalance_Aggr_Prop)), by=Date][,2])==0) %?% 
       'SAFE: New feature [InstalmentToBalance_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [InstalmentToBalance_Aggr_Prop] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$InstalmentToBalance_Aggr_Prop); 
plot(unique(datCredit_smp$Date),unique(datCredit_smp$InstalmentToBalance_Aggr_Prop), type="b")
if (doDescribe) describe(datCredit_smp$ArrearsToBalance_Aggr_Prop); 
plot(unique(datCredit_smp$ArrearsToBalance_Aggr_Prop), type="b")
### RESULTS [InstalmentToBalance_Aggr_Prop]: Variable has high volatility around 2010 as seen through the graphical plot. Mean of 0.01228 vs median of 0.012,
#            bounded by [0.01088, 0.01422] for 5%-95% percentiles; no outliers
# [ArrearsToBalance_Aggr_Prop]: Variable has  mean of 0.0006134 vs median of 0.0004895,
#            bounded by [0.0003790, 0.0015335] for 5%-95% percentiles


# - Proportion of curing loans across performing/default spell type
datCredit_smp[, CuringEvents_Aggr_Prop := sum(PerfSpell_Counter==1 & PerfSpell_Num>=2, na.rm=T)/.N, by=list(Date)]
cat( (sum(datCredit_smp[, sum(is.na(CuringEvents_Aggr_Prop)), by=Date][,2])==0) %?% 
       'SAFE: New feature [CuringEvents_Aggr_Prop] has logical values.\n' %:% 
       'WARNING: New feature [CuringEvents_Aggr_Prop] has illogical values \n' )
if (doDescribe)describe(datCredit_smp$CuringEvents_Aggr_Prop); plot(unique(datCredit_smp$CuringEvents_Aggr_Prop), type="b")
### RESULTS: Variable has mean of 0.001373 vs median of 0.0012615,
# bounded by [0.0006567, 0.0026194] for 5%-95% percentiles; no outliers


# - Aggregated age-to-term of portfolio over time, i.e., percentage-based maturity
datCredit_smp[, AgeToTerm_Aggr_Mean := mean(Age_Adj/Term, na.rm=T), by=Date]
cat( (sum(datCredit_smp[, sum(is.na(AgeToTerm_Aggr_Mean)), by=Date][,2])==0) %?% 
       'SAFE: New feature [AgeToTerm_Aggr_Mean] has logical values.\n' %:% 
       'WARNING: New feature [AgeToTerm_Aggr_Mean] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$AgeToTerm_Aggr_Mean); plot(unique(datCredit_smp$AgeToTerm_Aggr_Mean), type="b")
### RESULTS: Variable behaves as expected, i.e., increases as the loan portfolio matures. Has mean 0.3621 and median 0.3878
# bounded by [0.2568, 0.4006] for 5%-95% percentiles; no outliers


# - Aggregate maturity of performance spell ages over time
datCredit_smp[, PerfSpell_Maturity_Aggr_Mean := mean(PerfSpell_Age, na.rm=T), by=Date]
cat( (sum(datCredit_smp[, sum(is.na(PerfSpell_Maturity_Aggr_Mean)), by=Date][,2])==0) %?% 
       'SAFE: New feature [PerfSpell_Maturity_Aggr_Mean] has logical values.\n' %:% 
       'WARNING: New feature [Perf_SpellMaturity_Aggr_Mean] has illogical values \n' )
if (doDescribe) describe(datCredit_smp$PerfSpell_Maturity_Aggr_Mean); plot(unique(datCredit_smp$PerfSpell_Maturity_Aggr_Mean), type="b")
### RESULTS: Mean performance spell age seem to decrease over time. Has mean 135.1 and median 140.68;
# bounded by [93.54, 152.44] for 5%-95% percentiles; no outliers


# - Median-aggregated interest rate margin
# NOTE: The median is preferred over the mean since it resulted in a superior model, as investigated in the experimental script 3c(v)
# Creating an aggregated dataset
dat_IRM_Aggr <- datCredit_smp[, list(InterestRate_Margin_Aggr_Med = median(InterestRate_Margin_imputed_mean, na.rm=T)), by=list(Date)]
# Checking the time series of this variable
plot(dat_IRM_Aggr$InterestRate_Margin_Aggr_Med, type="b")
# Applying various lags
lags <- c(1,2,3,9) # Lags as found to be significant within the experimental script
# Create dataset for conducting sanity checks
dat_IRM_Aggr_Check1 <- data.table(Variable = NULL, Check = NULL)
ColNames <- colnames(dat_IRM_Aggr)[-1] # Names of the columns
for (i in seq_along(lags)){ # Looping over the specified lags and applying each to each of the specified columns
  for (j in seq_along(ColNames)){
    dat_IRM_Aggr[, (paste0(ColNames[j],"_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"),get(ColNames[j]))] # Impute NA's with non-lagged version of variable
  }
}
# [SANITY CHECK] Check whether the lags were created correctly
cat((anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_1]) | anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_2]) | anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_3]) | anyNA(dat_IRM_Aggr[,InterestRate_Margin_Aggr_Med_9])) %?%
      "WARNING: Missingness detected, [InterestRate_Margin_Aggr_Med_1], [InterestRate_Margin_Aggr_Med_2] and/or [InterestRate_Margin_Aggr_Med_3] compromised.\n" %:%
      "SAFE: No missingness, [InterestRate_Margin_Aggr_Med_1], [InterestRate_Margin_Aggr_Med_2] and [InterestRate_Margin_Aggr_Med_3] created successfully.\n")
### RESULTS: Safe, no missingness, hence continue with merge
# Merging the credit dataset with the aggregated dataset
datCredit_smp <- merge(datCredit_smp, dat_IRM_Aggr, by="Date", all.x=T)
# Validate merging success )by checking for missingness (should be zero)
list_merge_variables <- list(colnames(dat_IRM_Aggr))
results_missingness <- list()
for (i in 1:length(list_merge_variables)){
  output <- sum(is.na(datCredit_smp$list_merge_variables[i]))
  results_missingness[[i]] <- output
}
cat( (length(which(results_missingness > 0)) == 0) %?% "SAFE: No missingness, fusion with aggregated data is successful.\n" %:%
       "WARNING: Missingness in certain aggregated fields detected, fusion compromised.\n")
if (doDescribe) describe(datCredit_smp$InterestRate_Margin_Aggr_Med); plot(datCredit_smp[!duplicated(Date),InterestRate_Margin_Aggr_Med], type="b") # Only saving the base variable's descriptive statistics
### RESULTS: Variable follows a logical trend over time. Has mean -0.008077 and median -0.0085;
# bounded by [-0.012, -0.0040] for 5%-95% percentiles; no outliers


# - Save final snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1e"), datCredit_smp)

# Clean up
suppressWarnings(rm(dat_IRM_Aggr, dat_IRM_Aggr_Check1, list_merge_variables, results_missingness, output, lags, ColNames,varSLC_Info_Cat, 
                    varSLC_Info_Num, varCredit_Info_Cat, varCredit_Info_Num, check.fuse1, check.fuse3, check.fuse4, lookup_IDs,
                    Covariate_Info, lookup, lookup2, dat_g0_Delinq_Aggr, dat_DefaultRate, dat_Aggr)); gc()



# --- 9. Macroeconomic feature engineering

# - Confirm that required data objects are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp1e"), tempPath)
if (!exists('datMV')) unpack.ffdf(paste0(genPath,"datMV"), tempPath)

# - Lags of all MVs
# Specifying the lags (monthly) that should be applied
lags <- c(1,2,3,6,9,12)
# Creating a dataset with which to check if the lags are applied correctly to the macroeconomic variables
datMV_Check1 <- data.table(Variable = NULL, Check = NULL)
# Getting the column names with which to apply the lags
ColNames <- colnames(datMV)[-1]
# Looping over the specified lags and applying each to each of the specified columns
for (i in seq_along(lags)){
  for (j in seq_along(ColNames)){
    datMV[, (paste0(ColNames[j],"_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"), get(ColNames[j]))]
  }
}
# [SANITY CHECK] Check datMV for any missingness
cat( anyNA(datMV) %?% "WARNING: Missingness detected in the lagged macroeconomic variables.\n" %:%
       "SAFE: Lags applied successfully to the macroeconomic variables.\n")
### Results: Lagged variables created successfully, no missingness present


# - Rolling volatility/ standard deviations
# Specifying the rolling standard deviation windows (monthly) that should be applied
SD_windows <- c(4,5,6,9,12)
# Creating a dataset with which to check if the standard deviations are applied correctly to the macroeconomic variables
datMV_Check2 <- data.table(Variable = NULL,
                           Check = NULL)
# Getting the column names with which to apply the rolling standard deviations
ColNames <- colnames(datMV)[2:7]
# Looping over the specified lags and applying each to each of the specified columns
for (i in seq_along(SD_windows)){
  for (j in seq_along(ColNames)){
    datMV[, (paste0(ColNames[j],"_SD_",SD_windows[i])) := frollapply(get(ColNames[j]), n=SD_windows[i], FUN=sd, align="right")]
  }
}
# [SANITY CHECK] Check whether the lags were created correctly
cat( anyNA(datMV[Date>="2007-01-31"]) %?% "WARNING: Excessive missingness detected in the calculated SD macroeconomic variables.\n" %:%
       "SAFE: SD macroeconomic variables calculated and created successfully.\n")
### RESULTS: Variable created successfully. Note that there is missingness for some of the SD variables during 1980 because this is the start
# of the macroeconomic dataset and we can't reasonably calculate the SD for these dates since we don't have data from 1979 in the set.
# However, since we are only interested in the data from 2007 onwards, we do not have to worry about this issue and don't impute, since
# these dates will be discarded in any case, implying we would waste computing power if we did decide to impute.


# - Merging the macroeconomic information to the subsampled dataset
datCredit_smp <- merge(datCredit_smp, subset(datMV, select=colnames(datMV)[!(colnames(datMV) %in% ColNames)]), by = "Date", all.x = T)
# - Validate merging success )by checking for missingness (should be zero)
list_merge_variables <- list(colnames(datMV))
results_missingness <- list()
for (i in 1:length(list_merge_variables)){
  output <- sum(is.na(datCredit_smp$list_merge_variables[i]))
  results_missingness[[i]] <- output
}
cat( (length(which(results_missingness > 0)) == 0) %?% "SAFE: No missingness, fusion with macroeconomic data is successful.\n" %:%
       "WARNING: Missingness in certain macroecnomic fields detected, fusion compromised.\n")
### RESULTS: No missingness observed, continue with packing away the data

# - Cleanup
rm(datMV, list_merge_variables, results_missingness, datMV_Check1, datMV_Check2); gc()



# --- 10. Model-based feature engineering

# - Create binned version of TimeInPerfSpell for discrete-time hazard models
timeBinning <- function(x) {
  case_when(
    0 < x & x <= 3 ~ "01.[1,3]", 3 < x & x <= 6 ~ "02.(3,6]",
    6 < x & x <= 9 ~ "03.(6,9]", 9 < x & x <= 12 ~ "04.(9,12]",
    12 < x & x <= 18 ~ "05.(12,18]", 18 < x & x <= 24 ~ "06.(18,24]",
    24 < x & x <= 30 ~ "07.(24,30]", 30 < x & x <= 36 ~ "08.(30,36]",
    36 < x & x <= 48 ~ "09.(36,48]", 48 < x & x <= 60 ~ "10.(48,60]",
    60 < x & x <= 72 ~ "11.(60,72]", 72 < x & x <= 84 ~ "12.(72,84]",
    84 < x & x <= 96 ~ "13.(84,96]", 84 < x & x <= 96 ~ "14.(84,96]",
    96 < x & x <= 108 ~ "15.(96,108]", 108 < x & x <= 120 ~ "16.(108,120]",
    120 < x & x <= 144 ~ "17.(120,144]", 144 < x & x <= 168 ~ "18.(144,168]",
    168 < x & x <= 192 ~ "19.(168,192]",TRUE ~ "20.193+"
  )
}
datCredit_smp[, Time_Binned := timeBinning(TimeInPerfSpell)]
table(datCredit_smp$Time_Binned) %>% prop.table()
### RESULTS: Between 3% and 7% of observations in each bin; deemed appropriate, particulalarly in
# capturing the earlier time periods (shorter interval bins).

# Lag g0-delinq with appropriate period for discrete-time hazard model
datCredit_smp[, g0_Delinq_Lag_1 := shift(g0_Delinq,fill=0),by=LoanID]

# - Create start point variable
datCredit_smp[, Start := TimeInPerfSpell - 1]

# -- Save fused- and enriched subsampled dataset for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp2"), datCredit_smp)




# ------ 5. Apply a basic cross-validation clustered resampling scheme with possible n-way stratification

# - Confirm that required data objects are loaded into memory
if (!exists('datCredit_smp')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_smp2"), tempPath)
if (!exists('datKeys_sampled')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST_keys"), tempPath)

# - Implement the clustered (possibly stratified) resampling scheme by first randomly selecting loan IDs 
# using the given sampling fraction
set.seed(1, kind="Mersenne-Twister")
if (all(!is.na(stratifiers))){ # enforce Stratifiers
  dat_train_keys <- datKeys_sampled %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_frac) %>% as.data.table() 
} else { # No stratifiers
  dat_train_keys <- datKeys_sampled %>% slice_sample(prop=smp_frac) %>% as.data.table()
}

# - Extract the entire loan histories into the training set for those randomly select subject IDs
# Select only the first performing spell (given the model definition), while 
# the validation set deliberately includes multiple spells to test certain modelling assumptions
if (timeDef_TFD) {
  datCredit_train_TFD <- copy(datCredit_smp[get(clusVar) %in% dat_train_keys[, get(clusVar)],]) %>% 
    subset(PerfSpell_Num == 1)
} else {
  datCredit_train_PWPST <- copy(datCredit_smp[get(clusVar) %in% dat_train_keys[, get(clusVar)],])
}

# - Extract the entire loan histories into the validation set for those remaining subjects
datCredit_valid_PWPST <- copy(datCredit_smp[!(get(clusVar) %in% dat_train_keys[, get(clusVar)]),]);gc()

# - [SANITY CHECKS]
if (timeDef_TFD) {
  # Can subsample be reconstituted?
  check.1 <- datCredit_smp[,.N] == datCredit_train_TFD[,.N] + datCredit_valid_TFD[,.N] + datCredit_smp[get(clusVar_Spell) %in% vSpellKeys_MultiSpell,.N] # Should be TRUE
  # Does training set contain only first-time spells?
  check.2 <- datCredit_train_TFD[get(spellNum) == 1,.N] == datCredit_train_TFD[,.N] # Should be TRUE
  # Does validation spell contain spell numbers other than 1?
  (check.3 <- datCredit_valid_TFD[get(spellNum) != 1,.N] > 0) # Should be TRUE
} else {
  # Can subsample be reconstituted?
  check.1 <- datCredit_smp[,.N] == datCredit_train_PWPST[,.N] + datCredit_valid_PWPST[,.N]
  # Does training set contain only first-time spells?
  check.2 <- T # Irrelevant for this time definition, so assign default
  # Does validation spell contain spell numbers other than 1?
  check.3 <- T # Irrelevant for this time definition, so assign default
}
cat((check.1 %?% "SAFE: Training and validation datasets succcessfully reconstitute the subsampled dataset. \n" %:% 
       'WARNING: Training and validation datasets do not reconstitue the subsampled dataset. \n' ))
cat((check.2 %?% paste0("SAFE: Spells in the training dataset are selected as desired; First-spells only? ", timeDef_TFD, ".\n" ) %:%
       paste0("WARNING: Spells in the training dataset are not selected as desired; First-spells only? ", timeDef_TFD, ".\n" )))
cat((check.3 %?% paste0("SAFE: Spells in the validation dataset are selected as desired; Multi-spells? ", check.3, ".\n" ) %:%
       paste0("WARNING: Spells in the validation dataset are not selected as desired; Multi-spells? ", check.3, ".\n" )))

# - Clean up
suppressWarnings(rm(smp_perc, datKeys, datKeys_sampled, dat_train_keys, check.2, check.3, datCredit_smp_old_n, 
                    datCredit_smp_prior, dat_keys_exc, class_type, excCond, excCond2, datCredit_smp))


# --- 5.2 Saving the cross-validation scheme
# - Training dataset
pack.ffdf(paste0(genPath,"creditdata_train_PWPST"), datCredit_train_PWPST)

# - Validation dataset
pack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), datCredit_valid_PWPST)


# --- 5.3 Clean up
suppressWarnings(rm(dat_keys_smp_perf, dat_keys_smp_perf,  dat_train_keys_perf, dat_train_keys_def, datCredit_train_perf, datCredit_train_def,  datCredit_valid_perf, datCredit_valid_def,
                    check.4_a, check.4_b, check.4_c, check.5_a, check.5_b, datCredit_smp, datStrata_smp_min, datCredit_train_PWPST, datCredit_valid_PWPST));gc()

