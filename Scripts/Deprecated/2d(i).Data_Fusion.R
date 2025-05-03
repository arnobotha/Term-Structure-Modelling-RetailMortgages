# ===================================== DATA FUSION =====================================
# Fuse the prepared and previously enriched credit dataset with the base macroeconomic
# variables from the macroeconomic data and apply feature engineering to create some new
# variables.
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Esmerelda Oberholzer, Roelinde Bester, Marcel Muller, Roland Breedt
# ---------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2b.Data_Preparation_Credit.R
#   - 2c.Data_Enrich.R

# -- Inputs:
#   - datCredit_real | prepared credit data from script 2c
#   - various parameters set in the setup script 0
#
# -- Outputs:
#   - datCredit_finala | Subsetted dataset, with exclusions applied and some feature engineering
# ---------------------------------------------------------------------------------------
# NOTE: This script predominantly comes from another project (Kasmeer).
# =======================================================================================




# ------- 1. Remove some unnecessary variables to increase memory

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final3"), tempPath)

# - Subset the dataset to exclude unnecessary variables
datCredit_real <- subset(datCredit_real, select = -c(Age, AccountStatus, DelinqState_g0,
                                                     ZeroBal_Start, NCA_CODE, STAT_CDE,
                                                     CLS_STAMP, BOND_IND, PerfSpell_Event, DefSpell_Event))




# ------- 2. Apply exclusions on the credit dataset to increase available memory

# - Check the impact of the exclusions from script 2d | RECORD-LEVEL
(exclusions_credit <- datCredit_real[ExclusionID != 0, .N] / datCredit_real[, .N] * 100)
# Exclusion's impact: 5.881418%

# - Check the combined impact for possible overlaps | RECORD-LEVEL
(exclusions_all <- datCredit_real[ExclusionID != 0 | is.na(PerfSpell_Num), .N] / datCredit_real[, .N] * 100)
# Total exclusions' impact:  10.43866%

# - Now apply the exclusions (but keep the default exposures)
datCredit_real <- subset(datCredit_real, ExclusionID == 0); gc()

# - Checks
sum(datCredit_real$ExclusionID > 0) == 0 # check - success

# - Remove fields that will not likely be used in the eventual analysis/ modelling of default risk, purely to save memory
datCredit_real <- subset(datCredit_real, select = -c(ExclusionID))




# ------ 3. Feature Engineering that requires entire/ complete loan histories
# --- Default indicators
# - Last-status approach
maxDate <- max(datCredit_real[,Date], na.rm=T) - years(1) # Dates larger than maxDate do not have 12-month default because of the end of the sampling window
datCredit_real[, DefaultStatus1_lead_12 := ifelse(Date<=maxDate,imputeLastKnown(shift(DefaultStatus1, n=12, type="lead")), NA), by=list(LoanID)]
datCredit_real$DefaultStatus1_lead_12 %>% table() %>% prop.table()
### RESUKTS: 93.28 % of observations have not defaulted in exactly 12 months since reporting date and 6.72% did.
# Relocate variable next to current default status variable
datCredit_real <- datCredit_real %>% relocate(DefaultStatus1_lead_12, .after=DefaultStatus1)

# - Worst-ever approach
datCredit_real[, DefaultStatus1_lead_12_max := ifelse(Date<=maxDate,imputeLastKnown(frollapply(x=DefaultStatus1, n=13, align="left", FUN=max)), NA), by=list(LoanID)]
datCredit_real$DefaultStatus1_lead_12_max %>% table() %>% prop.table()
### RESULTS: 91.84% of observations have not defaulted in the next 12 months from reporting date, whilst 8.16% of accounts had.
# Relocate variable next to current default status variable
datCredit_real <- datCredit_real %>% relocate(DefaultStatus1_lead_12_max, .after=DefaultStatus1_lead_12)


# --- Delinquency-themed variables on a loan-level
# - Embed previous defaults into a new Boolean-valued input variable
datCredit_real[, PrevDefaults := ifelse(all(is.na(PerfSpell_Num)), F, max(PerfSpell_Num,na.rm = T) > 1), by=list(LoanID)]
cat( (datCredit_real[is.na(PrevDefaults), .N] == 0) %?% "SAFE: No missingness, [PrevDefaults] created successfully.\n" %:%
       "WARNING: Missingness detected, [PrevDefaults] compromised.\n")
describe(datCredit_real$PrevDefaults); describe(datCredit_real[Counter==1, PrevDefaults])
### RESULTS: 11.6% of records had previous defaults, which is 6.9% of accounts

# - Spell-level indicator for when a shift occurs in the state of g0_Delinq (target event)
# NOTE: This is an intermediary field used in the creation of subsequent fields
datCredit_real[, g0_Delinq_Shift := ifelse(lag(g0_Delinq, n=1)==g0_Delinq,0,1), by=list(LoanID)]
datCredit_real[is.na(g0_Delinq_Shift), g0_Delinq_Shift := 0] # All first observations have g0_Delinq_Shift = NA; set these values to zero.
cat( (datCredit_real[is.na(g0_Delinq_Shift), .N] == 0) %?% "SAFE: No missingness, [g0_Delinq_Shift] created successfully.\n" %:%
       "WARNING: Missingness detected, [g0_Delinq_Shift] compromised.\n")
datCredit_real$g0_Delinq_Shift %>% table() %>% prop.table()
### RESULT: 96.05% of the records had no change in their delinquency level from their associated previous record.

# - Delinquency state number, where each change in g_0 denotes such a "state" that may span several periods
datCredit_real[, g0_Delinq_Num := cumsum(g0_Delinq_Shift) + 1, by=list(LoanID)] # Assign state numbers over the entire loan history (add one to ensure that there are no delinquency spell numbers equal to zero)
cat( (datCredit_real[is.na(g0_Delinq_Num), .N] == 0) %?% "SAFE: No missingness, [g0_Delinq_Num] created successfully.\n" %:%
       "WARNING: Missingness detected, [g0_Delinq_Num] compromised.\n")
describe(datCredit_real$g0_Delinq_Num)
### RESULT: Mean state number of 3.3 across all rows; median: 1; max of 100. 
# This high max suggests outlier-accounts with rapid and frequent changes in g0

# - Account-level standard deviation of the delinquency state
datCredit_real[, g0_Delinq_SD := sd(g0_Delinq, na.rm=T), by=list(LoanID)]
datCredit_real[is.na(g0_Delinq_SD), g0_Delinq_SD := 0] # Some missing values exist at loan accounts originating at the end of the sampling period | Assign zero values to these
cat( (datCredit_real[is.na(g0_Delinq_SD), .N] == 0) %?% "SAFE: No missingness, [g0_Delinq_SD] created successfully.\n" %:%
       "WARNING: Missingness detected, [g0_Delinq_SD] compromised.\n")
describe(datCredit_real[, list(g0_Delinq_SD=mean(g0_Delinq_SD, na.rm=T)), by=list(LoanID)]$g0_Delinq_SD)
### RESULT: mean account-level SD in delinquency states of 0.21; median: 0, but 95%-percentile of 1.2
# This suggests that most accounts do not vary significantly in their delinquency states over loan life, which is sensible

# - 4-,5-,6-,9- and 12 month rolling state standard deviation
# NOTE: Usefulness of each time window length will yet be determined during prototyping/modelling
SD_LoanLevel<-datCredit_real[,list(SD_Loans=sd(g0_Delinq,na.rm=TRUE)),by=list(LoanID)] # Create standard deviation variable for each loan account
SD_Overall<-mean(SD_LoanLevel[,SD_Loans],na.rm=TRUE) # Obtain mean SD over loan accounts for imputation of NA's at the beginning of each loan's history

# frollapply function lags the variable across some fixed window
datCredit_real[, g0_Delinq_SD_12 := frollapply(g0_Delinq, n=12, FUN=sd, align="right",fill=SD_Overall), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_9 := frollapply(g0_Delinq, n=9, FUN=sd, align="right",fill=SD_Overall), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_6 := frollapply(g0_Delinq, n=6, FUN=sd, align="right",fill=SD_Overall), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_5 := frollapply(g0_Delinq, n=5, FUN=sd, align="right",fill=SD_Overall), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_4 := frollapply(g0_Delinq, n=4, FUN=sd, align="right",fill=SD_Overall), by=list(LoanID)]

# [SANITY CHECK] Check for missingness in engineered variables
cat((anyNA(datCredit_real[,g0_Delinq_SD_12]) | anyNA(datCredit_real[,g0_Delinq_SD_9]) | anyNA(datCredit_real[,g0_Delinq_SD_6])
     | anyNA(datCredit_real[,g0_Delinq_SD_5]) | anyNA(datCredit_real[,g0_Delinq_SD_4])) %?% "WARNING: Excessive missingness detected, [g0_Delinq_SD_4], [g0_Delinq_SD_5], [g0_Delinq_SD_6], [g0_Delinq_SD_9], and/or [g0_Delinq_SD_12] compromised.\n" %:%
      "SAFE: No missingness, [g0_Delinq_SD_4], [g0_Delinq_SD_5], [g0_Delinq_SD_6], [g0_Delinq_SD_9], and [g0_Delinq_SD_12] created successfully.\n")

#Clean-up
rm(SD_LoanLevel,SD_Overall)

# - Proportion of new loans vs existing portfolio over time
# NOTE: we therefore measure credit demand within market, underlying market conditions, and the implicit effect of bank policies)
# Creating an aggregated dataset
dat_NewLoans_Aggr <- datCredit_real[, list(NewLoans_Aggr_Prop = sum(Age_Adj==1, na.rm=T)/.N), by=list(Date)]
# Applying various lags
lags <- c(1,3,4,5) # Lags
ColNames <- colnames(dat_NewLoans_Aggr)[-1] # Names of the columns
for (i in seq_along(lags)){ # Looping over the specified lags and applying each to each of the specified columns
  for (j in seq_along(ColNames)){
    dat_NewLoans_Aggr[, (paste0(ColNames[j],"_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"),get(ColNames[j]))] # Impute NA's with the non lagged value
  }
}
# [SANITY CHECK] Check whether the lags were created correctly
cat((anyNA(dat_NewLoans_Aggr[,NewLoans_Aggr_Prop_1]) | anyNA(dat_NewLoans_Aggr[,NewLoans_Aggr_Prop_3]) | anyNA(dat_NewLoans_Aggr[,NewLoans_Aggr_Prop_4])
     | anyNA(dat_NewLoans_Aggr[,NewLoans_Aggr_Prop_5])) %?% "WARNING: Missingness detected, [NewLoans_Aggr_Prop_1], [NewLoans_Aggr_Prop_3], [NewLoans_Aggr_Prop_4], and/or [NewLoans_Aggr_Prop_5] compromised.\n" %:%
      "SAFE: No missingness, [NewLoans_Aggr_Prop_1], [NewLoans_Aggr_Prop_3], [NewLoans_Aggr_Prop_4], and [NewLoans_Aggr_Prop_5] created successfully.\n")
### RESULTS: Variables successfully created without any missingness

# Merging the credit dataset with the aggregated dataset
datCredit_real <- merge(datCredit_real, dat_NewLoans_Aggr, by="Date", all.x=T)
# Validate merging success by checking for missingness (should be zero)
list_merge_variables <- list(colnames(dat_NewLoans_Aggr))
results_missingness <- list()
for (i in 1:length(list_merge_variables)){
  output <- sum(is.na(datCredit_real$list_merge_variables[i]))
  results_missingness[[i]] <- output
}
cat( (length(which(results_missingness > 0)) == 0) %?% "SAFE: No missingness, fusion with aggregated data is successful.\n" %:%
       "WARNING: Missingness in certain aggregated fields detected, fusion compromised.\n")
describe(datCredit_real$NewLoans_Aggr_Prop); plot(unique(datCredit_real$NewLoans_Aggr_Prop), type="b")
### RESULTS: Variable has mean of 0.008 vs median of 0.007,
# bounded by [0.003, 0.014] for 5%-95% percentiles; no major outliers


# - Time in delinquency state
# NOTE: This variable is conceptually different to [TimeInPerfSpell].
# A performance spell starts when a loan is not in default and ends when it is in default.
# A delinquency spell starts when a loan "shifts" to a new delinquency level and ends the immediate period preceding the next shift to a different delinquency level.
datCredit_real[, TimeInDelinqState := 1:.N, by=list(LoanID, g0_Delinq_Num)]
cat( (datCredit_real[is.na(TimeInDelinqState), .N] == 0) %?% "SAFE: No missingness detected, [TimeInDelinqState] created successfully.\n" %:%
       "WARNING: Missingness detected, [TimeInDelinqState] compromised.\n")


# --- Delinquency-themed variables on a performance spell-level
# - Delinquency state number, where each change in g_0 denotes such a "state" that may span several periods during a performance spell
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_g0_Delinq_Num := cumsum(g0_Delinq_Shift) + 1, by=list(PerfSpell_Key)] # Assign state numbers over each performance spell
# [SANITY CHECK] Check new feature for illogical values
cat( ( datCredit_real[is.na(PerfSpell_g0_Delinq_Num),.N]==datCredit_real[is.na(PerfSpell_Key),.N]) %?% 
       'SAFE: New feature [PerfSpell_g0_Delinq_Num] has logical values.\n' %:% 
       'WARNING: New feature [PerfSpell_g0_Delinq_Num] has illogical values \n' )

# - State standard deviation on the performance spell level
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_g0_Delinq_SD := sd(g0_Delinq), by=list(PerfSpell_Key)]
datCredit_real[!is.na(PerfSpell_Key) & is.na(PerfSpell_g0_Delinq_SD), PerfSpell_g0_Delinq_SD := 0] # Assigning an standard deviation of zero to those performance spells that have an single observation
# [SANITY CHECK] Check new feature for illogical values
cat( ( datCredit_real[is.na(PerfSpell_g0_Delinq_SD),.N]==datCredit_real[is.na(PerfSpell_Key),.N]) %?% 
       'SAFE: New feature [PerfSpell_g0_Delinq_SD] has logical values.\n' %:% 
       'WARNING: New feature [PerfSpell_g0_Delinq_SD] has illogical values \n')


# 4. ------ Saving the final dataset and doing some housekeeping

# - remove intermediary fields, as a memory enhancement
datCredit_real[, g0_Delinq_Shift := NULL]

# --- Save snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final4"), datCredit_real)

# --- Housekeeping
suppressWarnings(rm(exclusions_all, exclusions_credit)); gc()






