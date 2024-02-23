# ===================================== DATA FUSION =====================================
# Fuse the prepared and previously enriched credit dataset with the base macroeconomic
# variables from the macroeconomic data and apply feature engineering to create some new
# variables.
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Esmerelda Oberholzer, Roelinde Bester, Marcel Muller
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
datCredit_real[, DefaultStatus1_lead_12 := imputeLastKnown(shift(DefaultStatus1, n=12, type="lead")), by=list(LoanID)]
datCredit_real$DefaultStatus1_lead_12 %>% table() %>% prop.table()
### RESUKTS: 93.35 % of observations have not defaulted in exactly 12 months since reporting date and 6.75% did.
# Relocate variable next to current default status variable
datCredit_real <- datCredit_real %>% relocate(DefaultStatus1_lead_12, .after=DefaultStatus1)

# - Worst-ever approach
datCredit_real[, DefaultStatus1_lead_12_max := imputeLastKnown(frollapply(x=DefaultStatus1, n=13, align="left", FUN=max)), by=list(LoanID)]
datCredit_real$DefaultStatus1_lead_12_max %>% table() %>% prop.table()
### RESULTS: 91.93% of observations have not defaulted in the next 12 months from reporting date, whilst 7.96% of accounts had.
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
datCredit_real[, g0_Delinq_SD_12 := frollapply(g0_Delinq, n=12, FUN=sd, align="right"), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_9 := frollapply(g0_Delinq, n=9, FUN=sd, align="right"), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_6 := frollapply(g0_Delinq, n=6, FUN=sd, align="right"), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_5 := frollapply(g0_Delinq, n=5, FUN=sd, align="right"), by=list(LoanID)]
datCredit_real[, g0_Delinq_SD_4 := frollapply(g0_Delinq, n=4, FUN=sd, align="right"), by=list(LoanID)]
cat( ((datCredit_real[is.na(g0_Delinq_SD_4), .N] == datCredit_real[Counter<4,.N]) &
        (datCredit_real[is.na(g0_Delinq_SD_5), .N] == datCredit_real[Counter<5,.N]) &
        (datCredit_real[is.na(g0_Delinq_SD_6), .N] == datCredit_real[Counter<6,.N]) &
        (datCredit_real[is.na(g0_Delinq_SD_9), .N] == datCredit_real[Counter<9,.N]) &
        (datCredit_real[is.na(g0_Delinq_SD_12), .N] == datCredit_real[Counter<12,.N])) %?% "SAFE: No excessive missingness, [g0_Delinq_SD_4], [g0_Delinq_SD_5], [g0_Delinq_SD_6], [g0_Delinq_SD_9], and [g0_Delinq_SD_12] created successfully.\n" %:%
       "WARNING: Excessive missingness detected, [g0_Delinq_SD_4], [g0_Delinq_SD_5], [g0_Delinq_SD_6], [g0_Delinq_SD_9], and/or [g0_Delinq_SD_12] compromised.\n")

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




# 4. ------ Saving the final dataset and doing some housekeeping
# --- Save snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final4"), datCredit_real)

# --- Housekeeping
suppressWarnings(rm(exclusions_all, exclusions_credit)); gc()






