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
datCredit_real$DefaultStatus1_lead_12 %>% table() %>% prop.table() # 93.35 % of observations have not defaulted in exactly 12 months since reporting date and 6.75% did.
# Relocate variable next to current default status variable
datCredit_real <- datCredit_real %>% relocate(DefaultStatus1_lead_12, .after=DefaultStatus1)

# - Worst-ever approach
datCredit_real[, DefaultStatus1_lead_12_max := imputeLastKnown(frollapply(x=DefaultStatus1, n=13, align="left", FUN=max)), by=list(LoanID)]
datCredit_real$DefaultStatus1_lead_12_max %>% table() %>% prop.table() # 91.93% of observations have not defaulted in the next 12 months from reporting date, whilst 7.96% of accounts had.
# Relocate variable next to current default status variable
datCredit_real <- datCredit_real %>% relocate(DefaultStatus1_lead_12_max, .after=DefaultStatus1_lead_12)


# --- Creating delinquency-spell level input variables
# - Indicator for when a shift in the state of g0_Delinq occurs (target event) --- Loan level
datCredit_real[, g0_Delinq_Shift := ifelse(lag(g0_Delinq, n=1)==g0_Delinq,0,1), by=list(LoanID)]
datCredit_real[is.na(g0_Delinq_Shift), g0_Delinq_Shift := 0] # All first observations have g0_Delinq_Shift = NA; set these values to zero.
datCredit_real$g0_Delinq_Shift %>% table() %>% prop.table() # [g0_Delinq_Shift] has 2 levels and no missing values. 96.07% of observations have a shift in delinquency and  3.95% of observations have no shifts in delinquency.

# - State number --- Loan level
datCredit_real[, g0_Delinq_Num := cumsum(g0_Delinq_Shift), by=list(LoanID)] # Assign state numbers over the entire loan history
datCredit_real[, g0_Delinq_Num := g0_Delinq_Num + 1] # Small adjustment to ensure that spell numbers can't be zero.
hist(datCredit_real$g0_Delinq_Num, breaks='FD')

# - State volatility --- Loan level
datCredit_real[, g0_Delinq_SD := sd(g0_Delinq), by=LoanID]
hist(datCredit_real$g0_Delinq_SD, breaks='FD')

# - [SANITY CHECK]
lookup <- datCredit_real[LoanID=="3000003205066", list(LoanID, Date, PerfSpell_Key, DefSpell_Key,
                                                      g0_Delinq, g0_Delinq_Shift, PerfSpell_g0_Delinq_Shift, DefSpell_g0_Delinq_Shift,
                                                      g0_Delinq_Num, DefSpell_g0_Delinq_Num, PerfSpell_g0_Delinq_Num,
                                                      g0_Delinq_SD, PerfSpell_g0_Delinq_SD, DefSpell_g0_Delinq_SD)]
lookup2 <- datCredit_real[LoanID=="3000005788425", list(LoanID, Date, PerfSpell_Key, DefSpell_Key,
                                                       g0_Delinq, g0_Delinq_Shift, PerfSpell_g0_Delinq_Shift, DefSpell_g0_Delinq_Shift,
                                                       g0_Delinq_Num, DefSpell_g0_Delinq_Num, PerfSpell_g0_Delinq_Num,
                                                       g0_Delinq_SD, PerfSpell_g0_Delinq_SD, DefSpell_g0_Delinq_SD)]
### RESULTS:~ All variables created correctly, no problems detected.
###           Safe to proceed.

# - Clean up
rm(lookup, lookup2)

# - Save snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final4"), datCredit_real)

# - Housekeeping
suppressWarnings(rm(exclusions_all, exclusions_credit)); gc()
