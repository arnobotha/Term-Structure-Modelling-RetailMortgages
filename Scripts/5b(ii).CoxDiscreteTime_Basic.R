# ======================================= INPUT SPACE: Cox PO (PWP) ADVANCED============================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the TPWPST model
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2f.Data_Fusion1.R
#   - 3b.Data_Fusion2_PWP_ST.R

# -- Inputs:
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b
#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------





# ------ 1. Final model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]
# remove previous objects from memory
rm(datCredit_train_PWPST, datCredit_valid_PWPST); gc()

# - Fit an "empty" model as a performance gain, used within some diagnostic functions
modLR_base <- glm(PerfSpell_Event ~ 1, data=datCredit_train, family="binomial")

# - Final variables
# Selection based on expert judgement alone, where the initial list being based on thematic selection
vars_basic <- c("-1", "Time_Binned", "log(TimeInPerfSpell):PerfSpell_Num_binned",
                "Arrears", "InterestRate_Nom", "M_Inflation_Growth_6")
modLR_basic <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars_basic, collapse = " + "))),
              data=datCredit_train, family="binomial")
summary(modLR_basic);
evalLR(modLR_basic, modLR_base, datCredit_train, targetFld="PerfSpell_Event", predClass=1)
### RESULTS: AIC:  128,556;  McFadden R^2:  27.78%; AUC:  93.80%.

# - Test goodness-of-fit using AIC-measure over single-factor models
aicTable_CoxDisc_basic <- aicTable(datCredit_train, vars_basic, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Top variables: Arrears, Time_Binned, InterestRate_Nom, M_Inflation_Growth_6

# Test accuracy using c-statistic over single-factor models
concTable_CoxDisc_basic <- concTable(datCredit_train, datCredit_valid, vars_basic, TimeDef=c("Cox_Discrete","PerfSpell_Event"), genPath=genObjPath, modelType="Cox_Discrete")
# Top variables: Arrears, Time_Binned, InterestRate_Nom, M_Inflation_Growth_6 

# - Combine results into a single object
Table_CoxDisc_basic <- concTable_CoxDisc_basic[,1:2] %>% left_join(aicTable_CoxDisc_basic, by ="Variable")

# Save objects
pack.ffdf(paste0(genObjPath,"CoxDisc_basic_fits"), Table_CoxDisc_basic)

