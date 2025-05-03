# ================================== DATA DIAGNOSTICS ===================================
# Conduct a set of diagnostics on stress-testing the structure of a given dataset
# according to set of assumptions
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha

# DESCRIPTION: 
# Diagnostics are recorded tests (statistical and logical) whose results are reported
# in the eventual PDF-report
# ---------------------------------------------------------------------------------------
# -- Script dependencies:


# -- Inputs:
#   - macro_data_hist | prepared macroeconomic data from script 2b(i)
#   - datCredit_real | prepared credit data from script 2b(ii)
#   - various parameters set in the setup script 0
#
# -- Outputs:
#   - <Diagnostics>
# ---------------------------------------------------------------------------------------
# NOTE: This script predominantly comes from another project (Kasmeer).
# =======================================================================================



# ------ 1. Portfolio diagnostics: general distributions and prevalence rates

# - 1a) Starting/stopping time of observing historical portfolio and of forecast-parts, including loan production  | Portfolio diagnostics
# Historical (incomplete) portfolio
test <- subset(datCredit_real, ExclusionID==0 & Partition == "Observed")
diag1a_real <- list(min(test$Date), max(test$Date), 
               interval(min(test$Date), max(test$Date)) %/% years(1),
               interval(min(test$Date), max(test$Date)) %/% months(1),
               test[,.N])
diag1b_real <- list(min(test[, Date_Origination]), max(test[, Date_Origination]),
               interval(min(test[, Date_Origination]), max(test[, Date_Origination])) %/% years(1),
               interval(min(test[, Date_Origination]), max(test[, Date_Origination])) %/% months(1))
# Forecast-part of completed portfolio
test <- subset(datCredit_real, ExclusionID==0 & Partition == "Forecast")
if (test[, .N] > 0) {
  diag1c_real <- list(min(test$Date), max(test$Date), 
                      interval(min(test$Date), max(test$Date)) %/% years(1),
                      interval(min(test$Date), max(test$Date)) %/% months(1),
                      test[,.N])
  diag1d_real <- list(min(test[, Date_Origination]), max(test[, Date_Origination]),
                      interval(min(test[, Date_Origination]), max(test[, Date_Origination])) %/% years(1),
                      interval(min(test[, Date_Origination]), max(test[, Date_Origination])) %/% months(1))  
} else {
  diag1c_real <- NA
  diag1d_real <- NA
}
gc() # Memory optimization


# - 1c) What is the distribution of event types and right-censoring? | Incomplete portfolio | Portfolio diagnostics
test <- subset(datCredit_real, ExclusionID==0 & Counter == Max_Counter & Partition == "Observed")
test[, Event_Type2 := case_when(
  Event_Type == "ACTIVE" ~ "4) Censored",
  Event_Type == "COMPLETED" ~ "3) Paid-up",
  Event_Type == "SETTLE" ~ "1) Early settlements",
  Event_Type == "WOFF" ~ "2) Write-offs")]
(diag2b_real <- table(test$Event_Type2) %>% proportions())
describe(test$Event_Type2)
### RESULT: 2% Completed, 54% Settled, 5% Written-off, 40% Censored


# - 1d) Prevalence of right-censored active cases across the whole dataset as at most recent observation point?
test <- datCredit_real[ExclusionID==0 & Partition == "Observed" & Counter == Max_Counter, 
                         list(Maturity_Desc = ifelse(Event_Type %in% c("SETTLE","WOFF","COMPLETED"), 
                                                     "Settled/Written-off/Paid-up", "Active/right-censored"))]
(diag3_real <- table(test$Maturity_Desc) %>% proportions())
describe(as.factor(test$Maturity_Desc))
(diag3b_real <- c( ifelse(diag3_real[1] > 0.5, "\\textbf{immature} (as expected)", 
                          "\\textbf{\\textcolor{red}{mature}} (unexpected)"),
              ifelse(diag3_real[1] > 0.5, "extensive", "little")
))
### RESULT: 60% completed, the rest right-censored/active cases


# - 1e) How prevalent were settled cases that were in default before settlement? Are all of these allocated as 'cures'?
samp1 <- subset(datCredit_real, ExclusionID==0 & Event_Type == "SETTLE")
defPriorSettled_real <- sum(samp1[Counter==Max_Counter, DefaultStatus1]) / samp1[Counter==Max_Counter, .N] * 100
samp2 <- subset(datCredit_real, ExclusionID==0 & Event_Type == "SETTLE" & Counter==Max_Counter & DefaultStatus1 == 1)
defPriorSettled_cured_real <- samp2[DefSpellResol_Type_Hist == "Cured", .N] / samp2[, .N] * 100
(diag4_real <- c(defPriorSettled_real, defPriorSettled_cured_real))
diag4_real[2] == 100; rm(samp1); gc()
### RESULTS: If TRUE, then PASSED the test. 
# ~9% of the settled loans were in default at the settlement point. 100% of these cases
# are correctly designated as cures. Logic confirmed.


# - 1f) Default incidence over time (historical) | Portfolio diagnostics
outcome_period <- 12; sample_period <- 60 # months
test <- datCredit_real[ExclusionID==0, list(Date, DefaultStatus1, DefaultStatus_lead = 
                           shift(DefaultStatus1, n=outcome_period, type="lead")), 
                  by=list(LoanID)][DefaultStatus1 == 0, 
                                   list(Count=.N, Defaults = sum(DefaultStatus_lead, na.rm=T),
                                        DefaultRate = sum(DefaultStatus_lead, na.rm=T) / .N), 
                                   by=list(Date)][Date>=rollback(maxDate_observed-months(sample_period), roll_to_first = T) & 
                                                    Date<=rollback(maxDate_observed-months(outcome_period), roll_to_first = T),]
ggplot(test, aes(x=Date, y=DefaultRate)) + geom_line() + theme_minimal()
describe(test$DefaultRate); hist(test$DefaultRate, breaks=2*test[, .N]^(1/3))
stat.test <- shapiro.test(test$DefaultRate)
(diag5_real <- c(outcome_period, sample_period, mean(test$DefaultRate, na.rm=T), median(test$DefaultRate, na.rm=T), 
            min(test$DefaultRate, na.rm=T), max(test$DefaultRate, na.rm=T), stat.test$p.value,
            skewness(test$DefaultRate), 1.96*sd(test$DefaultRate)/sqrt(test[,.N])
))
diag5_real[7] > 0.05
### RESULTS: If test is TRUE, then default rate is normally distributed, therefore stable.
# Overall, default rate averages at 2%


# - 1g) Monthly loan production statistics | Portfolio diagnostics
# Incomplete portfolio
test <- subset(datCredit_real, ExclusionID==0 & Counter == 1 & Partition == "Observed" & New_Ind == 1) %>% 
  group_by(Date) %>% summarise(Count=n())
describe(test$Count); hist(test$Count);
stat.test <- shapiro.test(test$Count)
(diag6_real <- c(mean(test$Count, na.rm=T), median(test$Count, na.rm=T),
            min(test$Count, na.rm=T), max(test$Count, na.rm=T), stat.test$p.value,
            skewness(test$Count), sd(test$Count)
))
diag6_real[5] > 0.05
### RESULTS: If test is TRUE, then loan production is normally distributed, therefore stable.





# ------ 2. Default spell diagnostics

# - 2a) What is the sample variance of spell-level loss rates?
# (Loss rate datCredit_real to be calculated as a single time-invariant value discounted to default point)
samp1 <- datCredit_real[ExclusionID==0, list(Var = var(LossRate_Real)),by=list(LoanID, DefSpell_Num)]
(diag.defspell1_real <- samp1[Var>0,.N] / samp1[,.N])
diag.defspell1_real == 0
### RESULT: If TRUE, then PASSED the test. We expect zero cases where default spells had a varying loss rate


# - 2b) Uniqueness of TimeInDefSpell as a filter mechanism?
(diag.defspell3_real <- datCredit_real[ExclusionID==0 & DefSpell_Num >= 1 & TimeInDefSpell==1, list(Age_Adj), 
                                       by=list(LoanID)][,.N] / 
    datCredit_real[ExclusionID==0 & DefSpell_Num >= 1 & DefSpell_LeftTrunc==0, list(Freq = .N), 
                   by=list(LoanID, DefSpell_Num)][,.N] * 100)
diag.defspell3_real == 100
### RESULT: If TRUE, then test is PASSED and structure/grain is confirmed.


# - 2c) Distribution of inferred number of default episodes (expected to be right-skewed) | Complete portfolio
test <- subset(datCredit_real, ExclusionID==0 & DefSpell_Num >= 1 & TimeInDefSpell==1)[,list(Max_Epi = max(DefSpell_Num)), by=list(LoanID)]
describe(test$Max_Epi); hist(test$Max_Epi, breaks=2*length(test$Max_Epi)^(1/3))
(diag.defspell4a_real <- table(test$Max_Epi) %>% proportions())
(diag.defspell4b_real <- c(round(mean(test$Max_Epi, na.rm=T),digits=2), round(median(test$Max_Epi, na.rm=T),digits=2), 
                      max(test$Max_Epi, na.rm=T),
                      round(skewness(test$Max_Epi),digits=4), case_when(skewness(test$Max_Epi) >= 0.1 ~ "right-skewed", 
                                                                        skewness(test$Max_Epi) < 0.1 & skewness(test$Max_Epi) > -0.1 ~ "symmetric",
                                                                        skewness(test$Max_Epi) <= -0.1 ~ "left-skewed") ) ) 
### RESULT: Confirmed, highly right-skewed. Serves as evidence of redefault, including the degree thereof


# - 2d) Distribution of inferred number of default episodes (expected to be right-skewed) | Incomplete portfolio
test <- subset(datCredit_real, ExclusionID==0 & DefSpell_Num >= 1 & TimeInDefSpell==1 & Partition == "Observed")[,list(Max_Epi = max(DefSpell_Num)), by=list(LoanID)]
describe(test$Max_Epi); hist(test$Max_Epi, breaks=2*length(test$Max_Epi)^(1/3))
(diag.defspell4c_real <- table(test$Max_Epi) %>% proportions())
(diag.defspell4d_real <- c(round(mean(test$Max_Epi, na.rm=T),digits=2), round(median(test$Max_Epi, na.rm=T),digits=2), 
                      max(test$Max_Epi, na.rm=T),
                      round(skewness(test$Max_Epi),digits=4), case_when(skewness(test$Max_Epi) >= 0.1 ~ "right-skewed", 
                                                                        skewness(test$Max_Epi) < 0.1 & skewness(test$Max_Epi) > -0.1 ~ "symmetric",
                                                                        skewness(test$Max_Epi) <= -0.1 ~ "left-skewed") ) )
### RESULT: Confirmed, highly right-skewed. Serves as evidence of redefault, including the degree thereof


# - 2f) Integrity tests of important control fields of default spells
# [DefSpell_Age] 
# Spell age should be capped at previous month date in 100% of right-censored cases
(diag.defspell6b_real <- datCredit_real[ExclusionID==0 & DefaultStatus1==1 & Date == maxDate_observed, .N] / 
    datCredit_real[ExclusionID==0 & DefaultStatus1==1 & Date == maxDate_observed & TimeInDefSpell == DefSpell_Age, .N] * 100)
# No records expected beyond recorded age of a default spell in the incomplete portfolio
(diag.defspell6c_real <- datCredit_real[ExclusionID==0 & DefaultStatus1==1 & Partition == "Observed" & TimeInDefSpell > DefSpell_Age, .N])
# Are all spells completely observed, incorporating right-censoring
(diag.defspell6d_real <- datCredit_real[ExclusionID==0 & DefaultStatus1==1 & Partition == "Observed" & 
                                          TimeInDefSpell == DefSpell_Age & DefSpell_LeftTrunc==0, .N] / 
    datCredit_real[ExclusionID==0 & DefaultStatus1==1 & Partition == "Observed"  & 
                     TimeInDefSpell == 1 & DefSpell_LeftTrunc==0, .N] * 100)
# This should be 100
diag.defspell6b_real == 100 & diag.defspell6c_real == 0 & diag.defspell6d_real == 100
### RESULTS: If combined test is TRUE, then integrity tests of [DefSpell_Age] passed

# [Account_Censored]
# Diagnostic for Account_Censored regarding terminal events (Shouldn't be 1)
(diag.defspell6e_real <- datCredit_real[ExclusionID==0 & Date == maxDate_observed & Age_Adj == Event_Time &
                                          Account_Censored == 1 & Event_Type %in% c("WOFF", "SETTLE", "COMPLETED"), .N] / 
    datCredit_real[ExclusionID==0 & Age_Adj == Event_Time & Date == maxDate_observed & 
                     Event_Type %in% c("WOFF", "SETTLE", "COMPLETED"), .N] * 100 )
diag.defspell6e_real == 0
### RESULTS: Proportion should be zero, confirmed if test is TRUE

# [DefSpell_Censored] + [DefSpell_Event]
# Does event and censoring indicators coincide (shouldn't)
(diag.defspell6f_real <- datCredit_real[ExclusionID==0 & DefaultStatus1==1 & TimeInDefSpell==DefSpell_Age & 
                                          DefSpell_Event == 1 & DefSpell_Censored == 1, .N])
diag.defspell6f_real == 0
### RESULT: if TRUE, then combined test is PASSED since we expect that event and censoring indicaotrs to never coincide

# [DefSpellResol_Type_Hist]
# Does resolution type correspond to censoring field?
(diag.defspell6g_real <- datCredit_real[ExclusionID==0 & DefaultStatus1==1 & TimeInDefSpell==DefSpell_Age & DefSpell_Censored == 1 & DefSpellResol_Type_Hist != "WOFF", .N] /
    datCredit_real[ExclusionID==0 & DefaultStatus1==1 & TimeInDefSpell==DefSpell_Age & DefSpell_Censored == 1, .N] * 100)
diag.defspell6g_real == 100
### RESULTS: If TRUE, then test is PASSED, thereby confirming that all [DefSpellResol_Type_Hist]-marked-as-censored records 
# correspond 100% to [DefSpell_Censored]-marked records


# - 2h) Spell-level distribution of resolution types (WOFF, Cured) for resolved defaults | Incomplete portfolio
test <- subset(datCredit_real, ExclusionID==0 & DefaultStatus1==1 & TimeInDefSpell == 1 & Partition == "Observed" & DefSpell_LeftTrunc == 0) %>% 
  mutate(ResolutionType = case_when(
    DefSpellResol_Type_Hist == "Cured" ~ "1) Cures",
    DefSpellResol_Type_Hist == "WOFF" ~ "2) Write-offs",
    DefSpellResol_Type_Hist == "Censored" ~ "3) Right-censored"))
(diag.defspell7b_real <- table(test$ResolutionType) %>% proportions())
describe(test$ResolutionType); test[is.na(ResolutionType), .N] == 0
### RESULT: Sensible split amongst outcomes: Censored: 9%, Cures: 70%, WOFF: 21%
# No missings expected, confirmed if test is TRUE, thereby confirming correct allocation of [DefSpellResol_Type_Hist]


# - 2i) Spell-level effect of right-censoring on default spells + Field logic integrity test?
(diag.defspell7c_real <- datCredit_real[ExclusionID==0 & DefaultStatus1 == 1 & TimeInDefSpell == 1 & 
                                          Partition == "Observed" & DefSpellResol_Type_Hist != "Censored", 
                                        list(Freq=.N),by=list(LoanID,DefSpell_Num)][,.N] / 
    datCredit_real[ExclusionID==0 & DefaultStatus1 == 1 & Partition == "Observed" & DefSpell_LeftTrunc == 0, 
                   list(Freq=.N),by=list(LoanID,DefSpell_Num)][,.N] * 100)
(diag.defspell7d_real <- diag.defspell7c_real == (1-diag.defspell7b_real[3])*100 ) # should be equal
### RESULT: ~ 91% of default spells are uncensored and available for LGD-modelling
# If test is TRUE, then test is PASSED, confirming integrity of co-derived fields: [DefSpellResol_Type_Hist]


# - 2j) Portfolio-level default experiences | Portfolio diagnostics
# Incomplete portfolio
test <- subset(datCredit_real, ExclusionID==0 & Partition == "Observed", 
               select=c("LoanID", "DefSpell_Num"))
test[, HasDefaults := ifelse(max(DefSpell_Num, na.rm=T) > 0, "Accounts defaulted", "Accounts never defaulted"), 
     by=list(LoanID)]
(diag.defspell8_real <- table(test$HasDefaults) %>% proportions()); rm(test); gc()
### RESULT: ~20% of accounts in the portfolio has defaulted over their lifetime, the rest hasn't


### CONTINUE: After survival modelling examination, amend both this part and the perfoming spell part
# I suspect we'll need to shift to period-level Kaplan-Meier (counting process-style) instead of spell-level,
# purely to incorporate left-truncated parts correctly.
### POSTHOC: No need, we can still adopt counting style process (thereby incorporating left-truncation) at the spell-level
# They yield equivalent estimates, evidenced previously.


# - 2k) Time-to-event analysis of default resolution (cures/write-offs)
test <- subset(datCredit_real, ExclusionID==0 & DefaultStatus1 == 1 & TimeInDefSpell==DefSpell_Age & 
                 Partition == "Observed")
# Subset again for resolved historical loans (cures/write-offs)
testb <- subset(test, DefSpellResol_Type_Hist != "Censored")
# naive summary statistics overall| time to resolution
describe(testb$DefSpell_Age); hist(testb$DefSpell_Age, breaks=2*testb[,.N]^(1/3))
# naive summary statistics by resolution type: time to write/off/cure
testb_woff <- subset(testb, DefSpellResol_Type_Hist=="WOFF")
testb_cure <- subset(testb, DefSpellResol_Type_Hist=="Cured")
describe(testb_cure$DefSpell_Age); hist(testb_cure$DefSpell_Age, breaks=2*testb_cure[,.N]^(1/3))
describe(testb_woff$DefSpell_Age); hist(testb_woff$DefSpell_Age, breaks=2*testb_woff[,.N]^(1/3))
### RESULT: Highly right-skewed distribution as expected in both cases, with longer tails for write-off 
# (mean: 30, median: 21) than for curing (mean: 19, median: 10). This makes sense given the effort
# of working-out defaults and realising the collateral

# - compute Kaplan-Meier survival estimates (product-limit) for two competing resolution types
# 1) overall
km_all_real <- survfit(Surv(time=DefSpell_Age, event=DefSpellResol_Type_Hist != "Censored", type="right") ~ 1, data=test)
summary(km_all_real)$table
# 2) Write-off
km_woff_real <- survfit(Surv(time=DefSpell_Age, event=DefSpellResol_Type_Hist=="WOFF", type="right") ~ 1, data=test)
summary(km_woff_real)$table
# 3) Cure
km_cure_real <- survfit(Surv(time=DefSpell_Age, event=DefSpellResol_Type_Hist=="Cured", type="right") ~ 1, data=test)
summary(km_cure_real)$table

# note that 1.96*sd(testb$DefSpell_Age, na.rm=T) / sqrt(testb[,.N]) constructs a 95% z-score for a confidence interval around the mean
# this corresponds to t.test(testb$DefSpell_Age)$"conf.int" 
# Adequate even for non-normal data: see https://stats.stackexchange.com/questions/112829/how-do-i-calculate-confidence-intervals-for-a-non-normal-distribution/112909
(diag.defspell9_real <- c("Right-censored %"=test[DefSpellResol_Type_Hist == "Censored",.N] / test[,.N], 
                     # 1a) naive estimates: time to resolution (cure/write-off)
                     "mean" = mean(testb$DefSpell_Age, na.rm=T), "95% CI"=1.96*sd(testb$DefSpell_Age, na.rm=T) / sqrt(testb[,.N]),
                     "median" = median(testb$DefSpell_Age, na.rm=T), "skew"=skewness(testb$DefSpell_Age, na.rm=T),
                     # 1b) survival estimates: time to resolution (cure/write-off)
                     summary(km_all_real)$table[5], "95% CI"=1.96*summary(km_all_real)$table[6], summary(km_all_real)$table[7],
                     # 2a) naive estimates: time to write-off
                     "mean write-off time" = mean(testb_woff$DefSpell_Age, na.rm=T), "95% CI (write-off)"=1.96*sd(testb_woff$DefSpell_Age, na.rm=T) / sqrt(testb_woff[,.N]),
                     "median write-off time" = median(testb_woff$DefSpell_Age, na.rm=T), "skew (write-off)"=skewness(testb_woff$DefSpell_Age, na.rm=T),
                     # 2b) survival estimates: time to write-off
                     "survival write-off"=summary(km_woff_real)$table[5], "survival write-off"=1.96*summary(km_woff_real)$table[6], "survival writeoff"=summary(km_woff_real)$table[7],
                     # 3a) naive estimates: time to cure
                     "mean cure time" = mean(testb_cure$DefSpell_Age, na.rm=T), "95% CI (cure)"=1.96*sd(testb_cure$DefSpell_Age, na.rm=T) / sqrt(testb_cure[,.N]),
                     "median cure time" = median(testb_cure$DefSpell_Age, na.rm=T), "skew (cure)"=skewness(testb_cure$DefSpell_Age, na.rm=T),
                     # 3b) survival estimates: time to cure
                     "survival cure"=summary(km_cure_real)$table[5], "survival cure"=1.96*summary(km_cure_real)$table[6], "survival cure"=summary(km_cure_real)$table[7])
)

# - Cleanup
rm(test, testb, testb_woff, testb_cure, km_all_real, km_woff_real, km_cure_real); gc()





# ------ 3. Performance spell diagnostics

# - 3a) Uniqueness of TimeInPerfSpell as a filter mechanism?
(diag.perfspell1_real <- datCredit_real[ExclusionID==0 & TimeInPerfSpell==1, list(Age_Adj), by=list(LoanID)][,.N] / 
    datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & PerfSpell_LeftTrunc==0, list(Freq = .N), by=list(LoanID, PerfSpell_Num)][,.N] * 100)
diag.perfspell1_real == 100
### RESULT: If TRUE, then test is PASSED and structure/grain is confirmed.


# - 3b) Distribution of inferred number of performance spells (expected to be right-skewed) | Incomplete portfolio
test <- subset(datCredit_real, ExclusionID==0 & !is.na(PerfSpell_Num) & 
                 TimeInPerfSpell==PerfSpell_Age & Partition == "Observed")[,list(Max_Epi = max(PerfSpell_Num)), by=list(LoanID)]
describe(test$Max_Epi); hist(test$Max_Epi, breaks=2*length(test$Max_Epi)^(1/3))
(diag.perfspell2a_real <- table(test$Max_Epi) %>% proportions())
(diag.perfspell2b_real <- c(round(mean(test$Max_Epi, na.rm=T),digits=2), round(median(test$Max_Epi, na.rm=T),digits=2), 
                       max(test$Max_Epi, na.rm=T),
                       round(skewness(test$Max_Epi),digits=4), case_when(skewness(test$Max_Epi) >= 0.1 ~ "right-skewed", 
                                                                         skewness(test$Max_Epi) < 0.1 & skewness(test$Max_Epi) > -0.1 ~ "symmetric",
                                                                         skewness(test$Max_Epi) <= -0.1 ~ "left-skewed") ) ) 
### RESULT: Confirmed, highly right-skewed. Serves as indirect evidence of redefault


# - 3c) Integrity tests of important control fields of performing spells
# [PerfSpell_Censored] + [PerfSpell_Event]
# Are all spells resolved within a completed portfolio using censored & main event fields as controllers?
(diag.perfspell3a_real <- datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & PerfSpell_LeftTrunc == 0 &
                                      (PerfSpell_Censored == 1 | PerfSpell_Event == 1), 
                                    list(Freq = .N), by=list(LoanID, PerfSpell_Num)][,.N] / 
    datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & TimeInPerfSpell == 1, list(Freq = .N), 
                   by=list(LoanID, PerfSpell_Num)][,.N] * 100)
# Does event and censoring indicators coincide (shouldn't)
(diag.perfspell3e_real <- datCredit_real[ExclusionID==0 &!is.na(PerfSpell_Num) & TimeInPerfSpell==PerfSpell_Age & 
                                      PerfSpell_Event == 1 & PerfSpell_Censored == 1, .N])
diag.perfspell3a_real == 100 & diag.perfspell3e_real == 0
### RESULT: if TRUE, then combined test is PASSED since we expect a 100% resolution rate over entire dataset and
# we expect that event and censoring indicators to never coincide.

# [TimeInPerfSpell] 
# Spell age should be capped at previous month date in 100% of right-censored cases
(diag.perfspell3b_real <- datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & Date == maxDate_observed, .N] / 
    datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & Date == maxDate_observed & TimeInPerfSpell == PerfSpell_Age, .N] * 100)
# No records expected beyond recorded age of a performing spell in the incomplete portfolio
(diag.perfspell3c_real <- datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & 
                                           Partition == "Observed" & TimeInPerfSpell > PerfSpell_Age, .N])
# Are all spells completely observed, incorporating right-censoring
(diag.perfspell3d_real <- datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & Partition == "Observed" &
                                           TimeInPerfSpell == PerfSpell_Age & PerfSpell_LeftTrunc==0, .N] / 
    datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & Partition == "Observed"  & TimeInPerfSpell == 1 &
                     PerfSpell_LeftTrunc==0, .N] * 100)
# This should be 100
diag.perfspell3b_real == 100 & diag.perfspell3c_real == 0 & diag.perfspell3d_real == 100
### RESULTS: If combined test is TRUE, then integrity tests of [DefSpell_Age] passed

# [PerfSpellResol_Type_Hist]
# Does resolution type correspond to censoring field?
(diag.perfspell3f_real <- datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & TimeInPerfSpell==PerfSpell_Age & 
                                           PerfSpell_Censored == 1 & PerfSpellResol_Type_Hist != "Defaulted", .N] /
    datCredit_real[ExclusionID==0 & !is.na(PerfSpell_Num) & TimeInPerfSpell==PerfSpell_Age & PerfSpell_Censored == 1, .N] * 100 )
diag.perfspell3f_real == 100
### RESULTS: If TRUE, then test is PASSED, thereby confirming that all [PerfSpellResol_Type_Hist]-marked-as-censored records 
# correspond 100% to [PerfSpell_Censored]-marked records


# - 3e) Spell-level distribution of resolution types (Default, Settlement) for resolved performances | Incomplete portfolio
test <- subset(datCredit_real, ExclusionID==0 & !is.na(PerfSpell_Num) & TimeInPerfSpell == PerfSpell_Age & 
                 Partition == "Observed") %>% 
  mutate(ResolutionType = case_when(
    PerfSpellResol_Type_Hist == "Defaulted"~ "1) Defaults",
    PerfSpellResol_Type_Hist == "Settled" ~ "2) Early settlements",
    PerfSpellResol_Type_Hist == "Paid-up" ~ "3) Paid-up",
    PerfSpellResol_Type_Hist == "Written-off" ~ "4) Write-offs",
    PerfSpellResol_Type_Hist == "Censored" ~ "5) Censored"))
describe(test$ResolutionType)
(diag.perfspell4b_real <- table(test$ResolutionType) %>% proportions())
test[is.na(ResolutionType), .N] == 0; rm(test); gc()
### RESULT: If TRUE, then test past and all resolution types have been allocated
# Defaults (19%); Early settlements (44%); Censored (35%), Repaid (2%), Write-offs (<0%)


# - 3f) Time-to-event analysis of performance resolution (default, settlement)
test <- subset(datCredit_real, ExclusionID==0 & !is.na(PerfSpell_Num) & TimeInPerfSpell==PerfSpell_Age & Partition == "Observed")
describe(test$PerfSpellResol_Type_Hist)
# %-prop of performance spells that are non-recurrent
test[PerfSpell_Num == 1, .N] / test[, .N] * 100
### RESULT: 91% of spells are first occurrences [including left-censored spells] vs 82% [excluding left-censoring]

# - naive estimates of time to resolution
# Subset again for resolved historical loans (cures/write-offs)
testb <- subset(test, PerfSpellResol_Type_Hist != "Censored")
# naive summary statistics overall| time to resolution
describe(testb$PerfSpell_Age); hist(testb[, PerfSpell_Age], breaks=2*testb[,.N]^(1/3))
# naive summary statistics by resolution type: time to resolution
testb_def <- subset(testb, PerfSpellResol_Type_Hist=="Defaulted")
testb_settle <- subset(testb, PerfSpellResol_Type_Hist=="Settled")
describe(testb_def$PerfSpell_Age); hist(testb_def$PerfSpell_Age, breaks=2*testb_def[,.N]^(1/3))
describe(testb_settle$PerfSpell_Age); hist(testb_settle$PerfSpell_Age, breaks=2*testb_settle[,.N]^(1/3))
### RESULT: Right-skewed distributional shape, as expected, Default distribution seems more exponentially distributed

# - compute Kaplan-Meier survival estimates (product-limit) for two competing resolution types
# 1) overall
kmPerf_all_real <- survfit(Surv(time=PerfSpell_Age, event=PerfSpellResol_Type_Hist != "Censored", type="right") ~ 1, data=test)
summary(kmPerf_all_real)$table
# 2) default
kmPerf_default_real <- survfit(Surv(time=PerfSpell_Age, event=PerfSpellResol_Type_Hist=="Defaulted", type="right") ~ 1, data=test)
summary(kmPerf_default)$table
# 3) Early Settlement
kmPerf_settle_real <- survfit(Surv(time=PerfSpell_Age, event=PerfSpellResol_Type_Hist=="Settled", type="right") ~ 1, data=test)
summary(kmPerf_settle_real)$table

# note that 1.96*sd(testb$PerfSpell_Age, na.rm=T) / sqrt(testb[,.N]) constructs a 95% z-score for a confidence interval around the mean
# this corresponds to t.test(testb$PerfSpell_Age)$"conf.int" 
# Adequate even for non-normal data: see https://stats.stackexchange.com/questions/112829/how-do-i-calculate-confidence-intervals-for-a-non-normal-distribution/112909
(diag.perfspell5_real <- c("Right-censored %"=test[PerfSpellResol_Type_Hist == "Censored",.N] / test[,.N], 
                      # 1a) naive estimates: time to resolution (default, settle)
                      "mean" = mean(testb$PerfSpell_Age, na.rm=T), "95% CI"=1.96*sd(testb$PerfSpell_Age, na.rm=T) / sqrt(testb[,.N]),
                      "median" = median(testb$PerfSpell_Age, na.rm=T), "skew"=skewness(testb$PerfSpell_Age, na.rm=T),
                      # 1b) survival estimates: time to resolution (default, settle)
                      summary(kmPerf_all_real)$table[5], "95% CI"=1.96*summary(kmPerf_all_real)$table[6], summary(kmPerf_all_real)$table[7],
                      # 2a) naive estimates: time to default
                      "mean default time" = mean(testb_def$PerfSpell_Age, na.rm=T), "95% CI (default)"=1.96*sd(testb_def$PerfSpell_Age, na.rm=T) / sqrt(testb_def[,.N]),
                      "median default time" = median(testb_def$PerfSpell_Age, na.rm=T), "skew (default)"=skewness(testb_def$PerfSpell_Age, na.rm=T),
                      # 2b) survival estimates: time to default
                      "survival default"=summary(kmPerf_default_real)$table[5], "survival default"=1.96*summary(kmPerf_default_real)$table[6], "survival default"=summary(kmPerf_default_real)$table[7],
                      # 3a) naive estimates: time to settle
                      "mean settle time" = mean(testb_settle$PerfSpell_Age, na.rm=T), "95% CI (settle)"=1.96*sd(testb_settle$PerfSpell_Age, na.rm=T) / sqrt(testb_settle[,.N]),
                      "median settle time" = median(testb_settle$PerfSpell_Age, na.rm=T), "skew (settle)"=skewness(testb_settle$PerfSpell_Age, na.rm=T),
                      # 3b) survival estimates: time to settle
                      "survival settle"=summary(kmPerf_settle_real)$table[5], "survival settle"=1.96*summary(kmPerf_settle_real)$table[6], "survival settle"=summary(kmPerf_settle_real)$table[7])
)

# - Cleanup
rm(test, testb, testb_def, testb_settle, kmPerf_all_real, kmPerf_default_real, kmPerf_settle_real)





# ------ 4. Loss rate diagnostics

# - 4a) Distributional analysis of realised loss rate (discounted to default point) | Incomplete Portfolio
testa <- subset(datCredit_real, ExclusionID==0 & DefaultStatus1==1 & TimeInDefSpell == 1 & DefSpellResol_Type_Hist == "WOFF")
describe(testa$LossRate_Real); hist(testa[LossRate_Real >=-0.1,LossRate_Real], breaks=2*testa[,.N]^(1/3))
if (testa[,.N] <= 5000){
  stat.test <- shapiro.test(testa$LossRate_Real)
} else {
  # p-values of SW-test cannot really be trusted when sample size exceeds 5000, therefore subsample as a quick-fix, despite the slight perturbation this introduces
  # see https://stats.stackexchange.com/questions/446262/can-a-sample-larger-than-5-000-data-points-be-tested-for-normality-using-shapiro
  stat.test <- shapiro.test(sample(testa$LossRate_Real, size=5000))
}
  
(diag.loss1a_real <- c(round(mean(testa$LossRate_Real, na.rm=T)*100,digits=2), round(median(testa$LossRate_Real, na.rm=T)*100,digits=2),
                  # check for negative minimum loss rates (dangerous)
                  ifelse(min(testa$LossRate_Real, na.rm=T) >= 0, 
                         # positive loss rate [black]
                         paste0(round(min(testa$LossRate_Real, na.rm=T)*100,digits=2), "\\%"),
                         # negative loss rate [red]
                         paste0("\\textcolor{red}{",round(min(testa$LossRate_Real, na.rm=T)*100,digits=2), "\\%}")
                  ),
                  # check for max loss rate > 100% (dangerous)
                  ifelse(max(testa$LossRate_Real, na.rm=T) <= 1, 
                         # positive loss rate [black]
                         paste0(round(max(testa$LossRate_Real, na.rm=T)*100,digits=2), "\\%"),
                         # negative loss rate [red]
                         paste0("\\textcolor{red}{",round(max(testa$LossRate_Real, na.rm=T)*100,digits=2), "\\%}")
                  ),
                  round(skewness(testa$LossRate_Real),digits=4), case_when(skewness(testa$LossRate_Real) >= 0.1 ~ "right-skewed", 
                                                                           skewness(testa$LossRate_Real) < 0.1 & skewness(testa$LossRate_Real) > -0.1 ~ "symmetric",
                                                                           skewness(testa$LossRate_Real) <= -0.1 ~ "left-skewed"),
                  "LossRate\\_Real", stat.test$p.value, round(1.96*sd(testa$LossRate_Real)/sqrt(testa[,.N])*100,digits=2)
))


# - 4d) Field integrity: Are all cases marked as "WOFF" at the write-off point
# NOTE: There are cases where write-off occurs without being predicated by a default spell
(diag.loss3a_real <- c(round(datCredit_real[ExclusionID==0 & WOff_Ind==1 & 
                                              (DefSpellResol_Type_Hist == "WOFF" | DefaultStatus1==0), .N] / 
                               datCredit_real[ExclusionID==0 & WOff_Ind==1, .N],digits=1) * 100,
                  "WOff\\_Ind", "DefSpellResol\\_Type\\_Hist"))
diag.loss3a_real[1] == 100
### RESULT: If test is TRUE, then test is PASSED, thereby confirming field integrity between WOff_Ind & DefSpellResol_Type

### CONTINUE: May need to refine loss definition for certain cases
# See Data Experiment 11.









### CONTINUE HERE. run/amend simulated data's diagnostics & analytics
# Then mirror on real data.



### CONTINUE: STUFF TO FIX
# 1) Minimum value in [Date_Origiation] = 1990-01-06. SAS Date convention for missings sometime


# - R Object cleanup
rm(testa, stat.test)