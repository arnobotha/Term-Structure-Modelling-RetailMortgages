# ============================ DATA EXPERIMENTS (Ancillary) =============================
# Formulate and conduct specific experiments to test hypotheses and other logic
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Credit Risk Forecasting Platform (Kasmeer-ERM)
# SCRIPT AUTHOR(S): Dr Arno Botha
# DESCRIPTION:
# Where the [QUANTS]-label is mentioned, offload experiment's outcome to quants for
# further investigation in pursuing an appropriate data treatment (or upstream fix)
# ---------------------------------------------------------------------------------------
# INPUT:
#   1) datCredit_real | raw monthly loan performance data imported in script 1b
#   2) various parameters set in the setup script 0
#
# OUTPUT:
#   1) <insight>
# =======================================================================================
# NOTE: This script predominantly comes from another project (Kasmeer).
# =======================================================================================




# ------------ Basic Data Analyses

# - Term
describe(dat.raw$LOAN_TERM); hist(dat.raw$LOAN_TERM)
### RESULTS: non-zero values only, with missingness, centered at 237.5 (median: 240), 
#   max of 970 months.

# - LN_TPE
describe(dat.raw[Counter == 1, LN_TPE])
### RESULTS: CHL (94%) "Consumer Home Loans", WHL (6%) "Wealth Home Loans"

# - status
describe(datCredit_real$status)





# ------------ Data Experiments

# ------ Data Experiment 1: Trailing zero balances | Exclusions/Treatments applied
# How prevalent is it that accounts have trailing zero-balances at the end of their 
#   histories? For terminal events (write-off, settlement, closure), cases with 
#   trailing zero-valued balances can severely contaminate the timing of such events.
#   This is material for any time-to-event modelling, e.g., default survival modelling.


# --- 0. Preliminaries

# - Sample data for testing purposes without contaminating the original dataset
test <- subset(datCredit_real, ExclusionID==0, 
               select=c("LoanID", "LN_TPE", "Date_Origination", "Date", "Counter", "Max_Counter", "Principal", 
                        "Balance", "Arrears", "Instalment", "HasWOff", "HasSettle",
                        "FurtherLoan_Ind", "HasFurtherLoan", "HasRedraw", "ExclusionID", "Age",
                        "Age_Adj", "New_Ind"))

# - Create some account-level aggregates to help with investigation
test[Principal > 0, Principal_Ratio := Balance / Principal]

# - Create lead and last-observed balance fields, whilst imputing missing-elements
# with the last known nonmissing value, using the custom function "imputeLastKnown()"
# as defined in script 0
test[ExclusionID == 0, Balance_lead1 := imputeLastKnown(shift(Balance,type="lead",n=1)), by=list(LoanID)]
test[ExclusionID == 0, Balance_last := Balance[.N], by=list(LoanID)]
test[ExclusionID == 0, Balance_secondlast := Balance[.N-1], by=list(LoanID)]
test[ExclusionID == 0, Principal_Ratio_lead1 := imputeLastKnown(shift(Principal_Ratio,type="lead",n=1)), by=list(LoanID)]
test[ExclusionID == 0, Principal_Ratio_last := Principal_Ratio[.N], by=list(LoanID)]
test[ExclusionID == 0, Principal_Ratio_secondlast := Principal_Ratio[.N-1], by=list(LoanID)]


# --- 1. Systematic threshold-based analysis: absolute vs relative thresholds for [Balance]

# - Assumption 1: Given threshold-values for [Principal_Ratio] (relative) and [Balance] (absolute),
# an account that truly has zero-valued trailing zero balances at its end should have a low 
# mean balances across that period. If an account is affected, then t_z must exist at which
# point the B(t) should be near zero across t=t_z,...,T for an account's lifetime T.


# - Define a function to iterate test logic across given thresholds
zeroBalance_Exp1_Job <- function(test, vecZeroBal_Ratio=0.01, vecZeroBal_Abs=100, describeDist=F, returnData=F) {
  
  if (length(vecZeroBal_Abs) > 1) verbose <- T else verbose <- F
  
  # - Iterate for each given threshold and collect results afterward
  for (i in 1:length(vecZeroBal_Abs)) {
    
    # - Ensure fields that will be created in this step do not exist yet
    suppressWarnings( test[, `:=`(ZeroBal_Ind1 = NULL, ZeroBal_Ind2 = NULL, ZeroBal_Ind = NULL,
                ZeroBal_Ind_lag = NULL, ZeroBal_Start = NULL, HasTrailingZeroBalances = NULL,
                Balance_Mean = NULL, ZeroBal_Length = NULL, NonZeroBal_Start = NULL,
                NonZeroBal_Mean = NULL, Contamination_Degree = NULL)] )
      
    # - Set this iteration's relative threshold beneath which an account will
    # have a so-called 'zero-valued' balance
    currThresh <- vecZeroBal_Abs[i]
    
    # - Indicate persisting zero balances across credit histories | rational/relative basis
    if (vecZeroBal_Ratio[1] > 0) {
      test[ExclusionID==0, ZeroBal_Ind1 := ifelse(Principal_Ratio <= vecZeroBal_Ratio[1] & 
                                                    Principal_Ratio_lead1 <= vecZeroBal_Ratio[1] & 
                                                    Principal_Ratio_secondlast <= vecZeroBal_Ratio[1] &
                                                    Principal_Ratio_last <= vecZeroBal_Ratio[1]
                                                  ,1,0), by=list(LoanID)] 
    }
    
    # - Indicate persisting zero balances across credit histories | absolute basis
    test[ExclusionID==0, ZeroBal_Ind2 := ifelse(Balance <= currThresh & 
                                                  Balance_lead1 <= currThresh & 
                                                  Balance_secondlast <= currThresh &
                                                  Balance_last <= currThresh
                                                  ,1,0), by=list(LoanID)]
    
    # - Create a final persisting zero-valued balance indicator by taking the min across
    # both of the previous indicators, effectively an "AND"-restriction.
    if (vecZeroBal_Ratio[1] > 0) {
      test[ExclusionID==0, ZeroBal_Ind := pmin(ZeroBal_Ind1, ZeroBal_Ind2, na.rm=T), 
           by=list(LoanID)]
    } else {
      test[ExclusionID==0, ZeroBal_Ind := ZeroBal_Ind2, by=list(LoanID)]
    }
    
    # - Lag the future values of given vector by n periods | n = 1
    # This will help in finding definitive 'regimes' of when a non-zero balance switches to
    # a zero-valued balance, as controlled by this iteration's thresholds for a 'zero-valued' balance.
    test[ExclusionID==0, ZeroBal_Ind_lag := shift(ZeroBal_Ind, type="lag", n=1), by=list(LoanID)]
    
    # - Find the starting point of zero-valued balances 
    # The logic "ZeroBal_Ind_lag==0 & ZeroBal_Ind==1" achieves this idea, particularly by
    # taking the index of where this logic is first true. However, there can be multiple 
    # such 'regimes'. Therefore, reverse the which(..), and take the index of the first TRUE-element,
    # which should correspond to the last such regime. This regime should then be truly trailing 
    # zero-valued balances.
    test[ExclusionID==0, ZeroBal_Start := ifelse(all(!is.na(ZeroBal_Ind)),
            min(coalesce(rev(which(ZeroBal_Ind_lag==0 & ZeroBal_Ind==1))[1] + 1 , 0), Max_Counter), 0),
         by=list(LoanID)]
    
    # - Indicate affected accounts for easier traversal
    test[ExclusionID==0, HasTrailingZeroBalances := 
           ifelse(ZeroBal_Start > 0, TRUE, FALSE),
         by=list(LoanID)]
    
    # - Calculate mean balance from the true ending point up to the available ending point
    test[ExclusionID==0 & HasTrailingZeroBalances==T & Counter>=ZeroBal_Start, 
         Balance_Mean := mean(Balance,na.rm=T), by=list(LoanID)]
    
    # - Calculate length of trailing zero-valued balances
    test[ExclusionID==0 & HasTrailingZeroBalances==T & Counter>=ZeroBal_Start, 
         ZeroBal_Length := Max_Counter - ZeroBal_Start + 1, by=list(LoanID)]
    
    # - Calculate the starting point of a k-month period that precedes the 
    # supposed starting point of trailing "zero-valued" balances | k = 6
    test[ExclusionID==0 & HasTrailingZeroBalances==T, 
         NonZeroBal_Start := max(1, ZeroBal_Start-6), by=list(LoanID)]
    
    # - Calculate mean balance over a k-month period that precedes the supposed starting 
    # point of trailing "zero-valued" balances.
    test[ExclusionID==0 & HasTrailingZeroBalances==T, 
         NonZeroBal_Mean := mean(Balance[max(NonZeroBal_Start,na.rm=t):max(ZeroBal_Start,na.rm=t)-1], na.rm=T), 
         by=list(LoanID)]
    
    # - Calculate contamination degree between [Balance_Mean] and [NonZeroBal_Mean] in measuring
    # aggressive zero-point adjustment
    test[ExclusionID==0 & HasTrailingZeroBalances==T,
          Contamination_Degree := Balance_Mean[max(ZeroBal_Start,na.rm=t)] / 
           (Balance_Mean[max(ZeroBal_Start,na.rm=t)] + NonZeroBal_Mean[1]),
         by=list(LoanID)]
    
    # - Mini experiment in creating an appropriate measure of "pulling zero-point back too far" (see previous)
    # let x be the zero-valued trailing balance part, and y the non-zero balances that precedes x
    #x <- mean(c(0,0,0,0)); y <- mean(c(30000,25000,20000,15000)) # ideal case
    #x <- mean(c(1000,0,0,0)); y <- mean(c(30000,25000,20000,15000,1000)) # some noise
    #x <- mean(c(7500,5000,1000,0,0,0)); y <- mean(c(30000,25000,20000,15000,7000,7500)) # some more noise
    #x <- mean(c(10000,7500,5000,1000,0,0,0)); y <- mean(c(30000,25000,20000,15000,7000,7500)) # some more noise
    #y / (x+y) # 'mixing' of mean balances ratio (between zero-valued part and non-zero-valued part) | 1 for ideal
    #x / (x+y) # degree of contamination in mean balance due to non-zero-valued balances | 0 for ideal
    
    # -- Testing cases
    # [LOOKUP] Specific non-qualifying case where balance dipped momentarily below threshold but recovered
    # Lookup <- subset(test, LoanID == 86710)
    # [LOOKUP] Specific non-qualifying case with extreme balance movements.
    # Lookup <- subset(test, LoanID == 155721)
    # [LOOKUP] Specific semi-qualifying case with extreme balance movements and noticeable "regimes"
    # Lookup <- subset(test, LoanID == 249565)
    # [LOOKUP] Specific semi-qualifying case with noticeable 'regimes' where missingness in some vectors can make life difficult
    # Lookup <- subset(test, LoanID == 155705)
    # [LOOKUP] Case that has an indicated flag
    # Lookup <- subset(test, LoanID == unique(test[HasTrailingZeroBalances==T, LoanID])[1])
    # [LOOKUP] Case without an indicated flag
    # Lookup <- subset(test, LoanID == unique(test[HasTrailingZeroBalances==F, LoanID])[1])
    # [LOOKUP] Case not detected by trailing zero indicator and where relative threshold failed
    # Lookup <- subset(test, LoanID == unique(test[Counter==(Max_Counter-1) & Balance <= 100 & HasTrailingZeroBalances==F, LoanID])[1])
    # [LOOKUP] Case with a large mean balance, contrary to 'zero-valued' trailing balances
    # Lookup <- subset(test, LoanID == unique(test[Counter==ZeroBal_Start & HasTrailingZeroBalances==T & Balance_Mean >=30000, LoanID])[1])
    
    # [DIAGNOSTIC] Prevalence of trailing zero-valued balances at the account-level?
    # Note: Counter's value doesn't matter, since ZeroBal_Start and Max_Counter does not 
    #   vary over time per loan. As long as we evaluate the expression at only a 
    #   single record / time point per loan.
    inner.diag1a <- test[ExclusionID==0 & Counter==1 & HasTrailingZeroBalances==T, .N] / 
      test[ExclusionID==0 & Counter == 1, .N]
    if (verbose) {
      cat(paste0("Given a 'zero'-balance relative threshold of ", vecZeroBal_Ratio[1]*100, 
                 "% and absolute threshold of ZAR ", currThresh,": \n\t1a. ", 
                 round(inner.diag1a*100, digits=2), "% of accounts qualifies in ",
                 "having trailing zero-valued balances at their end.\n"))
    }
    
    # [DIAGNOSTIC] Cases seemingly undetected by indicator
    inner.diag1b <- test[ExclusionID==0 & Counter==(Max_Counter-1) & Balance <= currThresh &
           HasTrailingZeroBalances==F, .N] / test[ExclusionID==0 & Counter,.N]
    if (verbose) {
      cat(paste0("\t1b. However, ", round(inner.diag1b*100,digits=2),"% of accounts were seemingly undetected by ",
                 "the new 
                 trailing-zero-valued-balance indicator. A form of measurement error.\n" ) )
    }
    
    # - Sample at the start of trailing zeros
    test_samp <- subset(test, ExclusionID==0 & Counter==ZeroBal_Start & HasTrailingZeroBalances==T)

    # [DIAGNOSTIC] Summary statistics of [Age] at the start of start of trailing 'zero-valued'
    if (describeDist & verbose) {
      describe(test_samp$Age)
      hist(test_samp[Age <= 500, Age], breaks=2*test_samp[,.N]^(1/3))
    }
    inner.diag2a <- mean(test_samp$Age, na.rm=T)
    inner.diag2b <- median(test_samp$Age, na.rm=T)
    inner.diag2c <- min(test_samp$Age, na.rm=T)
    inner.diag2d <- max(test_samp$Age, na.rm=T)
    if (verbose) {
      cat(paste0("\t2.  At t_z, the start of trailing zero-valued balances, the mean loan age is ",round(inner.diag2a), " (median: ", inner.diag2b, 
                 "), ranging between [", inner.diag2c, ", ", inner.diag2d, "].\n") )
    }
    
    # [DIAGNOSTIC] Summary statistics of [Balance] at the start of start of trailing 'zero-valued'
    if (describeDist & verbose) {
      describe(test_samp$Balance)
      hist(test_samp[Balance <= 1000 & Balance >=0, Balance], breaks=2*test_samp[,.N]^(1/3))
    }
    inner.diag3a <- mean(test_samp$Balance, na.rm=T)
    inner.diag3b <- median(test_samp$Balance, na.rm=T)
    inner.diag3c <- min(test_samp$Balance, na.rm=T)
    inner.diag3d <- max(test_samp$Balance, na.rm=T)
    if (verbose) {
      cat(paste0("\t3.  At t_z, the mean balance is ",round(inner.diag3a), " (median: ", round(inner.diag3b), 
                 "), ranging between [", round(inner.diag3c), ", ", round(inner.diag3d), "].\n") )
    }
    
    # [DIAGNOSTIC] Summary statistics of [Balance_Mean] across affected period
    if (describeDist & verbose) {
      describe(test_samp$Balance_Mean)
      hist(test_samp[Balance_Mean <= 10000 & Balance_Mean >=0, Balance_Mean], breaks=2*test_samp[,.N]^(1/3))
    }
    inner.diag4a <- mean(test_samp$Balance_Mean, na.rm=T)
    inner.diag4b <- median(test_samp$Balance_Mean, na.rm=T)
    inner.diag4c <- min(test_samp$Balance_Mean, na.rm=T)
    inner.diag4d <- max(test_samp$Balance_Mean, na.rm=T)
    if (verbose) {
      cat(paste0("\t4.  From t_z up to the end, the mean of the account-level mean balance is ",round(inner.diag4a), 
                 " (median: ", round(inner.diag4b), "), ranging between [", round(inner.diag4c), ", ", 
                 round(inner.diag4d), "].\n") )
    }
    
    # [DIAGNOSTIC] Summary statistics of [ZeroBal_Length] across affected period
    if (describeDist & verbose) {
      describe(test_samp$ZeroBal_Length)
      hist(test_samp[ZeroBal_Length <= 100, ZeroBal_Length], breaks=2*test_samp[,.N]^(1/3))
    }
    inner.diag5a <- mean(test_samp$ZeroBal_Length, na.rm=T)
    inner.diag5b <- median(test_samp$ZeroBal_Length, na.rm=T)
    inner.diag5c <- min(test_samp$ZeroBal_Length, na.rm=T)
    inner.diag5d <- max(test_samp$ZeroBal_Length, na.rm=T)
    if (verbose) {
      cat(paste0("\t5.  From t_z up to the end, the mean of the length of trailing zero-valued balances is ",round(inner.diag5a), 
                 " (median: ", round(inner.diag5b), "), 
                 ranging between [", round(inner.diag5c), ", ", round(inner.diag5d), "].\n") )
    }
    
    # [DIAGNOSTIC] Summary statistics of [Principal_Ratio] at the start of start of trailing 'zero-valued'
    if (vecZeroBal_Ratio[1] > 0) {
      if (describeDist & verbose) {
        describe(test_samp$Principal_Ratio)
        hist(test_samp[Principal_Ratio <= vecZeroBal_Ratio[1] & Principal_Ratio >= 0,Principal_Ratio], 
             breaks=2*test_samp[,.N]^(1/3))
      }
      inner.diag6a <- mean(test_samp$Principal_Ratio, na.rm=T)
      inner.diag6b <- median(test_samp$Principal_Ratio, na.rm=T)
      inner.diag6c <- min(test_samp$Principal_Ratio, na.rm=T)
      inner.diag6d <- max(test_samp$Principal_Ratio, na.rm=T)
      if (verbose) {
        cat(paste0("\t6.  At t_z, the mean balance-to-principal value is ",round(inner.diag6a*100,digits=4), 
                   "% (median: ", inner.diag6b*100, "%), 
                   ranging in [", round(inner.diag6c*100,digits=4), "%, ", round(inner.diag6d*100,digits=4), "%].\n") )
      }
    } else{
      inner.diag6a <- NA; inner.diag6b <- NA
    }
    
    # [DIAGNOSTIC] Prevalence of terminal events for these accounts with seemingly
    #   trailing zero-valued balances
    inner.diag7a <- test_samp[HasWOff==T, .N] / test_samp[,.N]
    inner.diag7b <- test_samp[HasSettle==T, .N] / test_samp[,.N]
    if (verbose) {
      cat(paste0("\t7.  Overall, of the accounts with seemingly trailing zero-valued balances, ",
                 round((inner.diag7a+inner.diag7b)*100, digits=1), "% suffered terminal events, 
                 with write-offs numbering ", round(inner.diag7a*100,digits=2), "% and early settlements numbering ",
                 round(inner.diag7b*100, digits=2), "%.\n" )  )
    }
    
    # [DIAGNOSTIC] Summary statistics of [NonZeroBal_Mean] across affected period
    if (describeDist & verbose) {
      describe(test_samp$NonZeroBal_Mean)
      hist(test_samp[NonZeroBal_Mean <= 200000, NonZeroBal_Mean], breaks=2*test_samp[,.N]^(1/3))
    }
    inner.diag8a <- mean(test_samp$NonZeroBal_Mean, na.rm=T)
    inner.diag8b <- median(test_samp$NonZeroBal_Mean, na.rm=T)
    inner.diag8c <- min(test_samp$NonZeroBal_Mean, na.rm=T)
    inner.diag8d <- max(test_samp$NonZeroBal_Mean, na.rm=T)
    if (verbose) {
      cat(paste0("\t8.  From t_z-6 up to t_z-1, the mean of the account-level mean balance is ",round(inner.diag8a), 
                 " (median: ", round(inner.diag8b), "), ranging between [", round(inner.diag8c), ", ", 
                 round(inner.diag8d), "].\n") )
    }
    
    # [DIAGNOSTIC] Summary statistics of [Contamination_Degree] given the resulting split 
    # between zero-valued balances and non-zero-valued balances
    if (describeDist & verbose) {
      describe(test_samp$Contamination_Degree)
      hist(test_samp[,Contamination_Degree], breaks=2*test_samp[,.N]^(1/3))
    }
    inner.diag9a <- mean(test_samp$Contamination_Degree, na.rm=T)
    inner.diag9b <- median(test_samp$Contamination_Degree, na.rm=T)
    inner.diag9c <- min(test_samp$Contamination_Degree, na.rm=T)
    inner.diag9d <- max(test_samp$Contamination_Degree, na.rm=T)
    if (verbose) {
      cat(paste0("\t9.  Overall, the mean degree of contamination by {Balances > Threshold} 
                 within supposedly zero-valued balances is ",round(inner.diag9a*100,digits=4), 
                 "% (median: ", inner.diag9b*100, "%), ranging in [", round(inner.diag9c*100,digits=4), 
                 "%, ", round(inner.diag9d*100,digits=4), "%].\n") )
    }
    
    
    # - prepare result set
    datIterate <- data.table("Threshold_Relative"=vecZeroBal_Ratio[1], "Threshold_Absolute" = currThresh,
               "a_Samp_ZeroBal_Prevalence"=inner.diag1a, "b_Samp_ZeroBal_Undetected"=inner.diag1b,
               "c1_ZeroBal_Age_Mean"=inner.diag2a, "c2_ZeroBal_Age_Median"=inner.diag2b,
               "d1_ZeroBal_Balance_Mean"=inner.diag3a, "d2_ZeroBal_Balance_Median"=inner.diag3b,
               "e1_ZeroBal_MeanBalance_Mean"=inner.diag4a, "e2_ZeroBal_MeanBalance_Median"=inner.diag4b,
               "f1_ZeroBal_Length_Mean"=inner.diag5a, "f2_ZeroBal_Length_Median"=inner.diag5b,
               "g1_ZeroBal_Princ_Mean"=inner.diag6a, "g2_ZeroBal_Princ_Median"=inner.diag6b,
               "h1_ZeroBal_WOff"=inner.diag7a, "h2_ZeroBal_Settled"=inner.diag7b,
               "i1_MeanNonZeroBalance_Mean"=inner.diag8a, "i2_MeanNonZeroBalance_Median"=inner.diag8b,
               "j1_MeanNonZeroBalance_Contam_Mean"=inner.diag9a, "j2_MeanNonZeroBalance_Contam_Median"=inner.diag9b)
    
    if (i==1) datResults <- datIterate else datResults <- rbind(datResults, datIterate)
  }
  
  if (returnData==F) return(datResults) else return(test)
}


# - Vector of balance-to-principal thresholds below which an account is considered as 
# having a "zero-balance". Each value will be tested for prevalence in a custom procedure
( vecZeroBal_Ratio <- c(0,0.005, 0.01, 0.015, 0.02) )

# - Vector of aboslute ZAR-valued thresholds below which an account is considered as 
# having a "zero-balance". Each value will be tested for prevalence in a custom procedure
( vecZeroBal_Abs <- c(0,10,25,50,75,100,150,200,300,400,500,750,1000,1500,2000))

# - Execute job across thresholds
ptm <- proc.time() # for runtime calculations
datZeroBal <- zeroBalance_Exp1_Job(test, vecZeroBal_Ratio[1], vecZeroBal_Abs)
datZeroBal <- rbind(datZeroBal, zeroBalance_Exp1_Job(test, vecZeroBal_Ratio[2], vecZeroBal_Abs))
datZeroBal <- rbind(datZeroBal, zeroBalance_Exp1_Job(test, vecZeroBal_Ratio[3], vecZeroBal_Abs))
datZeroBal <- rbind(datZeroBal, zeroBalance_Exp1_Job(test, vecZeroBal_Ratio[4], vecZeroBal_Abs))
datZeroBal <- rbind(datZeroBal, zeroBalance_Exp1_Job(test, vecZeroBal_Ratio[5], vecZeroBal_Abs))
proc.time() - ptm

# - Confirm prepared credit data is loaded into memory
if (!exists('datZeroBal')) unpack.ffdf(paste0(genObjPath,"DatExp1-Results"), tempPath)

# - additional bit of preparation for graphing multiple series together
datZeroBal_prep <- pivot_longer(datZeroBal, cols=3:ncol(datZeroBal), names_to="Metric", 
                                values_to="Value") %>% 
  drop_na(Value) %>% as.data.table() 

datZeroBal_prep[, Segment1 := case_when(
  Metric == "c1_ZeroBal_Age_Mean" ~ "Mean Loan Age at t_z",
  Metric == "f1_ZeroBal_Length_Mean" ~ "Mean Length of Affected History across [t_z, T ]",
  Metric == "d1_ZeroBal_Balance_Mean" ~ "At zero-point (t_z)",
  Metric == "e1_ZeroBal_MeanBalance_Mean" ~ "Beyond zero-point (t>=t_z)",
  Metric == "i1_MeanNonZeroBalance_Mean" ~ "6-months before zero-point (<t_z)",
  Metric == "h1_ZeroBal_WOff" ~ "Write-off",
  Metric == "h2_ZeroBal_Settled" ~ "Early settlement")]

# - 1. Graphing prevalence across thresholds
chosenCutoff <- 250
(g1 <- ggplot(datZeroBal_prep[Metric %in% c("a_Samp_ZeroBal_Prevalence"),], 
              aes(x=Threshold_Absolute, y=Value, colour=factor(Threshold_Relative), 
                  shape=factor(Threshold_Relative))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5, linetype="dashed") + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    labs(y="Prevalence Rate (%) of Accounts with trailing zero-valued balances", x="Absolute Balance Threshold (ZAR)") + 
    theme(legend.position = "bottom") + 
    scale_colour_brewer(palette = "Dark2", name = "Relative Balance-to-Principal Threshold (%)") + 
    scale_shape_manual(name="Relative Balance-to-Principal Threshold (%)", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=percent) + scale_x_continuous(labels=comma)
)

# - 2. Graphing means of loan ages at t_z and means of overall period length for 
#  t=t_z,...,T, both expressed in months of affected loans with
# trailing zero-valued balances, given threshold-choices
(g2 <- ggplot(datZeroBal_prep[Metric %in% c("c1_ZeroBal_Age_Mean", "f1_ZeroBal_Length_Mean"),], 
              aes(x=Threshold_Absolute, y=Value, colour=factor(Threshold_Relative), 
                  shape=factor(Threshold_Relative))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5, linetype="dashed") + 
    facet_grid(Segment1~., scales="free") + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    labs(y="Mean Value (months)", x="Absolute Balance Threshold (ZAR)") + 
    theme(legend.position = "bottom") + 
    scale_colour_brewer(palette = "Dark2", name = "Relative Balance-to-Principal Threshold (%)") + 
    scale_shape_manual(name="Relative Balance-to-Principal Threshold (%)", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma)
)

# - 3. Graphing means of balances t_z and means of overall account-level
# mean balances across t=t_z,...,T, expressed using those affected loans with
# trailing zero-valued balances, given threshold-choices
(g3 <- ggplot(datZeroBal_prep[Metric %in% c("d1_ZeroBal_Balance_Mean", 
                                            "e1_ZeroBal_MeanBalance_Mean", "i1_MeanNonZeroBalance_Mean"),], 
              aes(x=Threshold_Absolute, y=Value, colour=factor(Threshold_Relative), 
                  shape=factor(Threshold_Relative))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5, linetype="dashed") + 
    facet_grid(Segment1~., scales="free") + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    labs(y="Mean of Account-level Mean Balance (ZAR)", x="Absolute Balance Threshold (ZAR)") + 
    theme(legend.position = "bottom") + 
    scale_colour_brewer(palette = "Dark2", name = "Relative Balance-to-Principal Threshold (%)") + 
    scale_shape_manual(name="Relative Balance-to-Principal Threshold (%)", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma)
)

# - 4. Graphing means of overall account-level mean balances across t=t_z-6,...,t_z-1, 
# which is the preceding supposedly non-zero-valued balances before the zero-point.
(g4 <- ggplot(datZeroBal_prep[Metric %in% c("j1_MeanNonZeroBalance_Contam_Mean"),], 
              aes(x=Threshold_Absolute, y=Value, colour=factor(Threshold_Relative), 
                  shape=factor(Threshold_Relative))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5, linetype="dashed") + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    annotate("text", x=800, y=0.015,label="x: grand mean of account-level mean balances across [t_z,T ]", size=4)+
    annotate("text", x=850, y=0.014,label="y: grand mean of 6-month account-level mean balances for t < t_z", size=4)+
    labs(y="Mean Contamination degree (%) by non-zero balances x/(x+y)", x="Absolute Balance Threshold (ZAR)") + 
    theme(legend.position = "bottom") + 
    scale_colour_brewer(palette = "Dark2", name = "Relative Balance-to-Principal Threshold (%)") + 
    scale_shape_manual(name="Relative Balance-to-Principal Threshold (%)", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=percent) + scale_x_continuous(labels=comma)
)

# - 5. Graphing incidence rates of terminal events.
(g5 <- ggplot(datZeroBal_prep[Metric %in% c("h1_ZeroBal_WOff", "h2_ZeroBal_Settled"),], 
              aes(x=Threshold_Absolute, y=Value, colour=factor(Threshold_Relative), 
                  shape=factor(Threshold_Relative))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5, linetype="dashed") + 
    facet_grid(Segment1~., scales="free") + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    labs(y="Incidence rate (%)", x="Absolute Balance Threshold (ZAR)") + 
    theme(legend.position = "bottom") + 
    scale_colour_brewer(palette = "Dark2", name = "Relative Balance-to-Principal Threshold (%)") + 
    scale_shape_manual(name="Relative Balance-to-Principal Threshold (%)", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=percent) + scale_x_continuous(labels=comma)
)

# - Save graphs
dpi <- 150
ggsave(g1, file=paste0("exp_graphs/Exp1a-Prevalence_ZeroBal.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(g2, file=paste0("exp_graphs/Exp1b-MeanAges_ZeroBal.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(g3, file=paste0("exp_graphs/Exp1c-GrandMeanBalances_ZeroBal.png"),width=1200/dpi, height=1500/dpi,dpi=dpi)
ggsave(g4, file=paste0("exp_graphs/Exp1d-Contamination_NonZeroBal.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(g5, file=paste0("exp_graphs/Exp1e-TerminalEvents_ZeroBal.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)

# - Store experimental objects | Memory optimisation
pack.ffdf(paste0(genObjPath,"DatExp1-Results"), datZeroBal); rm(datZeroBal, datZeroBal_prep)

### RESULTS:
# 0. Let absolute ZAR-valued thresholds for [Balance] be denoted by A, and let
#   rational/relative %-valued thresholds for [Principal_Ratio] (Balance-to-Principal)
#   be denoted by R. 'Case' or 'affected account' or 'event' or 'TZB' (trailing zero-valued
#   balances) refers to those accounts considered to have trailing zero-valued balances, 
#   as controlled by A and R.
# 1. An efficient frontier of sorts manifests across the whole A-range and R-rate.
#   Little differentiation amongst the R-groups, though R=0.5% seems distinct.
#   That said, this may suggest that imposing R=0.05% may be too strict, as explained
#   by the slightly lower TZB-prevalence. However, for A<=500, there seems to be
#   very little difference across all R-thresholds. Increases in prevalence seems 
#   to slow down for A>=300, though the frontier's peak lies between A\in[75, 300].
# 2a. Peak length is reached between A\in[250,500] with very little differentiation 
#   amongst R-thresholds. However, for A>500, R-thresholds start diverging from one
#   another, albeit at different speeds. Results suggest that mean period is very
#   sensitive for A<=250, after which the process gets decidedly insensitive.
# 2b. Peak mean age reached across all R-thresholds at A=25 and massive decrease
#   that is only arrested (in speed) for A>=75. Drastic changes in the mean age
#   for A<=50 suggest volatile age distribution, supported by relatively low
#   prevalence in Graph1.
# 3a. The shape resembles a transposed (around y=0) version of Graph2b, in that
#   a trough is reached for very low A-thresholds (A=25), and increases in
#   the grand mean of 6-month mean balances preceding the zero-point (t_z)
#   are gradual from A>=75. Very little differentiation amongst R-thresholds,
#   except perhaps for R=0.5% again, as in Graph1.
# 3b. Mean balances are very low across all R-thresholds and for most of A-thresholds,
#   during which it increases monotonically. However, the speed of increases seem
#   to hasten for A=>500. However, if one supposes that small (but legitimate) balances
#   start contaminating the 'zero-valued' balances at the zero-point (t_z), then the sudden increase
#   in mean of B(t_z) for A>=500 is sensible due to this contamination. Differentiation in R-thresholds
#   also become meaningful for A>=500.
# 3c. Grand means of account-level mean balances beyond the zero-point (t_z) are very
#   small across most A-thresholds, which is sensible. Having chosen a suitable
#   zero-point (t_z), we expect the remainder of balances beyond that time point to
#   be very small and close to zero. Similar to Graph3b, grand means increase
#   noticeably faster for A>=500, coupled with increasing divergence amongst R-thresholds.
#   For A<500, the grand means are not meaningfully differentiated by any R-threshold.
# 4.  Contamination degree (%) is volatile for very small for A<=25, which suggests low sampling
#   sizes, supported by Graph1. Contamination reaches a plateau of sorts for A\in[75, 300], after which
#   point, contamination increases gradually. For A>=500, contamination diverges across
#   R-thresholds while for A<400, there is no point in using R-thresholds (contamination remains the same).
# 5a. Within the TZB-cases, the early settlement rate increases astronomically from A=0 up to A=75, which
#   suggests a great deal of turbulence, supported by low prevalence in Graph1. A long plateau is reached
#   for A>=100 with very little differentiation amongst all R-thresholds.
# 5b. Write-off rates within TZB-cases also increase markedly for A<=100, after which it slows down.
#   Using write-off incidence as a guiding light is sensible given that TZB-cases will ultimately
#   affect the timing of recorded write-off. Little differentiation amongst R-thresholds for A<=300.
#   Platuaus are reached for each R-threshold for A>=750, after which the write-off rate doesn't change
#   maningfully, except for R=0.

### CONCLUSION:
# Specifying A=250 seems to isolate TZB-cases reasonably well in that those means are not
# overly contaminated by non-zero but bigger balances, which will otherwise be discarded
# if the zero-point is too big. Conversely, the chosen A-threshold will also not ignore
# persisting zero valued-balances and/or fail to discard truly zero-valued balances.
# This threshold is supported by supposedly mistimed write-off and early settlement incidence
# rates that would anyways plateau with greater A-thresholds. Lastly, the case can be made for 
# a larger A-threshold, though this will need to be accompanied by an R-threshold, which adds
# complexity and may be invalid for very large mortgages in any case. A sufficiently small A-threshold
# circumvents the need for an R-threshold since the results show rough equivalence between A and R-values
# for A<=250.



# --- 2. Threshold-sensitive distributional analyses given specific thresholds

# -- 0. Preliminaries
currThresh <- 250

# - Sample data for testing purposes without contaminating the original dataset
test <- subset(datCredit_real, ExclusionID==0, 
               select=c("LoanID", "LN_TPE", "Date_Origination", "Date", "Counter", "Max_Counter", "Principal", 
                        "Balance", "Arrears", "Instalment", "HasWOff", "HasSettle",
                        "FurtherLoan_Ind", "HasFurtherLoan", "HasRedraw", "ExclusionID", "Age",
                        "Age_Adj", "New_Ind"))

# - Create lead and last-observed balance fields, whilst imputing missing-elements
# with the last known nonmissing value, using the custom function "imputeLastKnown()"
# as defined in script 0
test[ExclusionID == 0, Balance_lead1 := imputeLastKnown(shift(Balance,type="lead",n=1)), by=list(LoanID)]
test[ExclusionID == 0, Balance_last := Balance[.N], by=list(LoanID)]
test[ExclusionID == 0, Balance_secondlast := Balance[.N-1], by=list(LoanID)]
test[ExclusionID == 0, LastDate := Date[.N], by=list(LoanID)]

# - Now impose selected thresholds
test <- zeroBalance_Exp1_Job(test, vecZeroBal_Ratio=0, vecZeroBal_Abs=currThresh, returnData=T)

# - Sample at the start of trailing zeros or end of the account
test_samp <- subset(test, ExclusionID==0 & Counter==ZeroBal_Start & HasTrailingZeroBalances==T)

# -- 1. Analyses

# [DIAGNOSTIC] Prevalence of trailing zero-valued balances at the account-level?
(diag.real1a <- test[ExclusionID==0 & Counter==1 & HasTrailingZeroBalances==T, .N] / 
  test[ExclusionID==0 & Counter == 1, .N] * 100)
### RESULTS: 22.43% of accounts have "Trailing Zero-valued Balances" (TZB).

# [DIAGNOSTIC] Cases seemingly undetected by indicator
(diag.real1b <- test[ExclusionID==0 & Counter==(Max_Counter-1) & Balance <= currThresh &
                      HasTrailingZeroBalances==F, .N] / test[ExclusionID==0 & Counter,.N] *100 )
### RESULTS: 0.015% of accounts were failed to be detected by the decision rule. Negligible.

# [DIAGNOSTIC] Prevalence of terminal events for these accounts with TZBs
(diag.real2a <- test_samp[HasWOff==T, .N] / test_samp[,.N] * 100 )
(diag.real2b <- test_samp[HasSettle==T, .N] / test_samp[,.N] * 100 )
(diag.real2a + diag.real2b)
### RESULTS: Of those TZB-cases, ~86% have a pending terminal-event,
# which gives partial assurance on the idea of cutting out these seemingly
# useless TZB-parts from credit histories.

# [DIAGNOSTIC] Right-censorship of those TZB-cases without a terminal event
test_samp_b <- subset(test_samp, HasWOff==F & HasSettle==F, select=c("LastDate"))
describe(test_samp_b$LastDate); barplot(prop.table(table(test_samp_b$LastDate)))
### RESULTS: The remaining TZB-cases without a terminal event are universally
# right-censored. This is sensible and assures the broader idea of 'shortening'
# the run-up to a terminal-event, that is otherwise introduced by lagged system processes.

# [DIAGNOSTIC] Summary statistics of [Age] at the start of start of TZB-history
describe(test_samp$Age)
hist(test_samp[Age <= 500, Age], breaks=2*test_samp[,.N]^(1/3))
### RESULTS: Irregular but largely right-skewed distribution, almost Weibull-shaped.
# Has a large spike in values at around 240. And of course, the vast majority of loan accounts
# are 240-month mortgages, which can corroborate this spike. Mean loan age at the start of TZB
# is 141 months (median: 136), with some large extreme outliers (max: 1426).

# [DIAGNOSTIC] Summary statistics of [Balance] at the start of start of trailing 'zero-valued'
describe(test_samp$Balance)
hist(test_samp[Balance <= 1000 & Balance >=0, Balance], breaks=2*test_samp[,.N]^(1/3))
### RESULTS: Highly right-skewed distribution with mean balance of 5 (median: 0).
# Sensible for the point at which we consider the balance as 'zero-valued'.

# [DIAGNOSTIC] Summary statistics of [Balance_Mean] across affected period
describe(test_samp$Balance_Mean)
hist(test_samp[Balance_Mean <= 60 & Balance_Mean >=0, Balance_Mean], breaks=2*test_samp[,.N]^(1/3))
### RESULTS: Mean balances across the period of TZB-history are very small and 
# intensely right-skewed, which is to be expected. Grand mean of ~4 (median: 0), with some
# extreme outlying mean balances, though still very small in absolute terms (max:99).

# [DIAGNOSTIC] Summary statistics of [ZeroBal_Length] across affected period
describe(test_samp$ZeroBal_Length)
hist(test_samp[ZeroBal_Length <= 100, ZeroBal_Length], breaks=2*test_samp[,.N]^(1/3))
### RESULTS: Modulated right-skewed distribution of TZB-lengths with a mean of
# ~18 months (median: 5) that will effectively be removed from the tail-ends of
# credit-histories from those TZB-affected cases.

### CONCLUSION: Selected threshold deemed reasonable and safe.

# - cleanup
rm(test, test_samp, test_samp_b)



# --- 3. Towards a treatment: Treating trailing zero-valued balances

# - ZAR-valued threshold for [Balance]
currThresh <- 250 

cat(paste0("Detecting so-called Trailing Zero-valued Balances (TZB) cases with ZAR <= ",
           currThresh, "..\n"))

# - Sample data for testing purposes without contaminating the original dataset
test <- subset(datCredit_real, ExclusionID==0, 
               select=c("LoanID", "Date", "Counter", "Max_Counter", "Balance",
                        "HasWOff", "HasSettle", "ExclusionID", "TreatmentID"))

# - Create lead and last-observed balance fields, whilst imputing missing-elements
# with the last known nonmissing value, using the custom function "imputeLastKnown()"
# as defined in script 0
test[ExclusionID == 0, Balance_lead1 := 
       imputeLastKnown(shift(Balance,type="lead",n=1)), by=list(LoanID)]
test[ExclusionID == 0, Balance_last := Balance[.N], by=list(LoanID)]
test[ExclusionID == 0, Balance_secondlast := Balance[.N-1], by=list(LoanID)]
test[ExclusionID == 0, LastDate := Date[.N], by=list(LoanID)]

# - Indicate persisting zero balances across credit histories | absolute basis
test[ExclusionID==0, ZeroBal_Ind := 
       ifelse(Balance <= currThresh & Balance_lead1 <= currThresh & 
                Balance_secondlast <= currThresh &
                Balance_last <= currThresh,1,0), by=list(LoanID)]

# - Lag the future values of given vector by n periods | n = 1
# This will help in finding definitive 'regimes' of when a non-zero balance switches to
# a zero-valued balance, as controlled by this iteration's thresholds for a 'zero-valued' balance.
test[ExclusionID==0, ZeroBal_Ind_lag := 
       shift(ZeroBal_Ind, type="lag", n=1), by=list(LoanID)]

# - Find the starting point of where [ZeroBal_Remain_Ind] became (and stayed) 1, 
# The logic "ZeroBal_Ind_lag==0 & ZeroBal_Ind==1" achieves this idea, particularly by
# taking the index of where this logic is first true. However, there can be multiple 
# such 'regimes'. Therefore, reverse the which(..), and take the index of the first TRUE-element,
# which should correspond to the last such regime. This regime should then be truly trailing 
# zero-valued balances.
test[ExclusionID==0, ZeroBal_Start := 
       ifelse(all(!is.na(ZeroBal_Ind)),
              min(coalesce(rev(which(ZeroBal_Ind_lag==0 & ZeroBal_Ind==1))[1] + 1 , 0),
                  Max_Counter), 0), by=list(LoanID)]

# - Indicate affected accounts for easier traversal
test[ExclusionID==0, HasTrailingZeroBalances := 
       ifelse(ZeroBal_Start > 0, TRUE, FALSE),
     by=list(LoanID)]

# [DIAGNOSTIC] Prevalence of accounts with trailing zero-valued balances
# that suffered a terminal event later in loan life?
(diag.real1a <- test[ExclusionID==0 & Counter==1 & HasTrailingZeroBalances==T & 
                       (HasWOff==T | HasSettle==T), .N] / 
    test[ExclusionID==0 & Counter == 1, .N] * 100)
(diag.real1a_rec <- test[ExclusionID==0 & HasTrailingZeroBalances==T & (HasWOff==T | HasSettle==T) &
                           Counter>=ZeroBal_Start, .N] / test[ExclusionID==0, .N] * 100)

# - Conditional treatment
if (diag.real1a > 0) {
  
  # - Calculate mean balance from the true ending point up to the available ending point
  test[ExclusionID==0 & HasTrailingZeroBalances==T & Counter>=ZeroBal_Start, 
       Balance_Mean := mean(Balance,na.rm=T), by=list(LoanID)]
  
  # - Calculate length of trailing zero-valued balances
  test[ExclusionID==0 & HasTrailingZeroBalances==T & Counter>=ZeroBal_Start, 
       ZeroBal_Length := Max_Counter - ZeroBal_Start + 1, by=list(LoanID)]
  
  # - Sample at the start of trailing zeros or end of the account
  test_samp <- subset(test, ExclusionID==0 & Counter==ZeroBal_Start & HasTrailingZeroBalances==T)
  
  # [DIAGNOSTIC] Prevalence of terminal events for these accounts with TZBs
  (diag.real1_2a <- test_samp[HasWOff==T, .N] / test_samp[,.N] * 100 )
  (diag.real1_2b <- test_samp[HasSettle==T, .N] / test_samp[,.N] * 100 )
  diag.real1_2 <- diag.real1_2a + diag.real1_2b
  
  # [DIAGNOSTIC] Right-censorship of those TZB-cases without a terminal event
  (diag.real1_2c <- test_samp[HasWOff==F & HasSettle==F & LastDate == maxDate_observed, .N] / 
    test_samp[HasWOff==F & HasSettle==F, .N] * 100 )
  
  # [DIAGNOSTIC] Mean of [Balance_Mean], should be very low
  ( diag.real1_3a <- mean(test_samp$Balance_Mean, na.rm=T) )
  ( diag.real1_3b <- median(test_samp$Balance_Mean, na.rm=T) )
  
  # [DIAGNOSTIC] Mean of [ZeroBal_Length]
  ( diag.real1_4a <- mean(test_samp$ZeroBal_Length, na.rm=T) )
  ( diag.real1_4b <- median(test_samp$ZeroBal_Length, na.rm=T) )
  
  cat("TZB-cases detected in ", round(diag.real1a,digits=1), 
      "% of accounts, of which ", round(diag.real1_2,digits=1), "% suffered a terminal event.
      Of those accounts that did not terminate, ", round(diag.real1_2c,digits=1), 
      "% are right-censored at the study-end [", format(maxDate_observed, "%d-%b-%Y"), "]. 
      The grand mean balance across TZB-histories is ZAR", round(diag.real1_3a,digits=0), " (median:",
      round(diag.real1_3b,digits=0), ") and the mean length of TZB-histories 
      is", round(diag.real1_4a,digits=0),"months (median:", round(diag.real1_4b,digits=0), 
      "), which will be removed during treatment after moving terminal events earlier. 
      \tEffecting treatment ...\n")
  
  # - Cleanup in optimising memory use
  rm(test_samp)

  # - Fuse account-level indicator back into main longitudinal dataset | should be main in actual treatment
  testb <- merge(datCredit_real, test[,list(LoanID, Counter, HasTrailingZeroBalances, ZeroBal_Start)], 
                 by=c("LoanID", "Counter"), all.x=T)
  
  # [TREATMENT] Move terminal event incidence earlier | Write-off
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & HasWOff==T, 
        WOff_Ind := ifelse(Counter < (ZeroBal_Start-1), WOff_Ind, WOff_Ind[.N]), by=list(LoanID)]
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & HasWOff==T, 
        WriteOff_Amt := ifelse(Counter < (ZeroBal_Start-1), WriteOff_Amt, WriteOff_Amt[.N]), by=list(LoanID)]
  
  # [TREATMENT] Move terminal event incidence earlier | Early Settlement / Closure
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & HasSettle==T, 
        EarlySettle_Ind := ifelse(Counter < (ZeroBal_Start-1), EarlySettle_Ind, EarlySettle_Ind[.N]), by=list(LoanID)]
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & HasSettle==T, 
        EarlySettle_Amt := ifelse(Counter < (ZeroBal_Start-1), EarlySettle_Amt, EarlySettle_Amt[.N]), by=list(LoanID)]
  
  # [TREATMENT] Re-calculate affected aggregates
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & (HasWOff==T | HasSettle==T) & Counter < ZeroBal_Start, 
        Max_Counter := .N, by=list(LoanID)]

  # - Mark affected accounts with current treatment ID
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & (HasWOff==T | HasSettle==T), 
        TreatmentID := paste0(TreatmentID, ";5")]
  
  # [TREATMENT] Mark the affected histories as an Exclusion
  testb[ExclusionID==0 & HasTrailingZeroBalances==T & (HasWOff==T | HasSettle==T) & 
          Counter>=ZeroBal_Start, ExclusionID := 6]
  
  # - Create and add exclusion impact to a common table
  datExclusions <- rbind(datExclusions, 
                         data.table("Excl_ID"=5, 
                                    "Reason"=paste0("Trailing zero-valued balances [", currThresh, "]"),
                                    "Impact_Account" = diag.real1a, "Impact_Dataset" = diag.real1a_rec) )
  
  check_advTreat5 <- testb[ExclusionID==0 & HasTrailingZeroBalances==T & (HasWOff==T | HasSettle==T) & 
                              Counter==Max_Counter, .N] / testb[ExclusionID==0 & Counter == 1, .N] * 100 == diag.real1a
  cat( check_advTreat5 %?% paste0('SAFE: Accounts with trailing zero-valued balances (<=', currThresh, 
                                   ') and a terminal event were successfully treated.\n') %:% 
         paste0('SAFE: Failed to treat accounts with trailing zero-valued balances (<= ', currThresh, 
                ') and a terminal event.\n'))
  
  rm(testb, test)
}


# ------ END-RESULT: This experiment led to Advanced Data Treatment 5 and associated Exclusion 6





# ------ Data Experiment 2: Terminal event vectors & associated amounts | Exclusions applied

# --- 0. Preliminaries

# - Sample data for testing purposes without contaminating the original dataset
test <- subset(datCredit_real, ExclusionID ==0,
               select=c('LoanID', 'Date', 'Counter', 'Max_Counter', 'ExclusionID', 'TreatmentID', 
                        'Principal', 'Balance', 'Arrears', 'HasSettle', 'EarlySettle_Ind', 'EarlySettle_Amt',
                        'HasWOff', 'WOff_Ind', 'WOFF_DATE', 'WriteOff_Amt', 'CLS_STAMP','LN_TPE', "RMNG_TERM"))

# - Account-level aggregates & flags to facilitate specific tests
test[, HasWOff := max(WOff_Ind, na.rm=T), by=list(LoanID)]
test[, HasWOffDate := ifelse(any(!is.na(WOFF_DATE)), 1, 0) ,by=list(LoanID)]
test[, Max_WOff_Amt := max(WriteOff_Amt, na.rm=T), by=list(LoanID)]
test[, HasSettle := max(EarlySettle_Ind, na.rm=T), by=list(LoanID)]
test[, HasClosure := ifelse(any(!is.na(CLS_STAMP)), 1, 0), by=list(LoanID)]
test[, HasClosureAtLastRec := ifelse(!is.na(CLS_STAMP[.N]), 1, 0), by=list(LoanID)]


# --- 1. [WOff_Ind] & [WriteOff_Amt] & [WOFF_DATE]

# [DIAGNOSTIC] When 0 is indicated, then the associated amount must be 0, and vice versa. 
#   How often does logic prevail?
testb <- subset(test, WOff_Ind == 0)
describe(testb$WriteOff_Amt); hist(testb$WriteOff_Amt, breaks=2*testb[,.N]^(1/3))
(diag.real2_1a <- mean(testb$WriteOff_Amt, na.rm=T))
diag.real2_1a == 0 # Should be TRUE
### Results: Relatively small mean and zero median. But there are some very large outliers
# (only positive). Worrisome.
# [POST HOC] Issue resolved via cleaning

# [SANITY CHECK] How often is the write-off amount non-zero (when it should be zero)?
(diag.real2_1b <- datCredit_real[WOff_Ind == 0 & WriteOff_Amt != 0 & ExclusionID ==0, .N] / 
    datCredit_real[WOff_Ind == 0 & ExclusionID ==0, .N] * 100)
diag.real2_1b == 0 # Should be TRUE
### Results: ~ 4% of cases. Not material, but still warrants a bit of data cleaning
# [POST HOC] Issue resolved via cleaning

# [DIAGNOSTIC] When [WOff_Ind]==1, then the associated amount must > 0. Prevalence?
testc <- subset(datCredit_real, WOff_Ind == 1)
describe(testc$WriteOff_Amt); hist(testc$WriteOff_Amt)
(diag.real2_1c <- mean(testc$WriteOff_Amt, na.rm=T))
diag.real2_1c != 0 # Should be TRUE
### RESULTS: All cases are non-zero, as expected. Right-skewed distribution with a 
#   mean ~244k (median: 141k). Has some extreme outliers (expected), but curiously, 
#   has very small values (near 0) as well.

# [SANITY CHECK] How often is the write-off amount zero (when it should be non-zero) ?
(diag.real2_1d <-  datCredit_real[WOff_Ind == 1 & WriteOff_Amt == 0 & ExclusionID==0, .N] / 
    datCredit_real[WOff_Ind == 1 & ExclusionID==0, .N] * 100 )
diag.real2_1d == 0 # Should be TRUE
### RESULTS: 0% of cases, as expected.

# [DIAGNOSTIC] Account-level correlation between [WOff_Ind] (via [HasWOff]) and [HasWOffDate] ?
(diag.real2_1e <- prop.table(table(test[Counter == 1,list(HasWOff, HasWOffDate)])))
diag.real2_1e[1,2] == 0 & diag.real2_1e[2,1] == 0 # Should be TRUE
### RESULTS: 95% of accounts are zero-valued in both [HasWOff] and [HasWOffDate] whilst
# ~5% of accounts are 1-valued in both. Perfectly correlated. Safe.
# As such, [HasWOffDate] and the underlying [WOFF_DATE] can be dropped.


# --- 2. [EarlySettle_Ind] & [EarlySettle_Amt] & [CLS_STAMP]
# [DIAGNOSTIC] When 0 is indicated, then the associated amount must be 0, and vice versa. 
#   How often does this logic prevail?
testb <- subset(test, EarlySettle_Ind == 0)
describe(testb$EarlySettle_Amt); hist(testb$EarlySettle_Amt, breaks=2*testb[,.N]^(1/3))
(diag.real2_2a <- mean(testb$EarlySettle_Amt, na.rm=T))
diag.real2_2a == 0 # Should be TRUE
### Results: Universally zero, a good thing.

# [SANITY CHECK] How often is the settlement amount non-zero (when it should be zero)?
(diag.real2_2b <- datCredit_real[ExclusionID == 0 & EarlySettle_Ind == 0 & EarlySettle_Amt != 0, .N] /
    datCredit_real[ExclusionID == 0 & EarlySettle_Ind == 0, .N] * 100)
diag.real2_2b == 0 # Should be TRUE
### Results: 0% of cases, a good thing.

# [DIAGNOSTIC] When [EarlySettle_Ind]==1, then the settlement must > 0. Prevalence?
testc <- subset(test, EarlySettle_Ind == 1)
describe(testc$EarlySettle_Amt); hist(testc$EarlySettle_Amt)
(diag.real2_2c <- mean(testc$EarlySettle_Amt, na.rm=T))
diag.real2_2c != 0 # Should be TRUE
### RESULTS: Most cases are non-zero, as expected. Highly right-skewed amount.

# [DIAGNOSTIC] Account-level correlation between [HasSettle] and [HasClosure] ?
(diag.real2_2e <- prop.table(table(test[Counter == 1,list(HasSettle, HasClosure)])))
diag.real2_2e[1,2] == 0 & diag.real2_2e[2,1] == 0 # Should be TRUE
### RESULTS: 40% of accounts are zero-valued in both [HasSettle] and [HasClosure] 
#   but there  is ~5.5% of accounts that have a closure date missed by the main 
#   [EarlySettle_Ind] field. Requires refinement. On the other hand, ~! 55%% of 
#   accounts are 1-valued in both, with no mismatches.
testd <- subset(test, Counter == 1 & HasSettle == 0 & HasClosure == 1)
describe(testd$HasWOff)
### RESULTS: Of those 5.5% affected cases, ~91% have a write-off indicator, which suggests 
#   that [CLS_DATE] not only services early settlement but also cases of write-off, 
#   even though this does not hold universally.
# [POST HOC] Issue resolved. 100% were written-off post 

# [DIAGNOSTIC] Account-level prevalence of closures missed by [EarlySettle_Ind] or [WOff_Ind] flags?
(diag.real2_2f <- test[ExclusionID==0 & Counter==1 & HasClosure==1 & HasSettle==0 & HasWOff==0, .N] / 
    test[ExclusionID==0 & Counter==1, .N] *100 )
(diag.real2_2f_closures <- test[ExclusionID==0 & Counter==1 & HasClosure==1 & HasSettle==0 & HasWOff==0, .N] / 
    test[ExclusionID==0 & Counter==1 & (HasClosure==1 | HasSettle==1 | HasWOff==1), .N] *100 )
(diag.real2_2f_rec <- test[ExclusionID==0 & HasClosure==1 & HasSettle==0 & HasWOff==0, .N] / 
    test[ExclusionID==0, .N] *100 )
### RESULTS: ~ 0.5% of accounts are affected (0.8% of terminated accounts), but this affects 0.66% 
#   of all records, which is a bit more significant than the account-level prevalence. 
#   Given that these are rare terminal event-related cases, which we'll need to model, this issue 
#   warrants treatment rather than exclusion.
# [LOOKUP] Case with a closure data but no corresponding settlement or write-off flag | LoanID == 9787
# [POST HOC] issue resolved (zero-counts as expected) after cleaning.
lookup <- subset(datCredit_real, LoanID == 
                   unique(test[Counter==1 & HasClosure==1 & HasSettle==0 & HasWOff==0, LoanID])[1])

# [DIAGNOSTIC] Does [HasClosure] only coincide with the last record as one would expect?
(diag.real2_2g <- test[Counter==1 & HasClosureAtLastRec==T, .N] / 
    test[Counter==1 & HasClosure==T, .N] * 100)
### RESULTS: Yes, in most cases. Deemed immaterial
# [LOOKUP] Case with closure-signal earlier than last record
lookup <- subset(datCredit_real, LoanID == 
                   unique(test[Counter==1 & HasClosure==1 & HasClosureAtLastRec==0, LoanID])[1])

# [DIAGNOSTIC] Prevalence of early settlements coinciding with study-end?
(diag.real2_2g <- test[ExclusionID==0 & Counter==Max_Counter & EarlySettle_Ind == 1 & Date == maxDate_observed, .N] /
    test[ExclusionID==0 & Counter==Max_Counter & Date == maxDate_observed, .N] * 100)
### RESULTS: Only ~ 0.6% of presumably right-censored cases at the study-end were settled 'early'.
#   This reassures that the early settlement indicator was intelligently created


# --- 3. Towards a treatment: Closures unaccounted for by [EarlySettle_Ind] or [Woff_Ind]
# - Isolate affected cases
testd <- subset(test, Counter == Max_Counter & HasClosure == 1 & HasSettle == 0 & HasWOff == 0)

# - [DIAGNOSTIC-GROUP] Of the affected, how are some characteristics distributed?

# [LN_TPE]
describe(testd$LN_TPE)
### RESULTS: CHL 81% vs WHL 19%. Materially different distribution than population-level analysis

# [Arrears]
describe(testd$Arrears); hist(testd$Arrears, breaks=2*testd[,.N]^(1/3))
### RESULTS Majority of values are zeros (95% percentile still zero, e.g.)

# [Principal]
describe(testd$Principal); hist(testd$Principal, breaks=2*testd[,.N]^(1/3))
### RESULTS: Highly right-skewed distribution as expected. However, the degree of zero-valued 
#   principals is intriguing. Might this be correlated to zero-valued balances as well? 
#   (Can't use Pearson since data is clearly non-normal).

# [Balance]
describe(testd$Balance); hist(testd$Balance, breaks=2*testd[,.N]^(1/3))
testd[Balance <= 250, .N] / testd[,.N] * 100 # Same cut-off as in TZB-cases - see Data Experiment 1
testd[Balance <= 250, .N] / testd[is.na(RMNG_TERM),.N] * 100 # 98.6%
### RESULTS: Also highly right-skewed distribution.
# But ~80% of cases had a 'zero'-valued balance, which supports the idea of a 
# natural paid-up status instead of early settlement. However, when deselecting
# those with missingness in [RMNG_TERM], this proportion increases significantly to ~99%.

# [RMNG_TERM]
describe(testd$RMNG_TERM); hist(testd$RMNG_TERM, breaks=2*testd[,.N]^(1/3))
testd[!is.na(RMNG_TERM), .N] / testd[, .N] * 100 # 18% non-missing/usable
testd[RMNG_TERM <= 1, .N] / testd[!is.na(RMNG_TERM), .N] * 100 # ~96%
testd[RMNG_TERM <= 1 & Balance <= 250, .N] / testd[!is.na(RMNG_TERM) & RMNG_TERM <= 1, .N] * 100 # 100%
testd[(Balance <= 250 & is.na(RMNG_TERM)) | (!is.na(RMNG_TERM) & RMNG_TERM <= 1 & Date < maxDate_observed), .N] / 
  testd[, .N] * 100 # ~80%
### RESULTS: Also highly right-skewed distribution. Has missings.
# Missingness imply only partially using this field when available (but 18% of cases)
# But, positively, of the non-missing cases, ~96% of cases had [RMNG_TERM] <= 1
# Of this, 100% had a Balance <= cut-off. Both of these facts support the idea of
# a natural paid-up status instead of early settlement.
# Lastly, if testing for either zero-value balance or small remaining term,
# then ~80% of unflagged cases are feasibly paid-up events

# [Principal]-[Balance] correlation
ggscatter(testd, x="Principal", y="Balance", cor.method="kendall", cor.coef=T, add="reg.line", conf.int=T)
### RESULTS: Significant correlation


### CONCLUSION: Given the high prevalence of zero-valued arrears and balances, with 
#   very low values in remaining term, with no clear indication of 
#   write-off events, one can assume that the majority of these cases are natural
#   paid-up events.


# - Create repaid flag for unflagged account closures
test[, Repaid_Ind :=
     ifelse(Counter==Max_Counter & HasClosure==1 & HasSettle==0 & HasWOff == 0 &
              ( (Balance <= 250 & is.na(RMNG_TERM)) | 
              (!is.na(RMNG_TERM) & RMNG_TERM <= 1 & Date < maxDate_observed)), 1,0),
     by=list(LoanID)]
test[, HasRepaid := max(Repaid_Ind, na.rm=T), by=list(LoanID)]
test[Repaid_Ind == 1, .N] / test[Counter==Max_Counter & HasClosure==1 & HasSettle==0 & HasWOff == 0, .N] * 100
### RESULTS: Corresponds exactly to previous %. Correctly isolated therefore.

# [DIAGNOSTIC] Account-level prevalence of closures missed by [EarlySettle_Ind] or [WOff_Ind] or [Repaid_Ind] flags?
(diag.real2_2h <- test[ExclusionID==0 & Counter==1 & HasClosure==1 & HasSettle==0 & HasWOff==0 & HasRepaid==0, .N] / 
    test[ExclusionID==0 & Counter==1, .N] *100 )
(diag.real2_2h_closures <- test[ExclusionID==0 & Counter==1 & HasClosure==1 & HasSettle==0 & HasWOff==0 & HasRepaid==0, .N] / 
    test[ExclusionID==0 & Counter==1 & (HasClosure==1 | HasSettle==1 | HasWOff==1 | HasRepaid==1), .N] *100 )
(diag.real2_2h_rec <- test[ExclusionID==0 & HasClosure==1 & HasSettle==0 & HasWOff==0 & HasRepaid==0, .N] / 
    test[ExclusionID==0, .N] *100 )

# [TREATMENT] Mark the remaining unaccounted closures as early settlements
# while redistributing any non-zero balances
# NOTE: The associated [EarlySettle_Amt] and associated fields will be curated later
test[ExclusionID==0 & Counter==Max_Counter & HasClosure==1 & HasSettle==0 & 
       HasWOff == 0 & HasRepaid == 0, EarlySettle_Ind := 1]
test[ExclusionID==0 & Counter==Max_Counter & HasClosure==1 & HasSettle==0 & 
       HasWOff == 0 & HasRepaid == 0, EarlySettle_Amt := EarlySettle_Amt + Balance]
test[ExclusionID==0 & Counter==Max_Counter & HasClosure==1 & HasSettle==0 & 
       HasWOff == 0 & HasRepaid == 0, Balance := 0]

test[, HasSettle := max(EarlySettle_Ind, na.rm=T), by=list(LoanID)] # Re-calculate given treatment

# [SANITY CHECK] Account-level prevalence of closures missed by [EarlySettle_Ind] or [WOff_Ind] flags?
(check_advTreat1 <- test[ExclusionID==0 & Counter==1 & HasClosure==1 & HasSettle==0 & HasWOff==0 & HasRepaid==0, .N] / 
    test[ExclusionID==0 & Counter==1, .N] *100 == 0) # Should be true
cat( check_advTreat1 %?% 'SAFE: Unaccounted closures succesfully treated as early settlements.\n' %:% 
       'WARNING: Treating unaccounted closures as early settlements failed.\n')


# --- 4. Towards a treatment 8: Conforming event amounts to the balance at termination

# - Re-sample data for testing purposes without contaminating the original dataset
test <- subset(datCredit_real, ExclusionID ==0,
               select=c('LoanID', 'Date', 'Counter', 'Max_Counter', 'ExclusionID', 'TreatmentID', 
                        'Principal', 'Balance', 'Arrears', 'HasSettle', 'EarlySettle_Ind', 'EarlySettle_Amt',
                        'HasWOff', 'WOff_Ind', 'WOFF_DATE', 'WriteOff_Amt', 'CLS_STAMP','LN_TPE', "RMNG_TERM",
                        "Receipt_Inf", "HasRepaid"))

# - Isolate affected cases
teste <- subset(test, ExclusionID==0 & Counter==Max_Counter & (HasSettle == 1 | HasWOff == 1 | HasRepaid == 1))

# [DIAGNOSTIC] Distribution of balance at end-of-life? | Early Settlement
describe(teste[HasSettle==1, Balance]); hist(teste[HasSettle==1, Balance])
### RESULTS: Overwhelmingly close to zero as one would expect, but has some large extreme values
# such that the mean is skewed to 3k (median: 0).

# [TREATMENT] For early settlements, redistribute the 
# remaining balance to [Receipt]
teste[HasSettle==1, Receipt_Adj := Receipt_Inf + Balance]

# [DIAGNOSTIC] Distribution of balance at end-of-life? | Write-off
describe(teste[HasWOff==1, Balance]); hist(teste[HasWOff==1, Balance])
### RESULTS: Overwhelmingly close to zero as one would expect, has some minor
# 'extreme' values (max: ~300), with mean of ~2 (median: 0).

# [TREATMENT] For write-offs, redistribute the remaining
# balance to the Write-off Amount
teste[HasWOff==1, WriteOff_Amt_Adj := WriteOff_Amt + Balance]

# [DIAGNOSTIC] Distribution of balance at end-of-life? | Repaid
describe(teste[HasRepaid==1, Balance]); hist(teste[HasRepaid==1, Balance])
### RESULTS: By design, the balance is very small.

# [TREATMENT] For repaid accounts, redistribute the 
# remaining balance to [Receipt]
teste[HasRepaid==1, Receipt_Adj := Receipt_Inf + Balance]

# [DIAGNOSTIC] Prevalence of mismatched early settlement amounts observed
teste[HasSettle==1 & Counter==Max_Counter, Settle_Diff := EarlySettle_Amt - Receipt_Adj]
describe(teste[HasSettle==1, Settle_Diff])
hist(teste[HasSettle==1, Settle_Diff], breaks=2*teste[HasSettle==1, .N]^(1/3))
### RESULTS: Mostly zero-valued, as expected, but with some very large negative extreme 
# values that skews the mean to ~-49 (median: 0).
(diag.real2_3a <- teste[ExclusionID==0 & HasSettle==1 & Counter==Max_Counter & EarlySettle_Amt != Receipt_Adj, .N] / 
  teste[ExclusionID==0 & HasSettle==1 & Counter==Max_Counter, .N] * 100 )
(diag.real2_3b <- mean(teste[HasSettle==1, Settle_Diff], na.rm=T))
(diag.real2_3c <- median(teste[HasSettle==1, Settle_Diff], na.rm=T))
### RESULTS: In ~99%% of early settlements/closures, a mismatch occurred between the observed settlement amount,
#   itself deduced in the SAS-process, and the actual should-be settlement balance [Receipt_Inf] 
#   that would make the balance zero as expected.

# [DIAGNOSTIC] Prevalence of mismatched write-off amounts observed
teste[HasWOff==1 & Counter==Max_Counter, WOff_Diff := WriteOff_Amt - WriteOff_Amt_Adj]
describe(teste[HasWOff==1, WOff_Diff])
hist(teste[HasWOff==1, WOff_Diff], breaks=2*teste[HasWOff==1, .N]^(1/3))
### RESULTS: Mostly zero-valued, as expected, but with some large negative extreme 
# values that skews the mean to ~-1.5 (median: 0).
(diag.real2_3d <- teste[ExclusionID==0 & HasWOff==1 & Counter==Max_Counter & WriteOff_Amt != WriteOff_Amt_Adj, .N] / 
    teste[ExclusionID==0 & HasWOff==1 & Counter==Max_Counter, .N] * 100 )
(diag.real2_3e <- mean(teste[HasWOff==1, WOff_Diff], na.rm=T))
(diag.real2_3f <- median(teste[HasWOff==1, WOff_Diff], na.rm=T))
### RESULTS: In ~2.4% of write-offs, a mismatch occurred between the observed write-off amount
# and the actual should-be write-off amount that would make the balance zero as expected, when adding [Receipt]

# [DIAGNOSTIC] Prevalence of mismatched zero-valued balance observed
teste[HasRepaid==1 & Counter==Max_Counter, Repaid_Diff := Receipt_Inf - Receipt_Adj]
describe(teste[HasRepaid==1, Repaid_Diff])
hist(teste[HasRepaid==1, Repaid_Diff], breaks=2*teste[HasRepaid==1, .N]^(1/3))
### RESULTS: Mostly zero-valued as expected with extreme values capped at the ZAR cut-off as expected
# Mean of -1.4 (median: 0).
(diag.real2_3g <- teste[ExclusionID==0 & HasRepaid==1 & Counter==Max_Counter & Receipt_Inf != Receipt_Adj, .N] / 
    teste[ExclusionID==0 & HasRepaid==1 & Counter==Max_Counter, .N] * 100 )
(diag.real2_3h <- mean(teste[HasRepaid==1, Repaid_Diff], na.rm=T))
(diag.real2_3i <- median(teste[HasRepaid==1, Repaid_Diff], na.rm=T))

# - Conditional Treatment
if (diag.real2_3a > 0 | diag.real2_3d > 0 | diag.real2_3g > 0) {
  
  cat(paste0("DETECTED: Of early settlements / closures, ", round(diag.real2_3a,digits=1), 
             "% of accounts had event amounts that did not agree with calculated 
             event amounts such that [Balance]=0. Mean difference of ", comma(round(diag.real2_3b,digits=2)), 
             " (median: ", comma(round(diag.real2_3c,digits=2)), "). \n\tCorrecting by redistributing non-zero balance to [Receipt_Inf] ..\n"))
  
  cat(paste0("DETECTED: Of repaid accounts, ", round(diag.real2_3g,digits=1), 
             "% of accounts did not have [Balance]=0 at loan-end. Mean difference of ", comma(round(diag.real2_3h,digits=2)), 
             " (median: ", comma(round(diag.real2_3i,digits=2)), "). \n\tCorrecting by redistributing non-zero ",
             " balance to [Receipt_Inf] ..\n"))
  
  cat(paste0("DETECTED: Of write-offs, ", round(diag.real2_3d,digits=1), 
             "% of accounts had event amounts E_w that did not agree with calculated 
             event amounts such that 
                  B(t_w)' = B(t_w-1)*(1+i(t-1)/12) - R(t_w) - E_w  :=  0 where
             B(t) is the observed balance at time t, B(t)' is the calculated balance at time t, 
             t_w is the write-off point, i(t) is the nominal interest rate at time t, 
             R(t) is the inferred receipt at time t, E_w is the write-off event amount.
             Note: R(t) = B(t_w-1)*(1+i(t-1)/12)- B(t) - E_w by definition.
             
             Proposed correction: adjusting the write-off amount by E_w' = ( E_w + B(t_w) ), followed
                by setting B(t_w) = 0.
             
             Mean difference between E_w and E_w' of ZAR ", 
             comma(round(diag.real2_3e,digits=2)), " (median: ", comma(round(diag.real2_3f,digits=2)), 
             "). \n\tCorrecting ..\n"))
  
  # - Fuse treated accounts' event amounts with main dataset
  # To facilitate treatment success measurement later | should be main dataset
  test <- merge(test, teste[,list(LoanID,WriteOff_Amt_Adj, Receipt_Adj)], by=c("LoanID"), all.x=T)
  
  # [TREATMENT] For early settlements, redistribute the 
  # remaining balance to [Receipt] and correct Early Settlement Amount
  test[ExclusionID==0 & Counter==Max_Counter & HasSettle==1 & EarlySettle_Amt != Receipt_Adj, 
       Receipt_Inf := Receipt_Inf + Balance]
  test[ExclusionID==0 & Counter==Max_Counter & HasSettle==1 & EarlySettle_Amt != Receipt_Adj, 
       EarlySettle_Amt := Receipt_Inf]
  
  # [TREATMENT] For repaids, redistribute the remaining balance to [Receipt].
  test[ExclusionID==0 & Counter==Max_Counter & HasRepaid==1 & Receipt_Inf != Receipt_Adj, 
       Receipt_Inf := Receipt_Inf + Balance]
  
  # [TREATMENT] For write-offs, redistribute the remaining
  # balance to the Write-off Amount
  test[ExclusionID==0 & Counter==Max_Counter & HasWOff==1 & WriteOff_Amt != WriteOff_Amt_Adj, 
        WriteOff_Amt := WriteOff_Amt + Balance]
  
  # [TREATMENT] Set the balance explicitly to zero in either terminal case
  test[ExclusionID==0 & Counter==Max_Counter & (HasSettle==1 | HasWOff==1 | HasRepaid==1),
       Balance := 0]
  
  # - Mark the affected accounts with current treatment ID
  LoanIDs <- unique(teste[ExclusionID==0 & (WriteOff_Amt != WriteOff_Amt_Adj | 
                                              EarlySettle_Amt != Receipt_Adj | 
                                              Receipt_Inf != Receipt_Adj), LoanID])
  test[LoanID %in% LoanIDs & Counter==Max_Counter, TreatmentID := paste0(TreatmentID, ";8")]
  
  Lookup <- subset(datCredit_real, LoanID == unique(test[HasWOff == T, LoanID])[1])
  
  # [SANITY CHECK] Are terminal cases behaving as expected?
  (check_advTreat8a <- test[ExclusionID==0 & Counter==Max_Counter & HasSettle==1 &
                             Receipt_Inf != Receipt_Adj, .N] / 
      test[ExclusionID==0 & Counter==Max_Counter & HasSettle==1, .N] *100 == 0) # Should be true
  (check_advTreat8b <- test[ExclusionID==0 & Counter==Max_Counter & HasWOff==1 &
                              WriteOff_Amt != WriteOff_Amt_Adj, .N] / 
      test[ExclusionID==0 & Counter==Max_Counter & HasWOff==1, .N] *100 == 0) # Should be true
  (check_advTreat8c <- test[ExclusionID==0 & Counter==Max_Counter & HasRepaid==1 &
                              Receipt_Inf != Receipt_Adj, .N] / 
      test[ExclusionID==0 & Counter==Max_Counter & HasRepaid==1, .N] *100 == 0) # Should be true
  cat( (check_advTreat8a & check_advTreat8b & check_advTreat8c)  %?% 'SAFE: Terminal event amounts successfully corrected where necessary.\n' %:% 
         'WARNING: Failed to correct all broken terminal event amounts.\n')
  
  # - Cleanup
  suppressWarnings( test[, `:=`(Receipt_Adj = NULL, WriteOff_Amt_Adj = NULL)])
}


# - Cleanup
rm(test, testb, testc, testd, teste)

# ------ END-RESULT: This experiment led to basic cleaning tasks with associated checks. 
#   It also highlighted redundant fields whose information is already captured by other 
#   fields, thereby advocating the safe removal of the former. Moreover, this experiment 
#   led to Advanced Data Treatment 2 in amending the terminal event flags ([WOff_Ind] or 
#   [EarlySettle_Ind]) for those account closures that are unflagged as such.
#   Lastly, this experiment resulted in Advanced Data Treatment 8 in conforming event amounts
#   to zero-valued balances explicitly set at account termination





# ------ Data Experiment 3: Analyzing loan accounts at their starting period | Exclusions applied

# --- 0. Preliminaries
# - Sample data for testing purposes without contaminating the original dataset
test_full <- subset(datCredit_real, ExclusionID==0, 
                    select=c("LoanID", "Counter", "Max_Counter", "Age", "Date_Origination", "Date",
                             "Term", "Principal", "Balance", "Arrears", "ExclusionID", "TreatmentID",
                             "HasWOff", "HasSettle"))

# - Basic feature engineering to facilitate certain analyses
test_full[Principal > 0, Principal_Ratio := Balance / Principal]
test_full[, Term_Ratio := Age / Term]

# - Account-level aggregates & flags to facilitate specific tests
test_full[, HasZeroBalance_Start := max(ifelse(Counter == 1 & Balance == 0, 1, 0)), by=list(LoanID)]
test_full[, HasMissingOpenDate_all := all(ifelse(is.na(Date_Origination), TRUE, FALSE)), by=list(LoanID)]
test_full[, HasMissingOpenDate_any := any(ifelse(is.na(Date_Origination), TRUE, FALSE)), by=list(LoanID)]

# - Sample again at the first record
test <- subset(test_full, Counter==1)


# --- 1. Testing various fields at the first record

# - Distributional analyses on fields of interest
describe(test$Balance); hist(test$Balance, breaks=2*test[,.N]^(1/3))
### RESULTS: Highly righ-skewed distribution, centered at 585k (median: 476k), all non-negative values 
#   with extreme values. Has seemingly extensive zero-values.
describe(test$Principal); hist(test$Principal, breaks=2*test[,.N]^(1/3))
### RESULTS: Highly righ-skewed distribution, centered at 667k (median: 522k), all non-negative values
#   with extreme values. Has seemingly extensive zero-values.
describe(test$Principal_Ratio); hist(test$Principal_Ratio, breaks=2*test[,.N]^(1/3))
### RESULTS: U-shaped distribution, though the vast majority is centered at the upper end as expected,
#   with mean of 85% (median: 99%). Some very few extreme outliers > 100%, but more importantly are
#   those observations at the very low end of the spectrum towards 0%.
describe(test$Term); hist(test$Term, breaks=2*test[,.N]^(1/3))
### Results: Left-skewed distribution, overwhelmingly centered at 240-months as expected.
#   Has some extreme outliers to either side.
describe(test$Term_Ratio); hist(test$Term_Ratio, breaks=2*test[,.N]^(1/3))
# Meaningless histogram due to effects of extreme right-side outliers. Subsample first.
testb <- subset(test, Term_Ratio < 2); hist(testb$Term_Ratio, breaks=2*testb[,.N]^(1/3))
### RESULTS: Highly right-skewed distribution as expected (since truly new accounts should have values close to 0),
#   with a mean of 13% (median: 0%), but very large extreme right-side outliers.
describe(test$Age); hist(test$Age, breaks=2*test[,.N]^(1/3))
### Results: Very dramatically right-skewed distribution as expected, mean of 29.3 (median: 0),
#   with extreme outliers,though only non-negative values. Given the lengthy mortgage registration
#   process, we may need to deduce a cut-off of sorts below which accounts are  considered as new, 
#   and above which loans are old merely observed later in loan life. This will have implications 
#   for the [Age]-field any subsequent survival modelling.


# --- 2. Missingness in [Date_Origination] at first record

# [DIAGNOSTIC] Prevalence of missingness in [Date_Origination]
(diag.real3_1e <- test[HasMissingOpenDate_all==T, .N] / test[,.N] * 100 )
(diag.real3_1f <- test[HasMissingOpenDate_any==T, .N] / test[,.N] * 100 )
### RESULTS: Very minor missingness (~0.01% of accounts) that can perhaps be 
#   fixed by simply assigning the first observed date as the "open date". 
#   Also, missingness seems to be universally spread across the entire credit 
#   history of each affected account. This corroborates suggested treatment 
#   further with no fears of "overwriting" previous values in [Date_Origination].

# - Conditional treatment
if (diag.real3_1e > 0) {
  
  # [TREATMENT] Set first observed date as the "open date"
  test_full[HasMissingOpenDate_all==T, 
            Date_Origination := rollback(Date[1], roll_to_first = T)-days(1), by=list(LoanID)]
  
  # - Recalculate derived fields that depend on this treated field
  test_full[HasMissingOpenDate_all==T, Age := interval(Date_Origination, Date) %/% months(1)]
  
  # [LOOKUP] Test case
  Lookup <- subset(test_full, LoanID == unique(test_full[HasMissingOpenDate_all==T, LoanID])[1])
  
  # - Recalculate affected aggregates 
  test_full[, HasMissingOpenDate_all := all(ifelse(is.na(Date_Origination), TRUE, FALSE)), by=list(LoanID)]
  
  # [SANITY CHECK] Confirm successful treatment
  (check_treat2 <- test_full[Counter==1 & ExclusionID==0 & HasMissingOpenDate_all==T, .N] / 
      test_full[Counter==1 & ExclusionID==0, .N] * 100 == 0) # account-level
  cat( check_treat2 %?% 'SAFE: Missing value treatment on [Date_Origination] succeeded.\n' %:% 
         'WARNING: Missing value treatment on [Date_Origination] failed.\n')
}


# --- 3. Investigating and fixing zero-valued starting balances

# [DIAGNOSTIC] Prevalence of accounts that start with a zero-valued balance?
(diag.real3_1a <- test[Balance == 0, .N] / test[, .N] * 100 )
(diag.real3_1b <- test[Principal == 0, .N] / test[, .N] * 100)
(diag.real3_1c <- test[Balance == 0 & Principal == 0, .N] / test[, .N] * 100)
diag.real3_1c / diag.real3_1b
### RESULTS: ~ 2% of accounts start with a zero-balance, which is seemingly 
#   significant. However, only ~ 0.006% of accounts have a zero-valued starting 
#   principal, and all of these cases coincide with a zero-valued 
#   balance. It may be that the 2% of zero-valued starting balances are simply
#   older accounts whom we have not observed from the start in our sample. 
#   These zero-valued balances are perhaps simply behavioural.

# - Prepare conditional treatment
if (diag.real3_1a > 0) {
  
  # - Find the "true" start for the cases that start with a zero-balance
  test_full[HasZeroBalance_Start==1, Counter_TruStart1 := 
              max(Counter[1], which(Balance > 0)[1]),
            by=list(LoanID)]
  
  # [DIAGNOSTIC] For those cases with zero starting balances, can we find
  # a more appropriate starting point?
  (diag.real3_1d <- test_full[Counter==1 & Counter_TruStart1 != Counter, .N] / 
      test_full[Counter==1, .N] * 100)
  (diag.real3_1d_rec <- test_full[Counter<Counter_TruStart1, .N] / 
      test_full[, .N] * 100 )
  diag.real3_1d / diag.real3_1a # Should be 1
  ### RESULTS: Yes, 100% of cases have a non-zero balance later in life.
  #   We can therefore move the starting point for these cases later in life.
  # - What about the few accounts for which this is apparently not possible (having zero balances throughout) ?
  LoanIDs <- setdiff(test[Balance == 0, LoanID], 
                     test_full[Counter==1 & Counter_TruStart1 != Counter, LoanID])
  Lookup <- subset(datCredit_real, LoanID == LoanIDs[1])
  # [POST HOC] These few negligible accounts fell away after refining both Exclusion 1 and 2 a bit.
  
  # [DIAGNOSTIC] Assessing the potential impact of a possible treatment
  testb <- subset(test_full, Counter == 1 & Counter_TruStart1 != Counter)
  testb[, Degree := Counter_TruStart1 - Counter]
  describe(testb$Degree); hist(testb$Degree, breaks=2*testb[,.N]^(1/3))
  ### RESULTS: Bi-modal right-skewed distribution with non-negative values (>0) and a 
  #   mean of 29.98 (median: 12). Impact seems quite dramatic.
  # [DIAGNOSTIC] Proportion of these cases with zero starting balance for which the fix
  # will return the nth record below the maximum?
  testb[Counter_TruStart1 < Max_Counter, .N] / testb[,.N] * 100
  ### RESULTS: The envisaged treatment will work in the Vast majority of accounts (99.8%)
  # [LOOKUP] The breaking case?
  Lookup <- subset(test_full, LoanID == unique(testb[Counter_TruStart1==Max_Counter, LoanID])[1])
  # 2-record account that can be auto-excluded later (when rerunning Exclusion 1 after Treatment 2)
  # [LOOKUP] An ordinary case to serve as guidance for crafting a treatment
  Lookup <- subset(test_full, LoanID == unique(testb[Counter_TruStart1<Max_Counter, LoanID])[1])
  
  
  # - Towards a treatment using a Lookup-case
  # [TREATMENT] Mark the affected histories as an Exclusion
  test_full[ExclusionID == 0 & HasZeroBalance_Start==1 & Counter<Counter_TruStart1, ExclusionID := 5]
  
  # [TREATMENT] Re-assign the affected counter fields accordingly using the true starting point
  test_full[ExclusionID == 0 & HasZeroBalance_Start==1, Max_Counter:= .N, by=list(LoanID)]
  test_full[ExclusionID == 0 & HasZeroBalance_Start==1, Counter:= 1:.N, by=list(LoanID)]
  
  # - Re-assign [Date_Origination] to be the newly-observed first record, purely
  # to facilitate the recalculation of loan [Age] accordingly. Although this
  # contaminates the previously observed field somewhat, information of this
  # event is recorded in the [TreatmentID] field.
  test_full[ExclusionID == 0 & HasZeroBalance_Start==1, 
         Date_Origination := rollback(Date[1], roll_to_first = T)-days(1), by=list(LoanID)]
  test_full[ExclusionID == 0 & HasZeroBalance_Start==1, 
         Age := interval(Date_Origination, Date) %/% months(1)]
  
  # - Mark the affected accounts with current treatment ID
  LoanIDs <- unique(test_full[ExclusionID==0 & HasZeroBalance_Start==1, LoanID])
  test_full[LoanID %in% LoanIDs, TreatmentID := paste0(TreatmentID, ";1")]
  
  # - Recalculated other derived fields affected by this treated field
  test_full[ExclusionID==0 & HasZeroBalance_Start==1, Term_Ratio := Age / Term]
  
  # - Recalculate affected aggregates
  test_full[ExclusionID==0, HasZeroBalance_Start := max(ifelse(Counter == 1 & Balance == 0, 1, 0)), 
         by=list(LoanID)]

  # [SANITY CHECK] Account-level prevalence of accounts that start with a zero-valued balance?
  (check_advTreat2b <- test_full[ExclusionID==0 & Counter==1 & HasZeroBalance_Start == 1, .N] / 
      test_full[ExclusionID==0 & Counter==1, .N] *100 == 0) # Should be true
  cat( check_advTreat2b %?% 'SAFE: Zero-valued starting balances succesfully treated by adjusting 
       the true starting point.\n' %:% 'WARNING: Treating zero-valued starting balances failed.\n')  
}

# - Cleanup
rm(test_full, test, testb)

# ------ END-RESULT: This experiment led to Missing Data Treatment 2, basic cleaning tasks,
#   Advanced Data Treatment 1 with associated Exclusion 5, and 
#   Data Experiment 7 in trying to isolate truly 'new' accounts.





# ------ Data Experiment 4: Isolating single-record accounts for possible exclusion

# --- 1. Initial experiment

# [DIAGNOSTIC] Prevalence of single-record accounts?
# Put differently, on what proportion of accounts can we feasibly calculate variances?
(diag.real4a <- datCredit_real[Counter == 1 & Max_Counter == 1, .N] / 
    datCredit_real[Counter == 1, .N] * 100 )
### RESULTS: 11.5% of accounts only have a single record, which suggests their removal.

# - Recently disbursed loans may very well have a single record, not due to an intrinsic
#   system error, but simply due to right-censoring. These cases must be kept, 
#   lest we underestimate the number of new loans written recently.
(maxDate <- max(datCredit_real$Date, na.rm=T))

# [DIAGNOSTIC] Investigate age of affected cases via [Date_Origination]
test <- subset(datCredit_real, Counter == 1 & Max_Counter == 1)
describe(test$Date_Origination); barplot(prop.table(table(test$Date_Origination)))
### RESULTS: Fairly old, but distributed across the spectrum, with mean 2011-Apr-28 
#   (median: 2009-Apr-23), though a glut of observations appear towards the far right 
#   (very recent therefore) as well.
testb <- subset(test, Date_Origination >= "2021-05-01")
describe(testb$Date_Origination); barplot(prop.table(table(testb$Date_Origination)))
### RESULTS: It seems the glut starts one month earlier than max(datCredit_real$Date, na.rm=T).
# Upon reflection, this is sensible considering the extraction parameters of the
# underlying process (at the time of writing, the sample period end was set to June-2021).
# This suggests amending the maxDate accordingly

# [DIAGNOSTIC] Prevalence of single-record non-recently disbursed accounts?
# Amend maxDate accordingly to be the first of the last month in which loans were disbursed
(maxDate <- rollback(max(datCredit_real$Date_Origination, na.rm=T), roll_to_first = T))
(diag.real4b <-  datCredit_real[Counter == 1 & Max_Counter == 1 & 
                                  (Date_Origination < maxDate | is.na(Date_Origination)), .N] / 
    datCredit_real[Counter == 1, .N] * 100)
datCredit_real[Counter == 1 & Max_Counter == 1 & (Date_Origination < maxDate | is.na(Date_Origination)), .N] / 
  datCredit_real[Counter == 1 & Max_Counter == 1, .N] * 100
### RESULTS: The previous 11.5% decreased to 11.2%, which is to say ~ 98% of the 11.5% affected cases
# occurred before the latest date, thereby suggesting that 
# these accounts are definitely not recent disbursals and can therefore be excluded

# - cleanup
rm(test, testb)


# --- 2. Re-Isolating single-record accounts | Exclusions/treatments applied
# Advanced Data Treatment 1 affected [Counter] and [Max_Counter] fields,
# which may avail more cases of actual single-record cases for exclusion.

# [DIAGNOSTIC] Prevalence of single-record non-recently disbursed accounts?
# Amend maxDate accordingly to be the first of the last month in which loans were disbursed
(maxDate <- rollback(max(datCredit_real$Date_Origination, na.rm=T), roll_to_first = T))
(diag.real4b_2 <- datCredit_real[ExclusionID==0 & Counter==1 & Max_Counter==1 *
                                   (Date_Origination < maxDate | is.na(Date_Origination)), .N] / 
    datCredit_real[ExclusionID==0 & Counter == 1, .N] * 100 )
### RESULTS: ~ 0.003% of accounts have but a single record, having applied Advanced
#   Data Treatment 1. These cases would have been hidden in the original Data Experiment 4.


# ------ END-RESULT: This experiment led to Exclusion 1 and associated sanity checks,
#   as well as a post-treatment exclusion section





# ------ Data Experiment 5: Finding variances of certain fields over time | Exclusions applied
# Some fields are thought to remain fixed over time per loan, while others 
#   (e.g., Balance) are expected to vary. In order to verify some of these 
#   assumptions about the dataset's structure, one can calculate the sample 
#   variance of certain fields and validate accordingly.
# This experiment also examines partial or complete missingness in certain fields
#   towards crafting an appropriate treatment

# --- 0. Preliminaries
# - Sample data for testing purposes without contaminating the original dataset
selectedFields <- c("LoanID", "Date", "Counter", "Max_Counter" ,"Principal", "InterestRate_Nom",
                    "Balance", "Redraw_Ind", "FurtherLoan_Ind", "Term", "ExclusionID", "TreatmentID")
test <- subset(datCredit_real, ExclusionID == 0, select=selectedFields)

# - calculate  account-level variances of selected fields
test[, Var_Principal := var(Principal, na.rm=T), by=list(LoanID)]
test[, Var_Term := var(Term, na.rm=T), by=list(LoanID)]

# - Create helper flags (account-level)
test[, HasRedraw := max(Redraw_Ind, na.rm=T), by=list(LoanID)]
test[, HasFurtherLoan := max(FurtherLoan_Ind, na.rm=T), by=list(LoanID)]
test[, HasMissingTerm := any(is.na(Term)), by=list(LoanID)]
test[, HasMissingTerm_all := all(is.na(Term)), by=list(LoanID)]
test[, HasMissingPrincipal := any(is.na(Principal)), by=list(LoanID)]
test[, HasMissingPrincipal_all := all(is.na(Principal)), by=list(LoanID)]
test[, HasZeroPrincipal := any(Principal==0), by=list(LoanID)]
test[, HasZeroPrincipal_all := all(Principal==0), by=list(LoanID)]
test[, HasMissingInterest := any(is.na(InterestRate_Nom)), by=list(LoanID)]


# --- 1. Examining [Principal]
# [DIAGNOSTIC] How often is the variance non-zero for [Principal] where the
#   account has neither redraw nor further loan across its history?
(diag.real5_1a <- test[Counter == 1 & Var_Principal > 0 & HasRedraw == 0 & HasFurtherLoan == 0, .N] / 
    test[Counter == 1, .N] * 100 ) # account-level
(diag.real5_1a_rec <- test[Var_Principal > 0 & HasRedraw == 0 & HasFurtherLoan == 0, .N] / 
    test[, .N] * 100 ) # dataset-level
### RESULTS: ~ 25.4% of accounts (or ~ 15.8% of observations) had a varying principal
#   over time that were neither isolated by Redraw-logic nor FurtherLoan-logic. 
#   This is problematic. Enlist [QUANTS]-help in refining SAS data process.

# [DIAGNOSTIC] How prevalent is missingness in parts of [Principal] within the wider dataset
(diag.real5_1b <- test[Counter == 1 & HasMissingPrincipal==T, .N] / 
    test[Counter == 1, .N] * 100 ) # account-level
(diag.real5_1b_rec <- test[ is.na(Principal), .N] / test[, .N] * 100 ) # dataset-level
### RESULTS: No missingness.

# [DIAGNOSTIC] Prevalence of zero-valued Principals throughout the entire account history?
(diag.real5_1d <- test[ExclusionID==0 & Counter==1 & HasZeroPrincipal_all==T, .N] /
    test[ExclusionID==0 & Counter==1, .N] * 100 )
(diag.real5_1d_rec <- test[ExclusionID==0 & HasZeroPrincipal_all==T, .N] /
    test[ExclusionID==0, .N] * 100 )
### RESULTS: ~0.0007% of accounts (~0.0001% of records) have zero-valued Principals 
#   throughout their history. Possible exclusion
# [LOOKUP] Case with zeros in [Principal] throughout history
Lookup <- subset(test, LoanID == unique(test[HasZeroPrincipal_all==T, LoanID])[1])
# [POST HOC] Resolved with Exclusion 3

# [DIAGNOSTIC] How prevalent are zero-valued Principals in parts of [Principal] within the wider dataset
(diag.real5_1c <- test[ExclusionID==0 & Counter==1 & HasZeroPrincipal==T, .N] / 
    test[ExclusionID==0 & Counter==1 , .N] * 100 ) # account-level
(diag.real5_1c_rec <- test[ HasZeroPrincipal==T, .N] / test[, .N] * 100 ) # dataset-level
### RESULTS: ~ 0.09% of accounts (~0.06% of records) are partially affected by zero-valued principals.
# Given the small prevalence, perhaps a simple fix would be to assign the account-level mode of [Principal]
# to the zero-valued parts.
# [LOOKUP] Get an example case on which we'll prepare a treatment
Lookup <- subset(test, LoanID == unique(test[HasZeroPrincipal==T, LoanID])[1])
### RESULTS: There are cases with trailing zero-valued Principals. A simple fix
#   would be populate the last known non-zero value throughout the remainder of the 
#   affected account's history

# - Conditional treatment
if (diag.real5_1c > 0) {
  
  # [LOOKUP] Get an example case on which we'll prepare a treatment
  Lookup <- subset(test, LoanID == unique(test[HasZeroPrincipal==T, LoanID])[1])
  Lookup <- subset(test, LoanID == 400952) # interesting case that have zero-valued principals at both start and end
  
  # - Find the account-level position of the last known non-zero Principal-value
  test[ExclusionID==0 & HasZeroPrincipal==T, 
       Principal_LastNonzeroPos := rev(which(Principal > 0))[1], 
         by=list(LoanID)]
  
  # - Find the corresponding account-level last non-zero Principal-value itself
  test[ExclusionID==0 & HasZeroPrincipal==T, 
       Principal_LastNonzero := Principal[Principal_LastNonzeroPos], 
         by=list(LoanID)]
  
  # [TREATMENT] Assign last known non-zero Principal-value to trailing zero-valued [Principal]-cases
  test[ExclusionID==0 & HasZeroPrincipal==T & Counter>Principal_LastNonzeroPos,
         Principal := Principal_LastNonzero]
  
  # - Find the account-level position of the first known non-zero Principal-value
  test[ExclusionID==0 & HasZeroPrincipal==T, 
       Principal_FirstNonzeroPos := which(Principal > 0)[1], 
         by=list(LoanID)]
  
  # - Find the corresponding account-level first non-zero Principal-value itself
  test[ExclusionID==0 & HasZeroPrincipal==T, 
       Principal_FirstNonzero := Principal[Principal_FirstNonzeroPos], 
         by=list(LoanID)]
  
  # [TREATMENT] Assign first known non-zero Principal-value to starting zero-valued [Principal]-cases
  test[ExclusionID==0 & HasZeroPrincipal==T & Counter<Principal_FirstNonzeroPos,
         Principal := Principal_FirstNonzero]
  
  # [TREATMENT] Assign the account-level mode of [Principal] to those remaining zero-valued cases
  test[ExclusionID==0 & HasZeroPrincipal==T, 
       Principal := ifelse(Principal > 0, Principal,
                           getmode( .SD[Principal > 0, Principal]) # filter out zero-valued cases when getting mode
       ), 
       by=list(LoanID)]
  
  # - Recalculate affected account-level aggregates
  test[, HasZeroPrincipal := any(Principal==0), by=list(LoanID)]
  
  # [SANITY CHECK] Account-level prevalence of accounts that start with a zero-valued balance?
  (check_advTreat3a <- test[ExclusionID==0 & Counter==1 & HasZeroPrincipal==T, .N] / 
      test[ExclusionID==0 & Counter==1, .N] *100 == 0) # Should be true
  cat( check_advTreat3a %?% 'SAFE: Zero-valued [Principal]-values succesfully treated by interleaving the first
       and last known non-zero values accordingly, along with treating the remaining cases with the mode.\n' %:% 
         'WARNING: Treating zero-valued [Principal]-values failed.\n') 
}


# --- 3. Examining [InterestRate_Nom] | Should not be missing

# [DIAGNOSTIC] How prevalent is missingness in parts of [InterestRate_Nom] within the wider dataset
(diag.real5_3b <- test[Counter == 1 & HasMissingInterest==T, .N] / 
    test[Counter == 1, .N] * 100 ) # account-level
(diag.real5_3b_rec <- test[ HasMissingInterest==T, .N] / test[, .N] * 100 ) # dataset-level
### RESULTS: Some missingness in 0.005% of accounts (0.003% of records)

# - Conditional Treatment
if (diag.real5_3b > 0) {
  
  cat(paste0("DETECTED: Partially missing [InterestRate_Nom]-cases in ", round(diag.real5_3b,digits=3), 
             "% of accounts (", round(diag.real5_3b_rec,digits=3), " of records).",
             "\n\tTreating by assigning the first/last-known non-zero account-level [InterestRate_Nom]-value ...\n"))  
  
  # [LOOKUP] Specific case on which we'll prepare a treatment
  Lookup <- subset(test, LoanID == 125088)
  
  # - Mark the affected accounts with current treatment ID
  test[ExclusionID==0 & is.na(InterestRate_Nom), TreatmentID := paste0(TreatmentID, ";3b")]
  
  # [TREATMENT] Back-fill starting missingness with the first non-missing element (if available),
  # then treat remaining partial missingness by imputing with the last-known non-missing element 
  test[ExclusionID==0 & HasMissingInterest==T, InterestRate_Nom := 
         imputeLastKnown(imputeFirstKnown(InterestRate_Nom)), by=list(LoanID)]
  
  # - Recalculate affected account-levle aggregates
  test[, HasMissingInterest := any(is.na(InterestRate_Nom)), by=list(LoanID)]
  
  # [SANITY CHECK] Account-level prevalence of accounts that start with a zero-valued balance?
  check_advTreat3b <- test[ExclusionID==0 & Counter==1 & HasMissingInterest==T, .N] / 
    test[ExclusionID==0 & Counter==1, .N] *100 == 0 # Should be true
  cat( check_advTreat3b %?% 'SAFE: Missing [InterestRate_Nom]-values succesfully treated by 
    interleaving the first and last known non-missing values accordingly.\n' %:% 
      'WARNING: Failed to treat missing [InterestRate_Nom]-values.\n')
  
}


# --- 2. Examining [Term] | Varying when it should be fixed (unless behavioural event occurred)

# [DIAGNOSTIC] How often is the variance non-zero for [Term] where the account 
#   has neither redraw nor further loan across its history?
(diag.real5_2a <- test[Counter == 1 & Var_Term > 0 & HasRedraw == 0 & HasFurtherLoan == 0, .N] / 
    test[Counter == 1, .N] * 100 ) # account-level
(diag.real5_2a_rec <- test[Var_Term > 0 & HasRedraw == 0 & HasFurtherLoan == 0, .N] /
    test[, .N] * 100 ) # dataset-level
### RESULTS: ~ <0.057% of accounts (or ~ 0.038% of observations) had a varying 
#   Term over time that were neither isolated by Redraw-logic nor FurtherLoan-logic. 
#   Seemingly immaterial, perhaps a simple treatment can be created?
# [POST HOC] Issue resolved
# [LOOKUP] Case with varying term
Lookup <- subset(test, LoanID == unique(test[Var_Term>0, LoanID])[1])
### RESULTS: This looks like data errors. Given the immateriality, a simple fix
# would be to simply assign the account-level mode to those observations with varying [Term]-values

# - Conditional Treatment | Varying Term
if (diag.real5_2a > 0) {
  
  cat(paste0("Varying [Term]-values unaccounted for by Redraws or Further Loans found in ",
             round(diag.real5_2a,digits=4), "% of accounts. \n\tTreating by assigning the account-level mode .."))
  
  # - Calculate account-level mode of [Term]
  test[ExclusionID==0 & Var_Term>0 & HasRedraw==0 & HasFurtherLoan==0,
       Term_Mode := getmode(Term), by=list(LoanID)]
  
  # [TREATMENT] Assign the account-level mode to the non-equal parts
  test[ExclusionID==0 & Term!=Term_Mode & HasRedraw==0 & HasFurtherLoan==0,
       Term := Term_Mode[1], by=list(LoanID)]
  
  # - Recalculate affected aggregates
  test[, Var_Term := var(Term, na.rm=T), by=list(LoanID)]
  
  # - Mark affected accounts with current treatment ID
  test[ExclusionID==0 & Term!=Term_Mode & HasRedraw==0 & HasFurtherLoan==0, 
       TreatmentID := paste0(TreatmentID, ";6")]
  
  # [SANITY CHECK] Account-level prevalence of accounts that varying [Term]-values?
  (check_advTreat6 <- test[Counter == 1 & Var_Term > 0 & HasRedraw == 0 & HasFurtherLoan == 0, .N] / 
      test[Counter == 1, .N] * 100 == 0) # Should be true
  cat( check_advTreat6 %?% 'SAFE: Cases with varying [Term]-values successfully treated with the account-level mode.\n' %:% 
         'WARNING: Failed to treat cases with varying [Term]-values.\n') 
}


# -- Treating Missingness in [Term]
# [DIAGNOSTIC] How prevalent is missingness in [Term] within these affected cases??
testb <- subset(test, Counter == 1 & Var_Term > 0 & HasRedraw == 0 & HasFurtherLoan == 0)
testb[HasMissingTerm==T, .N] / testb[,.N] * 100
### RESULTS: 100% of accounts had some degree of missingness in the Term-vector
# Before treating varying terms, let's first clean-up some of the missingness.
# [POST-HOC] Issue resolved

# [DIAGNOSTIC] How prevalent is missingness in parts of [Term] within the wider dataset
(diag.real5_2b <- test[Counter == 1 & HasMissingTerm==T, .N] / 
    test[Counter == 1, .N] * 100 ) # account-level
(diag.real5_2b_rec <- test[ is.na(Term), .N] / test[, .N] * 100 ) # dataset-level
### RESULTS: has missingness in ~50% of accounts (or ~~2.2% of observations)

# [DIAGNOSTIC] How prevalent is missingness in all of [Term] within the wider dataset
(diag.real5_2c <- test[Counter == 1 & HasMissingTerm_all==T, .N] / 
    test[Counter == 1, .N] * 100 ) # account-level
### RESULTS: has complete missingness in ~1.6% of accounts. As such, a simple 
#   treatment that can suffice in the vast majority of cases suffering only from
#   partial missingness would be simply to assign the mode.

# [TREATMENT] Assign mode of account-level Term-vectors to its missing parts, 
#   using custom "getmode()" function defined in script 0
test[HasMissingTerm == T & HasMissingTerm_all == F, 
     Term := ifelse(is.na(Term), getmode(Term), Term),
     by=list(LoanID)]

# [DIAGNOSTIC] Distributional analysis on [Max_Counter] for those with missing [Term]-values
# Possible treatment for those cases where [HasMissingTerm_all] = T can include 
# assigning [Max_Counter].
describe(test[Counter == 1 & HasMissingTerm_all == T, Max_Counter]); 
hist(test[Counter == 1 & HasMissingTerm_all == T, Max_Counter])
### RESULTS: Right-skewed distribution, min of 1, mean of 15.21 (median: 3). 
#   Can be viable. However, this will introduce problems when forecasting a loan 
#   up to completion later on. Perhaps assign the mode Term from those unaffected 
#   cases, since distributional analyses showed that the majority value is 240-months.
#   Sensible for a mortgage portfolio anyways.

# [TREATMENT] For those cases with universal missing terms, use mode imputation
term_mode <- getmode(test[Counter == 1 & HasMissingTerm_all == F, Term])
test[HasMissingTerm_all == T, Term := term_mode]

# - cleanup
rm(test, testb)

# ------ END-RESULT: This experiment led to engineering some basic features/fields to 
#   facilitate checks, as well as Exclusion 3. 
#   Moreover, several related issues were identified whereby the logic 
#   for deducing further loans and/or redraws are deficient.
#   Partially zero-valued [Principal]-cases and missing [InterestRate_Nom] cases were 
#   found and treated in Advanced Data Treatment 3 (both 3a and 3b parts)
#   Missing [Term]-values were also found, leading to Missing Value Treatment 1
#   Varying [Term]-values were also found in a very small minority of accounts. This led to
#   Advanced Data Treatment 6.





# ------ Data Experiment 6: Zero-balance accounts | Exclusions applied

# - Preliminaries
test <- subset(datCredit_real, ExclusionID == 0, select=c("LoanID","Counter", "Balance"))

# - Distributional analysis on first record
describe(test[Counter == 1, Balance]); 
hist(test[Counter == 1, Balance], breaks=2*test[Counter == 1, .N]^(1/3))
### RESULTS: Highly right-skewed distribution, but no negatives, mean of 594k (median: 463k).
# However, seemingly large proportion of zero-values, based on histogram

# [DIAGNOSTIC] What proportion of first records have zero-balances?
test[Counter == 1 & Balance == 0, .N] / test[Counter == 1, .N] * 100
### RESULTS: 3.7% of accounts have zero-balances at the start. 
# Note: This is a bigger problem investigated in Data Experiment 3

# [DIAGNOSTIC] What proportion of account's credit histories have zero balances throughout?
test[ , HasZeroBalances := all(ifelse(Balance <= 0, TRUE, FALSE)), by=list(LoanID)]
(diag.real6a <- test[Counter == 1 & HasZeroBalances == TRUE, .N] / test[Counter == 1, .N] * 100)
### RESULTS: ~ 1.7% of accounts have zero-balances across time. This is a candidate exclusion

# - Cleanup
rm(test)

# ------ END-RESULT: This experiment led to Exclusion 2 and associated sanity checks.





# ------ Data Experiment 7: Identifying truly 'new' loan accounts | Exclusions applied

# --- 0. Preliminaries
# - Sample data for testing purposes without contaminating the original dataset
test_full <- subset(datCredit_real, ExclusionID==0, 
                    select=c("LoanID", "Counter", "Max_Counter", "Age", "Date_Origination", "Date",
                             "Term", "Principal", "Balance", "Arrears", "ExclusionID", "TreatmentID",
                             "HasWOff", "HasSettle"))

# - Basic feature engineering to facilitate analyses
test_full[Principal > 0, Principal_Ratio := Balance / Principal]
test_full[, Term_Ratio := Age / Term]
test_full[, LastRec_Date := max(Date), by=list(LoanID)]

# - Sample again at the first record
test <- subset(test_full, Counter==1 & ExclusionID==0 )

# OBJ: Need to find the "true" starting point for when credit risk measurement should start per account,
#   given the lengthy and highly variable mortgage origination process, further exacerbated by
#   data issues. To facilitate this investigation, use only newly-originated accounts within 
#   the sampling period. Any insight we derive ought to generalize to the broader dataset, e.g.,
#   [Age] starts at 1 for all first records within k=6 months of [Date_Origination]. Need to find k.
testc <- subset(test, Date_Origination >= minDate_observed); 
testc[,.N] / test[,.N] * 100 # Sampling fraction: ~55%


# --- 1. Preliminary first-record analysis: loan ages, balance-to-principal ratios

# [Diagnostic] Prevalence of zero-valued ages for seemingly new accounts that were
# disbursed during the sampling period?
(diag.real7a <- testc[ExclusionID==0 & Counter==1 & Date_Origination >= minDate_observed & Age==0,.N] / 
  testc[ExclusionID==0 & Counter==1 & Date_Origination >= minDate_observed, .N] * 100 )
### RESULTS: Majority of accounts will not need treatment (~93%). However,
# the remaining ones are too numerous to warrant exclusion, without sacrificing data,
# perhaps needlessly so.

# [DIAGNOSTIC] Explore the first record of more recently disbursed accounts
describe(testc$Age); hist(testc$Age, breaks=2*testc[,.N]^(1/3))
# Histogram not meaningful due to extreme outliers. Subsample again just for graphing purposes
testc_samp <- subset(testc, Age <= 10); hist(testc_samp$Age)
### RESULTS: Age distribution is highly right-skewed with some extreme positive outliers. 
#   Mean of 0.4 (median: 0), though a non-zero value only started occuring from 95th percentile.
#   Distributional analysis suggests a low cut-off. Also [Age] should not start at 0, given
#   that month-end balances are observed, which implies at least some time has lapsed (not 0)

# [DIAGNOSTIC] Is [Term_Ratio] negatively correlated with [Principal_Ratio] as one would expect
# for truly new loan accounts? Discard extreme vales for now and old loans.
testc_samp <- subset(testc, Term_Ratio <= 1 & Principal_Ratio <= 1 & !is.na(Principal_Ratio) & 
                       !is.na(Term_Ratio))
cor.test(x=testc_samp$Term_Ratio, y=testc_samp$Principal_Ratio, method="pearson")
ggscatter(testc_samp, x="Term_Ratio", y="Principal_Ratio", cor.method="pearson", cor.coef=T, add="reg.line", conf.int=T)
# cor.test(x=testb$Term_Ratio, y=testb$Principal_Ratio, method="kendall") # won't compute for some reason. R bombs out.
# ggscatter(test, x="Term_Ratio", y="Principal_Ratio", cor.method="kendall", cor.coef=T, add="reg.line", conf.int=T)
### RESULTS: Significant but very minor linear correlation (Pearson, flawed, but computable at least).
#   No clear pattern, other than observing a clump of observations towards the upper-left corner, i.e.,
#   Small values of [Term_Ratio] and big values of [Principal_Ratio]. This is to be the norm.

# - Various lookups to help bootstrap a solution
# [LOOKUP] Interesting case of a very old loan (but still originated within sampling period ) for which we miss an awful lot of observations
Lookup <- subset(test_full, LoanID == unique(testc[Term_Ratio >= 0.25, LoanID])[1] )
# [LOOKUP] Case for which we'll need to solve. Supposedly 'old' loan that is actually brand new and should have [Age]=1
Lookup <- subset(test_full, LoanID == unique(testc[Age >=6 & Age <=42, LoanID])[1] )
# [LOOKUP] Case with lowish balance_ratio
Lookup <- subset(test_full, LoanID == unique(testc[Principal_Ratio >= 0.8 & Principal_Ratio <= 0.9, LoanID])[1] )


# --- 2. Systematic threshold-based analysis: Age-to-term vs Balance-to-Principal

# - Assumption 1: Truly newly-disbursed loans should have an outstanding balance very close to principal
# i.e., Principal_Ratio ~ 100%. As we observe loans later in life, this value should theoretically
# decrease.

# - Assumption 2: Truly newly-disbursed loans should have the vast majority of their remaining lifetime
# still outstanding, i.e., Term_Ratio ~ 0%. As we observe loans later in life, this value should 
# theoretically increase.

# - Define a custom function by which we'll test a threshold (or set thereof) imposed
#   on age-to-term [Term_Ratio] and on balance-to-principal [Principal_Ratio]
#   to isolate truly new accounts from data/system difficulties that presents as 'new' accounts.
#  For each threshold, examine the prevalence of resulting 'new' accounts with 
#   associated starting age and balance-to-principal profiles. 
# output: data.table() object with distributional summaries per threshold-tuple.
testStartingPoint_exp3 <- function(Age_Thresh=0, Balance_Thresh=1, testc, singleTest=T) {
  
  if (singleTest==T) { # provide verbose analysis for a single set of thresholds
    
    # - Indicate whether an account qualifies as a newly disbursed loan, despite the confounding high [Age]-value
    testc[, New_Ind := ifelse(Term_Ratio <= Age_Thresh[1] & Principal_Ratio >= Balance_Thresh[1], 1, 0), by=list(LoanID)]
    
    # [DIAGNOSTIC] Prevalence of 'new' accounts given age threshold
    inner.diag1 <- testc[New_Ind == 1, .N] / testc[,.N]
    cat(paste0("Given thresholds: ", 
               round(inner.diag1*100, digits=2), "% of accounts qualifies as 'new'.\n"))
    
    # [DIAGNOSTIC] Summary statistics of [Age] for those accounts deemed as 'old'
    inner.diag2a <- mean(testc[New_Ind==0, Age], na.rm=T)
    inner.diag2b <- median(testc[New_Ind==0, Age], na.rm=T)
    inner.diag2c <- min(testc[New_Ind==0, Age], na.rm=T)
    inner.diag2d <- max(testc[New_Ind==0, Age], na.rm=T)
    cat(paste0("\tThe 'old' accounts range between ",inner.diag2c, " and ", inner.diag2d, " in [Age],
             with mean ",  round(inner.diag2a,digits=1), " (median: ", inner.diag2b, ").\n"))
    
    # [DIAGNOSTIC] Summary statistics of [Principal_Ratio] for those accounts deemed as 'old'
    inner.diag3a <- mean(testc[New_Ind==0, Principal_Ratio], na.rm=T)
    inner.diag3b <- median(testc[New_Ind==0, Principal_Ratio], na.rm=T)
    inner.diag3c <- min(testc[New_Ind==0, Principal_Ratio], na.rm=T)
    inner.diag3d <- max(testc[New_Ind==0, Principal_Ratio], na.rm=T)
    cat(paste0("\tThe 'old' accounts range between ", round(inner.diag3c*100,digits=2), "% and ", 
               round(inner.diag3d*100,digits=2), "% in [Principal_Ratio],
             with mean ",  round(inner.diag3a*100,digits=2), "% (median: ", 
             round(inner.diag3b*100,digits=2), "%).\n"))
    
    # [DIAGNOSTIC] Summary statistics of [Age] for those accounts deemed as 'new'
    inner.diag4a <- mean(testc[New_Ind==1, Age], na.rm=T)
    inner.diag4b <- median(testc[New_Ind==1, Age], na.rm=T)
    inner.diag4c <- min(testc[New_Ind==1, Age], na.rm=T)
    inner.diag4d <- max(testc[New_Ind==1, Age], na.rm=T)
    cat(paste0("\tThe 'new' accounts range between ",inner.diag4c, " and ", inner.diag4d, " in [Age],
             with mean ",  round(inner.diag4a,digits=1), " (median: ", inner.diag4b, ").\n"))
    
    # [DIAGNOSTIC] Summary statistics of [Principal_Ratio] for those accounts deemed as 'new'
    inner.diag5a <- mean(testc[New_Ind==1, Principal_Ratio], na.rm=T)
    inner.diag5b <- median(testc[New_Ind==1, Principal_Ratio], na.rm=T)
    inner.diag5c <- min(testc[New_Ind==1, Principal_Ratio], na.rm=T)
    inner.diag5d <- max(testc[New_Ind==1, Principal_Ratio], na.rm=T)
    cat(paste0("\tThe 'new' accounts range between ", round(inner.diag5c*100,digits=2), "% and ", 
               round(inner.diag5d*100,digits=2), "% in [Principal_Ratio],
             with mean ",  round(inner.diag5a*100,digits=2), "% (median: ", 
             round(inner.diag5b*100,digits=2), "%).\n"))
    
    return(testc)
    
  } else {
    
    # - Iterate for each given threshold and collect results afterward
    for (i in 1:length(Balance_Thresh)) {
      
      currThresh1 <- Balance_Thresh[i]
      
      # - Indicate whether an account qualifies as a newly disbursed loan, despite the confounding high [Age]-value
      testc[, New_Ind := ifelse(Term_Ratio <= Age_Thresh[1] & Principal_Ratio >= currThresh1, 1, 0), by=list(LoanID)]
      
      # [DIAGNOSTIC] Prevalence of 'new' accounts given age threshold
      inner.diag1 <- testc[New_Ind == 1, .N] / testc[,.N]
      
      # [DIAGNOSTIC] Summary statistics of [Age] for those accounts deemed as 'old'
      inner.diag2a <- mean(testc[New_Ind==0, Age], na.rm=T)
      inner.diag2b <- median(testc[New_Ind==0, Age], na.rm=T)
      inner.diag2c <- min(testc[New_Ind==0, Age], na.rm=T)
      inner.diag2d <- max(testc[New_Ind==0, Age], na.rm=T)
      
      # [DIAGNOSTIC] Summary statistics of [Principal_Ratio] for those accounts deemed as 'old'
      inner.diag3a <- mean(testc[New_Ind==0, Principal_Ratio], na.rm=T)
      inner.diag3b <- median(testc[New_Ind==0, Principal_Ratio], na.rm=T)
      inner.diag3c <- min(testc[New_Ind==0, Principal_Ratio], na.rm=T)
      inner.diag3d <- max(testc[New_Ind==0, Principal_Ratio], na.rm=T)
      
      # [DIAGNOSTIC] Summary statistics of [Age] for those accounts deemed as 'new'
      inner.diag4a <- mean(testc[New_Ind==1, Age], na.rm=T)
      inner.diag4b <- median(testc[New_Ind==1, Age], na.rm=T)
      inner.diag4c <- min(testc[New_Ind==1, Age], na.rm=T)
      inner.diag4d <- max(testc[New_Ind==1, Age], na.rm=T)
      
      # [DIAGNOSTIC] Summary statistics of [Principal_Ratio] for those accounts deemed as 'new'
      inner.diag5a <- mean(testc[New_Ind==1, Principal_Ratio], na.rm=T)
      inner.diag5b <- median(testc[New_Ind==1, Principal_Ratio], na.rm=T)
      inner.diag5c <- min(testc[New_Ind==1, Principal_Ratio], na.rm=T)
      inner.diag5d <- max(testc[New_Ind==1, Principal_Ratio], na.rm=T)
      
      # - Stamp away results
      datIterate <- data.table("Iteration"=i, "Threshold_Age"=Age_Thresh[1], 
                               "Threshold_Balance"= currThresh1, "1_New_Prevalence"=inner.diag1,
                               "2a_Old_Age_Mean"=inner.diag2a, "2b_Old_Age_Median"=inner.diag2b,
                               "3a_Old_Principal_Mean"=inner.diag3a, "3b_Old_Principal_Median"=inner.diag3b,
                               "4a_New_Age_Mean"=inner.diag4a, "4b_New_Age_Median"=inner.diag4b,
                               "5a_New_Principal_Mean"=inner.diag5a, "5b_New_Principal_Median"=inner.diag5b)
      
      if (i==1) datResults <- datIterate else datResults <- rbind(datResults, datIterate)
    }
    
    return(datResults)
  }
  
}

# - Age-to-term threshold vector below which an account is considered as 'new'
(Age_Thresh <- c(0.01, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3))

# - Balance-to-principal threshold vector above which an account is consider as 'new'
( Balance_Thresh <- c(seq(0.99,0.8, by=-0.01), seq(0.775, 0.5, by=-0.025), 0.25) )

# - Execute job across given thresholds
ptm <- proc.time() # for runtime calculations
datResults <- testStartingPoint_exp3(Age_Thresh=Age_Thresh[1], Balance_Thresh=Balance_Thresh, testc, singleTest = F)
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[2], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[3], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[4], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[5], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[6], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[7], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
datResults <- rbind(datResults, testStartingPoint_exp3(Age_Thresh=Age_Thresh[8], Balance_Thresh=Balance_Thresh, testc, singleTest = F ))
proc.time() - ptm

# - Confirm prepared credit data is loaded into memory
if (!exists('datResults')) {
  unpack.ffdf(paste0(genObjPath,"DatExp7-Results"), tempPath)
}

# - additional bit of preparation for graphing multiple series together
datResults_prep <- pivot_longer(datResults, cols=4:ncol(datResults), names_to="Metric", 
                                values_to="Value") %>% 
  drop_na(Value) %>% as.data.table() 
datResults_prep[, Segment := case_when(
  Metric %in% c("2a_Old_Age_Mean", "2b_Old_Age_Median", 
                "3a_Old_Principal_Mean", "3b_Old_Principal_Median") ~ "Old",
  Metric %in% c("4a_New_Age_Mean", "4b_New_Age_Median",
                "5a_New_Principal_Mean", "5b_New_Principal_Median", "1_New_Prevalence") ~ "New"
)]

# - 1. Graphing prevalence across thresholds
chosenCutoff <- 0.925
(g1 <- ggplot(datResults_prep[Metric %in% c("1_New_Prevalence"),], 
              aes(x=Threshold_Balance, y=Value, colour=factor(Threshold_Age), 
                  shape=factor(Threshold_Age))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5, linetype="dashed") + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    labs(y="Prevalence Rate (%) of 'New' Accounts", x="Balance Threshold (%) for 'New' accounts") + 
    theme(legend.position = "bottom") + 
    scale_colour_brewer(palette = "Dark2", name = "Age Threshold (%) for 'New' Accounts") + 
    scale_shape_manual(name="Age Threshold (%) for 'New' Accounts", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=percent) + scale_x_continuous(labels=percent)
)

# - 2. Graphing Age-related distributional summaries across thresholds
(g2 <- ggplot(datResults_prep[Metric %in% c("2a_Old_Age_Mean", "4a_New_Age_Mean"),], 
              aes(x=Threshold_Balance, y=Value, linetype=Metric, shape=factor(Threshold_Age), 
                  colour=factor(Threshold_Age))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5) + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    facet_grid(Segment~., scales="free") + 
    labs(y="Mean Age (months)", x="Balance Threshold (%) for 'New' accounts") + 
    theme(legend.position = "bottom") + scale_linetype_discrete(name="Age segment", labels=c("Old", "New")) + 
    scale_colour_brewer(palette = "Dark2", name = "Age Threshold (%) for 'New' Accounts") + 
    scale_shape_manual(name="Age Threshold (%) for 'New' Accounts", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=comma) + scale_x_continuous(labels=percent)
)

# - 3. Graphing Principal-related distributional summaries thresholds
(g3 <- ggplot(datResults_prep[Metric %in% c("3a_Old_Principal_Mean", "5a_New_Principal_Mean"),], 
              aes(x=Threshold_Balance, y=Value, linetype=Metric, shape=factor(Threshold_Age), 
                  colour=factor(Threshold_Age))) + theme_bw() + 
    geom_point(size=3) + geom_line(size=0.5) + 
    #geom_smooth(method='lm', se=F) + 
    geom_vline(xintercept=chosenCutoff, size=1) + 
    facet_grid(Segment~., scales="free") + 
    labs(y="Mean Balance Ratio (%)", x="Balance Threshold (%) for 'New' accounts") + 
    theme(legend.position = "bottom") + scale_linetype_discrete(name="Age segment", labels=c("Old", "New")) + 
    scale_colour_brewer(palette = "Dark2", name = "Age Threshold (%) for 'New' Accounts") + 
    scale_shape_manual(name="Age Threshold (%) for 'New' Accounts", 
                       values=c(16,17,15,3,7,8,1,13)) + 
    scale_y_continuous(labels=percent) + scale_x_continuous(labels=percent)
)

# - Save graphs
dpi <- 150
ggsave(g1, file=paste0("exp_graphs/Exp3a-PrevalenceNewLoans.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(g2, file=paste0("exp_graphs/Exp3b-MeanAges_NewOld.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(g3, file=paste0("exp_graphs/Exp3c-MeanBalanceRatios_NewOld.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)

# - Store experimental objects | Memory optimisation
pack.ffdf(paste0(genObjPath,"DatExp7-Results"), datResults); rm(datResults, datResults_prep)

### RESULTS: 
# 0. Notation: let balance thresholds (x-axis) be denoted by B. Likewise, let 
#   Age threshold groups (lines) be denoted by A.
# 1. Prevalence curves across A-thresholds in [10%,30%] are remarkably similar,
#   especially towards the preliminary chosen B=92.5% value, yielding 
#   prevalence of 'new' accounts at ~ 87.5%. What this implies is that the majority
#   of first-record cases (originated within the sampling period) will be
#   'new' cases given threshold-values. Sensible.
# 2. Mean_Age curves across B-thresholds are interesting. They seem fairly 
#   stable for 'new'-loans across B-thresholds for all A-groups,
#   except for B>=92.5%, at which point some  mean ages (A in [20%,30%]) dramatically decrease.
#   Moreover, all A-groups are very close to one another in mean age, ranging between
#   [0.05, 0.25]-months, which is sensible for seemingly 'new' loans.
#   Greater volatility exists in 'old'-loans where mean ages decrease fairly consistently
#   over all B-thresholds across all A-groups. But larger/stricter B-values imply greater
#   prevalence of 'older' loans (incorrectly so), which means the age-distribution gets
#   mixed all the more with lower loan ages from actual 'new' accounts, thereby explaining
#   the decreasing trend overall. Of note is the clumping/converging of A-groups in [10%, 20%].
#   Moreover, this decreasing trend seems to hasten after reaching B>=97.5%, which
#   suggests an upper bound of sorts. At B=90%, the range of age-means for A in [10%,17.5%]
#   spans [2.5, ~2.9] vs [2.4,~2.7] for B=92.5%. Since there is no difference in changing B=92.5 
#   to B=90% in the 'new'-loans panel, but a change in the 'old'-loans panel when doing so, 
#   this analysis proposes B=90%.
# 3. Mean_PrincipalRatio (or 'value') increases almost linearly across B-thresholds for all 
#   A-groups for 'new' accounts, spanning 97.25% at B=25% up to ~100% for B=100. Not a lot
#   of differentiation across B-values, and near identical across A-thresholds. Manifests 
#   Assumption 1 very well, especially at B>=90%, yielding values of >=99.6%, i.e.,
#   balance very close to principal. Sensible for seemingly 'new' loans.
#   For 'old'-loans, mean-values increase nonlinearly across B-thresholds for all A-groups,
#   dramatically so for B>=90%. Moreover, some differentiation amongst A-groups, especially
#   for B<90%, though mean values start converging together for B>=90 for most A-groups in 
#   [10%,30%]. Most importantly, at=90%, mean values range in [37%, 38%] for A in [10%, 17.5%],
#   which is sensible for seemingly 'old' loans since one expects a much lower value in
#   balance-to-principal.


# --- 3. Threshold-sensitive analysis on various aspects of loan accounts when first observed

# - Set thresholds accordingly
testc_samp <- testc %>% mutate(New_Ind = ifelse(Age <= 0 | (Term_Ratio <= 0.15 & Principal_Ratio >= 0.9), 1, 0))
#testc_samp <- testStartingPoint_exp3(Age_Thresh=0.15,Balance_Thresh=0.9, testc)
testc_samp_new <- subset(testc_samp, New_Ind == 1)
testc_samp_old <- subset(testc_samp, New_Ind == 0)
### Definition History (excluding Date_Origination >= minDate_observed) for simplcity due to prior 
#   sampling that includes this filter already: 
# Def_1: New_Ind = ifelse(Term_Ratio <= 0.15 & Principal_Ratio >= 0.9, 1, 0)
# Def_2: New_Ind = ifelse(Age == 0 | (Term_Ratio <= 0.15 & Principal_Ratio >= 0.9), 1, 0)
# Def_3: New_Ind = ifelse(Age <= 1 | (Term_Ratio <= 0.15 & Principal_Ratio >= 0.9), 1, 0)

# [DIAGNOSTIC] New loans prevalence?
testc_samp_new[,.N] / testc[,.N]
### RESULTS: Def_1: new-loan prevalence rate: 88%
# Def_2: New-loan prevalence rate: 95.6%
# Def_3: New-loan prevalence rate: 99.5%

# [DIAGNOSTIC] Age-analysis: new vs old.
describe(testc_samp_new$Age)
hist(testc_samp_new$Age, breaks=2*testc_samp_new[,.N]^(1/3))
### RESULTS: Def_1: Vast majority of cases have 0-values in Age, as expected with
#   a mean of 0.14 (median:0) and zeros in all major percentiles. Sensible and
#   manifests Assumption 2, even with some extreme outliers (max: 48).
# Def_2: Mean of 0.13 (median: 0). Distribution has not evolved significantly. Positive result.
# Def_3: Mean of 0.17 (median: 0). Distribution largely unchanged. Positive result.
describe(testc_samp_old$Age)
hist(testc_samp_old$Age, breaks=2*testc_samp_old[,.N]^(1/3))
### RESULTS: Def_1: Vast majority also have 0-valued ages, mean of 2.6 (median: 0),
#   but with remarkably greater extreme values (max: 172) than for new-cases.
#   Note also the prevalence of these 'old' cases is much lower (~12% of first-record cases).
#   However, zero-valued starting ages in 'old' loans are contradictory ..
# Def_2: Largely the same, though mean has shifted to 7.1 (median: 1), which is a positive result.
# Def_3: Distribution changed drastically. New mean of 58 (median: 48), still slightly right-skewed.
#   However, this definition better manifests Assumption 2 in that 'old' loans should be 'older'.
#   A better result.

# [DIAGNOSTIC] Principal Ratio Analysis: new vs old.
describe(testc_samp_new$Principal_Ratio)
hist(testc_samp_new[Principal_Ratio <=1, Principal_Ratio], breaks=2*testc_samp_new[,.N]^(1/3))
### RESULTS: Def_1: Vast majority of values lie in [98%,100%] across all major percentiles as expected, 
#   with a left-skewed distribution and mean of 99.65% (median: 100%). 
#   These facts alone are befitting to 'new' loans. Has some extreme values (min: 90%, max: 350%),
#   though these are presumably data errors.
# Def_2: Regrettably, principal-ratios now span a slightly larger range [96%,100%] for major 
#   percentiles 10% - 95%. Mean: ~96%, median: still 100%. This distributional change is to be expected 
#   though, since those loans with lower principal-ratios are predominantly starting with Age = 0, which
#   should objectively be 'new' loans anyways.
# Def_3: Distribution changed somewhat in left-tail, but this is understandable for the same reason
#   as in Def_2. Mean decreased to 92% (median: 100%), though the distribution now covers [99%,100%] at 
#   only 25% - 95% percentiles (used to be from 10%). Extreme values persist. A mediocre result.
# [LOOKUP] Extreme value lookup. Data error in first [balance]-confirmed.
Lookup <- subset(test_full, LoanID == unique(testc_samp_new[Principal_Ratio > 1.5, LoanID])[1])
describe(testc_samp_old$Principal_Ratio)
hist(testc_samp_old[Principal_Ratio <=1, Principal_Ratio], breaks=2*testc_samp_old[,.N]^(1/3))
### RESULTS: Def_1: Bi-modal almost U-shaped distribution, though the left-side mode is far greater 
#   than the right-side mode. Mean of ~38% (median: 33%), though zero-values are prevalent,
#   which is sensible for seemingly 'old' loans. Threshold-effect visible for values >= 90%.
# Def_2: Bi-modal shape starting to disappear, which manifests Assumption 1 better.
#   Mean of 17% (median: 0.9%).
# Def_3: Significant distributional change, largely left-skewed with mean of 77% (median: 89%).
#   Assumption 1 is more broken with 'old' loans with such a high starting balance-to-principal value.
#   A negative result.

# [DIAGNOSTIC] Analysis on the credit histories of these 'older' but 'new' disbursals
describe(testc_samp_old$Max_Counter); hist(testc_samp_old$Max_Counter)
### RESULTS: Def_1: Mean of ~71 (median: 55), which implies quite lengthy histories.
#   There is merit in somehow preserving these accounts.
# Def_2: Exhibits a more right-skewed distribution, with mean of 64 (median: 54) but retains
#   extreme valeus as expected.
# Def_3: Right-skewed distributional shape still, but mean decreased to 52 (median: 35).
describe(testc_samp_old$LastRec_Date); barplot(prop.table(table(testc_samp_old$LastRec_Date)))
### RESULTS: Def_1: Large left-skewed distributional shape, with mean of 2017-08-03 (median: 2018-11-30),
#   which suggests fairly recent last observation dates. Vast majority, however, are still actively
#   observed.
# Def_2: non-recent values are more frequent than before across the spectrum, though most are still observed
#   fairly recently. The former better supports the idea of truly 'old' accounts.
# Def_3: More similar to Def_1, mean of 2017-09-26 (median: 2018-08-31), greater proportion last observed
#   very recently. A negative result.

# [DIAGNOSTIC] Evidence of zero-valued ages in 'old' accounts?
# An appropriate heuristic may fix these threshold-related failures
(diag.real7b <- testc_samp_old[Age == 0, .N ] / testc[,.N] * 100 )
(diag.real7c <- testc_samp_old[Age == 0, .N ] / testc_samp_old[,.N] * 100 )
### RESULTS: About 7.5% of accounts may have been incorrectly classified as 'old',
#   though this constitutes ~63% of all 'old' accounts. Therefore material.
# [POST-HOC] Issule resolved, having refined the definition to Def_2.
describe(testc_samp_old[Age==0, Principal_Ratio])
hist(testc_samp_old[Age==0 & Principal_Ratio <=1, Principal_Ratio], 
     breaks=2*testc_samp_old[Age==0 & Principal_Ratio <=1,.N]^(1/3))
### RESULTS: Previous modal-shape is still manifested, which confounds results.
#   How can these be 'new' loans yet have a mean balance-to-principal of 50% 
#   (median: 53%), which are characteristic of 'old' loans just observed later 
#   during their loan lives?
# [POST-HOC] Issule resolved, having refined the definition to Def_2.

# [LOOKUP] Confounding new/old case
Lookup <- subset(test_full, LoanID == unique(testc_samp_old[Age==0 & Principal_Ratio < 0.6, LoanID])[2])

### CONCLUSION:
# Def_2 seems to exhibit both assumptions better than both Def_1 and Def_3 and is therefore adopted.

# - Cleanup
rm(test_full, test, testc, testc_samp, testc_samp_new, testc_samp_old)



# --- 4. Towards a treatment: Identifying truly 'new' loan accounts

# - Sample data for treatment purposes without contaminating the original dataset
treat_samp1 <- subset(datCredit_real, ExclusionID==0, 
                    select=c("LoanID", "Counter", "Age", "Date_Origination",
                             "Term", "Principal", "Balance", "ExclusionID", "TreatmentID"))

# - Basic feature engineering & account-level aggregates 
#   to facilitate analyses & treatments.
treat_samp1[Principal > 0, Principal_Ratio := Balance / Principal]
treat_samp1[, Term_Ratio := Age / Term]
treat_samp1[, First_Age := Age[1], by=list(LoanID)]

# [DIAGNOSTIC] Prevalence of accounts disbursed during sampling period?
(diag.real7d <- treat_samp1[ExclusionID==0 & Counter==1 & Date_Origination >= minDate_observed, .N]  / 
  treat_samp1[ExclusionID==0 & Counter==1, .N] * 100 )
### RESULTS: ~54.6% of accounts disbursed during observation as 'new' accounts.

# [DIAGNOSTIC] Prevalence of older ages for seemingly new accounts that were
# disbursed during the sampling period?
(diag.real7e <- treat_samp1[ExclusionID==0 & Counter==1 & Date_Origination >= minDate_observed & Age>0,.N] / 
    treat_samp1[ExclusionID==0 & Counter==1 & Date_Origination >= minDate_observed, .N] * 100 )
(diag.real7e_rec <- treat_samp1[ExclusionID==0 & Date_Origination >= minDate_observed & First_Age>0,.N] / 
    treat_samp1[ExclusionID==0 & Date_Origination >= minDate_observed, .N] * 100 )
### RESULTS: 7% of accounts (7.8% of records) disbursed within sampling period are not 
#   coded correctly as 'new', based on not having a zero-valued Age as expected.

# [TREATMENT] Create an account-level new/old indicator based on derived decision rule
treat_samp2 <- treat_samp1 %>% filter(Counter==1) %>%
  mutate(New_Ind = ifelse(Date_Origination >= minDate_observed & 
                            (Age == 0 | (Term_Ratio <= 0.15 & Principal_Ratio >= 0.9)),
                          1, 0))

# [DIAGNOSTIC] Prevalence rate of new loans: a logic comparison of 'new'
(diag.real7f_a <- treat_samp2[New_Ind == 1, .N] / treat_samp2[,.N] * 100 )
(diag.real7f_b <- treat_samp1[ExclusionID==0 & Counter==1 & Age==0,.N] / 
    treat_samp1[ExclusionID==0 & Counter==1, .N] * 100 )
(diag.real7f_b_rec <- treat_samp1[ExclusionID==0 & First_Age==0,.N] / 
    treat_samp1[ExclusionID==0 , .N] * 100 )
### RESULTS: 52% of accounts are now marked as 'new', whereas
#   50.7% of accounts (42% of records) qualify under crude logic rule, i.e.,
#   having zero-valued ages when first observed.

# [DIAGNOSTIC] Prevalence rate of new loans disbursed during sampling period?
(diag.real7g <- treat_samp2[New_Ind == 1 & Date_Origination >= minDate_observed, .N] / 
    treat_samp2[Date_Origination >= minDate_observed,.N] * 100 )
### RESULTS: ~95.6% of loans disbursed during sampling period are now marked as 'new'.

# - Fuse account-level indicator [New_Ind] back into temporary longitudinal dataset
#   using a left-join on [LoanID]. (Should be main dataset when productionalising this treatment)
treat_samp1 <- merge(treat_samp1, treat_samp2[,list(LoanID, New_Ind)], by=c("LoanID"), all.x=T)

# [SANITY CHECK] Treatment success?
check_advTreat4a <- treat_samp1[ExclusionID==0 & Counter==1 & New_Ind==1, .N] / 
  treat_samp1[ExclusionID==0 & Counter==1, .N] * 100 == diag.real7f_a
cat( check_advTreat4a %?% 'SAFE: Accounts disbursed during sampling period 
  successfully isolated as "new", despite non-zero starting [Age]-values.\n' %:% 
       'WARNING: Failed to isolate accounts disbursed during sampling period.\n') 

# - Recalculate a new [Age] field, adjusted to start at 1 for survival modelling.
treat_samp1[ExclusionID==0, Age_Adj := ifelse(New_Ind==0, Age, 1:.N),
            by=list(LoanID)]

# - Mark affected accounts with current treatment ID
LoanIDs <- subset(treat_samp1, New_Ind==1 & Counter==1 & Age != (Age_Adj), select="LoanID")
treat_samp1[LoanID %in% LoanIDs$LoanID, TreatmentID := paste0(TreatmentID, ";4")]

# [DIAGNOSTIC] Account-level effect of age adjustment?
( diag.real7h_a <- treat_samp1[ExclusionID==0 & Counter==1 & New_Ind==1 & Age != Age_Adj, .N] / 
  treat_samp1[ExclusionID==0 & Counter==1 & New_Ind==1, .N] * 100 )
( diag.real7h_b <- treat_samp1[ExclusionID==0 & Counter==1 & New_Ind==1 & Age != (Age_Adj-1), .N] / 
    treat_samp1[ExclusionID==0 & Counter==1 & New_Ind==1, .N] * 100 )
### RESULTs: ~98% of new accounts have an adjusted [Age]-field that differs
#   from the original. However, this is mostly due to shifting the [Age]-range to start
#   at 1. Only ~2.8% of new accounts truly have a 'fixed' [Age]-field that reflects the
#   fact that they are indeed 'new' accounts, and so should start at [Age]>=1.

# [SANITY CHECK] Are newly-disbursed accounts now correctly aged?
check_advTreat4b <- treat_samp1[ExclusionID==0 & Counter==1 & 
                                  Date_Origination >= minDate_observed & Age_Adj>1,.N] / 
  treat_samp1[ExclusionID==0 & Counter==1 & 
                Date_Origination >= minDate_observed, .N] * 100 < diag.real7e/10 # Should be true
cat( check_advTreat4b %?% 'SAFE: Newly disbursed accounts (during sampling period) are 
    now correctly aged for survival modelling; treatment reduced the original problem by 
    at least 1000%.\n' %:% 'WARNING: Failed to age newly disbursed accounts 
     (during sampling period).\n') 


# ------ END-RESULT: This experiment led to Advanced Data Treatment 4 and associated checks.




# ------ Data Experiment 8: Persistently Zero-valued Instalments | Exclusions/treatments applied
# Should [Receipt]>0 but [Instalment]==0 holds true for significant period,
# then it may suggest measurement error in [Instalment] itself. Such an error
# is further corroborated by the known fact that certain subportfolios came 
# from legacy computer systems that either corrupted the recording of certain 
# fields during migration, or never had the ability to record certain fields.

# --- 0. Preliminaries

# - Sample data for testing purposes without contaminating the original dataset
test <- subset(datCredit_real, ExclusionID==0, 
               select=c("LoanID", "LN_TPE", "Date_Origination", "Date", "Counter", "Max_Counter", "Principal", 
                        "Balance", "Arrears", "Receipt_Inf", "Instalment", "PRE_D7_ACC_STATUS", "STAT_CDE",
                        "HasWOff", "HasSettle",
                        "HasFurtherLoan", "HasRedraw", "FurtherLoan_Ind", "Redraw_Ind",  "Age_Adj", "New_Ind",
                        "ExclusionID", "TreatmentID")) %>% rename(Receipt=Receipt_Inf)


# --- 1. Examining zero-valued instalments/arrears throughout credit histories

# - Create account-level aggregates to facilitate specific tests
test[ExclusionID==0, HasZeroIns_All := ifelse(all(Instalment==0),T,F),by=list(LoanID)]
test[ExclusionID==0, HasZeroArrears_All := ifelse(all(Arrears==0),T,F),by=list(LoanID)]

# [DIAGNOSTIC] Prevalence of account-level zero-valued instalments
(diag.real8_1a <- test[ExclusionID==0 & Counter==1 & HasZeroIns_All==T, .N] / 
  test[ExclusionID==0 & Counter==1, .N] * 100)
(diag.real8_1a_abs <- test[ExclusionID==0 & HasZeroIns_All==T, .N] )
(diag.real8_1a_rec <- diag.real8_1a_abs / test[ExclusionID==0, .N] * 100)
### RESULTS: ~0.1% of accounts (~0.03% of records) have zero-valued instalments.
# Given immateriality, perhaps exclude these cases?

Lookup <- subset(datCredit_real, LoanID == unique(test[Counter==1 & HasZeroIns_All==T, LoanID])[1])

# [DIAGNOSTIC] Prevalence of account-level zero-valued arrears
(diag.real8_1b <- test[ExclusionID==0 & Counter==1 & HasZeroArrears_All==T, .N] / 
    test[ExclusionID==0 & Counter==1, .N] * 100)
(diag.real8_1b_rec <- test[ExclusionID==0 & HasZeroArrears_All==T, .N] /test[ExclusionID==0, .N] * 100)
### RESULTS: ~67% of accounts (~60% of records) have zero-valued arrears throughout
# their histories. However, this is to be expected of most loan portfolios, lest
# the bank becomes bankrupted.

# [DIAGNOSTIC] Prevalence of zero-valued arrears amongst zero-valued instalments?
(diag.real8_1c <- test[ExclusionID==0 & Counter==1 & HasZeroIns_All==T & HasZeroArrears_All==T, .N] /
  test[ExclusionID==0 & Counter==1 & HasZeroIns_All==T, .N] * 100)
### RESULTS: 100% of those with zero-valued instalments also have zero-valued arrears.
# This complicates delinquency measurement itself to the point of near-impossibility.
# Given immateriality, I recommend the exclusion of these cases instead of treatment.
# Perhaps tie this conditionally to diag.real8_1c >= 90% as part of defensive programming.



# --- 2. Examining partially zero-valued instalments throughout credit histories

# - Create lead fields, whilst imputing missing-elements
# with the last known non-missing value, using the custom function "imputeLastKnown()"
# as defined in script 0
test[ExclusionID == 0, Instalment_lead1 := imputeLastKnown(shift(Instalment,type="lead",n=1)), by=list(LoanID)]
test[ExclusionID == 0, Receipt_lead1 := imputeLastKnown(shift(Receipt,type="lead",n=1)), by=list(LoanID)]

# - Create differenced arrears field
test[ExclusionID == 0, Arrears_Diff := c(0,diff(Arrears)), by=list(LoanID)]

# - Indicate persisting zero-valued Instalments across credit histories | absolute basis
# Logic is based on finding consecutively zero-valued instalments, contrasted by non-zero-valued receipts or
# arrears progression - both of which implies that the expected instalment 
# could not have been zero-valued.
test[ExclusionID==0, ZeroIns_Ind := ifelse(Instalment <= 0 & Instalment_lead1 <= 0 & 
                                             ((Receipt > 0 & Receipt_lead1 > 0) | Arrears_Diff > 0)
                                            ,1,0), by=list(LoanID)]

# - Lag the future values of given vector by n periods | n = 1
# This will help in finding definitive 'regimes' of when a non-zero instalment switches to
# a zero-valued instalment.
test[ExclusionID==0, ZeroIns_Ind_lag := shift(ZeroIns_Ind, type="lag", n=1), by=list(LoanID)]

# - Find the starting point of zero-valued instalments | First Regime
# The logic "ZeroIns_Ind_lag==0 & ZeroIns_Ind==1" achieves this idea, particularly by
# taking the index of where this logic is first true. However, there can be multiple 
# such 'regimes'. Therefore, take the index of the first TRUE-element,
# which should correspond to the first such regime.
test[ExclusionID==0, ZeroIns_Start := 
       ifelse(all(!is.na(ZeroIns_Ind)),
              min(coalesce(which((is.na(ZeroIns_Ind_lag) | ZeroIns_Ind_lag==0) & ZeroIns_Ind==1)[1], 0), 
                  Max_Counter), 0),
     by=list(LoanID)]

# - Find the ending point of zero-valued instalments | First Regime
test[ExclusionID==0, ZeroIns_End := 
       ifelse(all(!is.na(ZeroIns_Ind)),
              min(coalesce(which(ZeroIns_Ind_lag==1 & ZeroIns_Ind==0)[1]), 
                  Max_Counter, na.rm=T), 0),
     by=list(LoanID)]

# - Indicate affected accounts for easier traversal
test[ExclusionID==0, HasZeroIns := ifelse(ZeroIns_Start > 0, TRUE, FALSE),
     by=list(LoanID)]

# - Calculate length of period wherein instalments are zero-valued
test[ExclusionID==0 & HasZeroIns==T, ZeroIns_Length := ZeroIns_End - ZeroIns_Start + 1]

# - Extract first two characters from "account status" to examine qualitative delinquency assessment
test <- test %>% mutate(AccStatus_extr = as.numeric(substr(PRE_D7_ACC_STATUS, 1,3))) %>%
  relocate(AccStatus_extr, .after=PRE_D7_ACC_STATUS)

# - Summarize qualitative delinquency status during period with zero-valued instalments
test[ExclusionID==0 & HasZeroIns==T,
     ZeroIns_AccStatus_Mean := mean(AccStatus_extr[ZeroIns_Start[1]:ZeroIns_End[1]], na.rm=T),
     by=list(LoanID)]

# - Find first non-zero instalment | simpler logic
test[ExclusionID==0, Instalment_NonZero_Pos := which(Instalment>0)[1], by=list(LoanID)]

# - Calculate length of period of zero-valued instalment-regimes | simpler logic
test[ExclusionID==0 & Instalment_NonZero_Pos > 1, ZeroIns_Length2 := Instalment_NonZero_Pos-1,
     by=list(LoanID)]


# --- 3. Analysis

# -- Testing cases
# [LOOKUP] Specific curious case found in Exp 1
# Lookup <- subset(test, LoanID == 29975)
# [LOOKUP] Case that has an indicated flag
# Lookup <- subset(test, LoanID == unique(test[HasZeroIns==T, LoanID])[1])
# [LOOKUP] CHL-Case that has an indicated flag
# Lookup <- subset(test, LoanID == unique(test[HasZeroIns==T & LN_TPE=="CHL", LoanID])[1])
# [LOOKUP] Case without an indicated flag
# Lookup <- subset(test, LoanID == unique(test[HasZeroIns==F, LoanID])[1])
# [LOOKUP] Non-qualifying specific case with last instalments being zero
# Lookup <- subset(test, LoanID == 3000000009711)
# [LOOKUP] Qualifying specific case with intermittent zero-valued receipt during zero-valued instalment-regime
# Lookup <- subset(test, LoanID == 18424)
# [LOOKUP] Case that has an indicated flag with zero-valued arrears
# Lookup <- subset(test, LoanID == unique(test[HasZeroIns==T & HasZeroArrears_All==T, LoanID])[1])
# [LOOKUP] Case that has an indicated flag with mean account status of 1
# Lookup <- subset(test, LoanID == unique(test[HasZeroIns==T & ZeroIns_AccStatus_Mean == 1, LoanID])[1])
# [LOOKUP] Case that has an indicated flag with mean account status of > 1
# Lookup <- subset(test, LoanID == unique(test[HasZeroIns==T & ZeroIns_AccStatus_Mean > 1, LoanID])[1])
# [LOOKUP] Case flagged by simpler logic
# Lookup <- subset(test, LoanID == unique(test[Counter==1 & Instalment_NonZero_Pos > 1, LoanID])[1])
# [LOOKUP] Case unflagged by more sophisticated logic but flagged by simpler logic
# Lookup <- subset(test, LoanID == unique(test[Counter==1 & HasZeroIns==F & Instalment_NonZero_Pos > 1, LoanID])[1])

# [DIAGNOSTIC] Prevalence of zero-valued instalments at the account-level?
(diag.real8_2a <- test[ExclusionID==0 & Counter==1 & HasZeroIns==T, .N] / 
  test[ExclusionID==0 & Counter == 1, .N] *100 )
(diag.real8_2a_rec <- sum(test[ExclusionID==0 & HasZeroIns==T & Counter==1, ZeroIns_Length]) / 
    test[ExclusionID==0, .N] *100 )
### RESULTS: 1.3% of accounts (1.2% of records) are affected by persistent zero-valued instalments.
# [POST HOC: Treatment 9] Reduced to ~0.31% of accounts (~0.06% of records).
# [POST HOC: [ZeroIns_Ind] Definition change] Slightly reduced to ~0.29% of accounts (~0.06% of records)
# [POST HOC: Treatment 10] Drastically reduced prevalence to ~0.05% of accounts (~0.001% of records)

# [DIAGNOSTIC] Prevalence of terminal event incidence?
(diag.real8_2b <- test[ExclusionID==0 & Counter==1 & HasZeroIns==T & (HasWOff==T | HasSettle==T), .N] / 
    test[ExclusionID==0 & Counter == 1 & HasZeroIns==T, .N] *100 )
### RESULTS: ~ 99% of the affected accounts have a terminal event, which precludes them from
# exclusion, given their valuable contribution to estimating the very same rare terminal events.
# [POST HOC: Treatment 9] Slightly less prevalent at ~95% of affected accounts having a terminal event
# [POST HOC: [ZeroIns_Ind] Definition change] Prevalence unchanged at ~95%
# [POST HOC: Treatment 10] Prevalence dropped to ~80%

# [DIAGNOSTIC] Prevalence of account-level zero-valued arrears balance
(diag.real8_2c <- test[ExclusionID==0 & Counter==1 & HasZeroIns==T & HasZeroArrears_All==T, .N] / 
    test[ExclusionID==0 & Counter == 1 & HasZeroIns==T, .N] *100 )
### RESULTS: ~64% of accounts have zero-valued arrears throughout their history, which will
# complicate delinquency measurement as well as possible treatments for the zero-valued instalment problem
# [POST HOC: Treatment 9] Significantly increased prevalence at ~89% of accounts having zero-valued arrears.
# [POST HOC: [ZeroIns_Ind] Definition change] Prevalence unchanged at ~89%
# [POST HOC: Treatment 10] Prevalence drastically dropped to ~35%

# [DIAGNOSTIC] Prevalence of subportfolios in those accounts with zero-valued instalments
describe(test[ExclusionID==0 & Counter == 1 & HasZeroIns==T, LN_TPE])
### RESULTS: Majority still CHL (72%), but split is markedly different than at the population-level
# [POST HOC: Treatment 9] Distribution reversed, with majority now being WHL (81%)
# [POST HOC: [ZeroIns_Ind] Definition change] Distribution skewed even more, with WHL numbering 88%
# [POST HOC: Treatment 10] Distribution changed back to majority in CHL (68%)

# [DIAGNOSTIC] Prevalence of zero-valued installments at the last record, perhaps coinciding with a terminal event?
(diag.real8_2d <- test[ExclusionID==0 & Counter==Max_Counter & HasZeroIns==T & Instalment==0, .N] / 
    test[ExclusionID==0 & Counter==Max_Counter & HasZeroIns==T, .N] *100 )
(diag.real8_2d <- test[ExclusionID==0 & Counter==Max_Counter & HasZeroIns==T & Instalment==0 & (HasWOff==T | HasSettle==T), .N] / 
    test[ExclusionID==0 & Counter==Max_Counter & HasZeroIns==T, .N] *100 )
### RESULTS: ~80% of cases had a terminal event coinciding with the zero-valued instalment.
# This suggests that the instalment may have been 'artificially' zeroed upon the account's closure, even 
# though receipts are still very much possible. This can be easily corrected simply by copying the previous instalment
# [POST HOC: Treatment 9] Significantly reduced prevalence of only ~11% of accounts having a zero-valued
# instalment at the end.
# [POST HOC: [ZeroIns_Ind] Definition change] Prevalence unchanged at ~10% of accounts
# [POST HOC: Treatment 10] Drastically changed to ~60% of accounts


# --- 4. Towards a Treatment 9: Treating zero-valued instalments at end-of-loan-life

# - Preliminaries
test[ExclusionID == 0, HasZeroIns_Last := ifelse(Instalment[.N] == 0, T, F), by=list(LoanID)]
test[ExclusionID == 0, HasZeroIns_SecondLast := ifelse(Instalment[.N-1] == 0, T, F), by=list(LoanID)]

# [DIAGNOSTIC] Prevalence of zero-valued instalments at the last record
diag.real8_3a <- test[ExclusionID==0 & Counter==1 & HasZeroIns_Last == T & HasZeroIns_SecondLast == F &
                         (HasWOff==T | HasSettle==T), .N] / test[ExclusionID==0 & Counter==1, .N] * 100 
diag.real8_3a_abs <- test[ExclusionID==0 & Counter==Max_Counter & HasZeroIns_Last == T & 
                            HasZeroIns_SecondLast == F & (HasWOff==T | HasSettle==T), .N]
diag.real8_3a_rec <- diag.real8_3a_abs / test[ExclusionID==0 & Counter==Max_Counter, .N] * 100

# - Conditional treatment
if (diag.real8_3a > 0) {
  
  cat("DETECTED: Zero-valued instalments at end-of-loan-life with terminal events.\n\tPrevalence: ",
      round(diag.real8_3a,digits=2), "% of accounts (", round(diag.real8_3a_rec,digits=2), 
      "% of records).\n")
  
  # [TREATMENT] Copy previous non-zero instalment to last record
  test[ExclusionID==0 & HasZeroIns_Last == T & HasZeroIns_SecondLast == F & (HasWOff==T | HasSettle==T), 
       Instalment := ifelse(Counter<Max_Counter, Instalment, Instalment[.N-1]),
       by=list(LoanID)]
  
  # - Mark affected accounts with current treatment ID
  test[ExclusionID==0 & Counter==Max_Counter & HasZeroIns_Last == T & HasZeroIns_SecondLast == F 
        & (HasWOff==T | HasSettle==T), TreatmentID := paste0(TreatmentID, ";9")]
  
  # - Recalculate affected aggregates
  test[ExclusionID == 0, HasZeroIns_Last := ifelse(Instalment[.N] == 0, T, F), by=list(LoanID)]
  
  # [SANITY CHECK] Are terminal cases behaving as expected?
  check_advTreat9 <- test[ExclusionID==0 & Counter==1 & HasZeroIns_Last == T & 
                            HasZeroIns_SecondLast == F & (HasWOff==T | HasSettle==T), .N] / 
    test[ExclusionID==0 & Counter==1, .N] *100 == 0 # Should be true
  cat( check_advTreat9  %?% 
         'SAFE: Terminated accounts with zero-valued instalments at the end were successfully treated.\n' %:% 
         'WARNING: Failed to treat terminated accounts with zero-valued instalments at their end.\n')
}


# --- 5. Analysis post Treatment 9

# - subsample those accounts marked as zero-valued instalments
test_samp <- subset(test, ExclusionID==0 & HasZeroIns==T & Counter == 1)

# [DIAGNOSTIC] Length of period wherein zero-valued instalments prevail
describe(test_samp$ZeroIns_Length); hist(test_samp$ZeroIns_Length, breaks=2*test_samp[,.N]^(1/3))
### RESULTS: Right-skewed distribution with mean 14.07 (median: 8) with a bi-modal shape, located
# at 0 and 48. t-shaped.

# [DIAGNOSTIC] Distribution of the mean account status during period of zero-valued instalments
describe(test_samp$ZeroIns_AccStatus_Mean); hist(test_samp$ZeroIns_AccStatus_Mean, breaks=2*test_samp[,.N]^(1/3))
### RESULTS: Heavily right-skewed distribution with mean 1.3 (median:1) and expected extreme value of 5 (NPL).
# The results suggest that the overwhelming majority of accounts with partially zero-valued instalments
# were not delinquent (even though the [Arrears]-field is incorrectly zero-valued)
# In turn, a simple fix can simply be to fix [Instalment]=[Receipt] for those periods where the 
# mean account status remained zero.

# [DIAGNOSTIC] Distribution of the starting point of zero-valued instalment regimes
describe(test_samp$ZeroIns_Start); hist(test_samp$ZeroIns_Start, breaks=2*test_samp[,.N]^(1/3))
### RESULT: Heavily right-skewed distribution mean starting point of 11.06 (median: 1).
# This suggests that a sizable part of the affected accounts start-off with no instalment, only
# to be filled later in life, presumably due to fixes effected on the underlying computer system

# [DIAGNOSTIC] Prevalence of accounts with starting zero-valued instalment-regimes
(diag.real8_4a <- test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, .N] / 
  test[ExclusionID==0 & Counter==1, .N] * 100)
(diag.real8_4b <- test[ExclusionID==0 & Counter == 1 & HasZeroIns==T & Instalment_NonZero_Pos > 1, .N] / 
    test[ExclusionID==0 & Counter==1 & HasZeroIns==T, .N] * 100)
(diag.real8_4c <- test[ExclusionID==0 & Counter == 1 & HasZeroIns==F & Instalment_NonZero_Pos > 1, .N] / 
    test[ExclusionID==0 & Counter==1 & Instalment_NonZero_Pos > 1, .N] * 100)
describe(test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, LN_TPE])
describe(test[ExclusionID==0 & Counter == 1 & HasZeroIns==F & Instalment_NonZero_Pos > 1, ZeroIns_Length2])
describe(test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, ZeroIns_Length2])
hist(test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, ZeroIns_Length2])
### RESULTS: ~0.36% of accounts affected, but more interestingly, 83% of those accounts identified
# as having zero-valued instalment-regimes are identified by a much simpler logic, which specifically
# looks at these regimes occurring at the start. Affected cases are also overwhelmingly WHL (~97%),
# which lends credence to the idea of a broken computer system having caused these zero-valued instalment regimes.
# However, 33.6% of accounts identified by the simpler logic were specifically excluded by the more sophisticated
# logic. That said, the overwhelming majority of these unflagged cases (>=90%) had a period length of 1 month.
# This suggests perhaps taking a lead-variant of [Instalment] with this simpler logic in
# fixing the first zero-valued instalment with the first non-zero element.


# --- 6. Towards a treatment 10: Retro-filling zero-valued instalment-regimes at the start

# - Preliminaries
test[ExclusionID==0, Instalment_NonZero_Pos := which(Instalment>0)[1], by=list(LoanID)]
test[ExclusionID==0, Arrears_Diff := c(0,diff(Arrears)), by=list(LoanID)]
test[ExclusionID==0 & Instalment_NonZero_Pos > 1, ZeroIns_Length2 := Instalment_NonZero_Pos-1,
     by=list(LoanID)]

# [DIAGNOSTIC] Prevalence of accounts with starting zero-valued instalment-regimes
diag.real8_4a <- test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, .N] / 
    test[ExclusionID==0 & Counter==1, .N] * 100
diag.real8_4a_abs <- test[ExclusionID==0 & Counter < Instalment_NonZero_Pos & Instalment_NonZero_Pos > 1, .N]
diag.real8_4a_rec <- diag.real8_4a_abs / test[ExclusionID==0,.N] * 100
diag.real8_4d <- mean(test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, ZeroIns_Length2])
diag.real8_4e <- median(test[ExclusionID==0 & Counter == 1 & Instalment_NonZero_Pos > 1, ZeroIns_Length2])

# - Conditional treatment
if (diag.real8_4a > 0) {
  
  cat("DETECTED: Zero-valued instalment-regimes found at the start of credit histories.\n\tPrevalence:",
      round(diag.real8_4a,digits=2), "% of accounts (", round(diag.real8_4a_rec,digits=3), 
      "% of records).\n\tMean regime length of",
      round(diag.real8_4d, digits=1), "months (median: ", round(diag.real8_4e,digits=1), 
      ").\n\tTreating by [Instalment] = [Receipt] where [Receipt] > 0,",
      "\n\totherwise assign the first non-zero [Instalment].\n")
  
  # [TREATMENT] Assign [Receipt] selectively, otherwise assign first non-zero [Instalment]
  # Note: The filter "Counter <= Instalment_NonZero_Pos, Instalment" is necessary to 
  # lookup the first non-zero [Instalment], which "Counter < Instalment_NonZero_Pos, Instalment" will
  # not achieve. Effective contamination is zero since the value at Counter == Instalment_NonZero_Pos
  # is 'overwritten' by itself.
  test[ExclusionID==0 & Instalment_NonZero_Pos > 1 & Counter <= Instalment_NonZero_Pos, Instalment := 
         (Arrears_Diff>=0 & Receipt>0)*Receipt + 
         (Arrears_Diff>=0 & Receipt==0)*Instalment[Instalment_NonZero_Pos[1]]
       ,by=list(LoanID)]
  
  # - Mark affected accounts with current treatment ID
  test[ExclusionID==0 & Instalment_NonZero_Pos > 1 & Counter < Instalment_NonZero_Pos, 
       TreatmentID := paste0(TreatmentID, ";10")]
  
  # - Recalculate affected aggregates
  test[ExclusionID==0, Instalment_NonZero_Pos := which(Instalment>0)[1], by=list(LoanID)]
  
  # [SANITY CHECK] Did the treatment succeed?
  check_advTreat10 <- test[ExclusionID==0 & Counter==1 & Instalment_NonZero_Pos > 1, .N] / 
    test[ExclusionID==0 & Counter==1, .N] *100 == 0 # Should be true
  cat( check_advTreat10  %?% 
         'SAFE: Accounts with zero-valued instalment-regimes at the start were successfully treated.\n' %:% 
         'WARNING: Failed to treat accounts with zero-valued instalment-regimes at the start.\n')
}

# ------ END-RESULT: This experiment led to Exclusion 4, Advanced Data Treatment 9 and Advanced Data Treatment 10.




# ------ Data Experiment 9: Comparison between competing delinquency measurements | Exclusions/treatments applied

# --- 0. Preliminaries

# -- Point to credit data and prepare some graphing-related fields
test <- subset(datCredit_real, ExclusionID == 0, 
                        select=c("LoanID", "Counter", "Max_Counter", "Age_Adj", "Date", "New_Ind", "AccountStatus", "g0_Delinq", "DefaultStatus1", 
                                 "g1_Delinq", "DefaultStatus2", "Event_Type", "Event_Time", "LN_TPE", "Principal", "Term",
                                 "InterestRate_Nom", "Balance", "Arrears", "Instalment", "Receipt", "Redraw_Ind", "Redrawn_Amt",
                                 'FurtherLoan_Ind', "FurtherLoan_Amt", "HasSettle", "EarlySettle_Ind", "EarlySettle_Amt", 
                                 "HasWOff", "WOff_Ind", "WriteOff_Amt")) %>%
  setkey(LoanID, Counter)

# - Create account-level flags to facilitate certain logic tests
test[, HasMismatch := any(DefaultStatus1 != DefaultStatus2), by=list(LoanID)]


# --- 1. Element-wise comparison

# [DIAGNOSTIC] Prevalence of mismatch between default indicators
(diag.real9_1a <- test[Counter==1 & HasMismatch==T, .N] / test[Counter==1, .N] * 100 )
(diag.real9_1a_rec <- test[DefaultStatus1 != DefaultStatus2, .N] / test[, .N] * 100 )
### RESULTS: ~62% of accounts (~49% of records) have mismatches in default status. Very significant

# [DIAGNOSTIC] Prevalence of mismatching indicators over subportfolio?
prop.table(table(test[DefaultStatus1 != DefaultStatus2, LN_TPE])) * 100
### RESULTS: No apparent difference in subportfolio distribution when compared to the unconditional case

# [LOOKUP] Look-up a specific mismatching CHL-case
Lookup <- subset(datCredit_real, LoanID == 3000000057511)
#write_xlsx(Lookup, "Mismatch_Delinquency.xlsx")


# --- 2. Deduce [Receipt] again, primarily based on the SAS-logic

# - Lagged balance and interest fields
test[, Balance_prev := shift(Balance), by=list(LoanID)]
test[Counter==1, Balance_prev := ifelse(New_Ind, Principal[1], Balance[1]), by=list(LoanID)]
test[Counter==1, Balance_prev := Balance[1], by=list(LoanID)]
test[, IntRate_prev := shift(InterestRate_Nom), by=list(LoanID)]
test[Counter==1, IntRate_prev := InterestRate_Nom[1], by=list(LoanID)]
test[, Arrears_prev := shift(Arrears), by=list(LoanID)]
test[Counter==1, Arrears_prev := 0, by=list(LoanID)]

# - Cash flow (Balance difference) between two successive points in time
test[, CashFlow := (1+IntRate_prev/12)*Balance_prev - Balance, by=list(LoanID)]

# - Arrears difference between two successive points in time
#test[, Arrears_Diff := c(0,diff(Arrears)), by=list(LoanID)]
test[, Arrears_Diff := (1+IntRate_prev/12)*Arrears_prev - Arrears, by=list(LoanID)]

# - Receipt inferral function
# NOTE: Since Balance is already treated for terminal events (See Advanced Data Treatment 8), 
# this logic should not correspond exactly to the SAS-logic (which does not treat Balance explicitly .. yet)
test[, Receipt2 := round(ifelse(EarlySettle_Ind==1, EarlySettle_Amt,
                                ifelse(WOff_Ind==1, (1+IntRate_prev/12)*Balance_prev - WriteOff_Amt,
                                       ifelse(Redraw_Ind == 1 | FurtherLoan_Ind == 1, 
                                              pmax(Instalment - pmax(Arrears_Diff,0)),
                                              pmax(CashFlow - pmax(Arrears_Diff,0), 0) ) ) ), digits=2),
     by=list(LoanID)]

# - Calculate difference, bearing in mind the the erstwhile [Receipt]-field will have had some treatments
test[, Receipt_Diff := round(Receipt - Receipt2,digits=2)]

# [DIAGNOSTIC] Prevalence of receipt differences
describe(test$Receipt_Diff); hist(test$Receipt_Diff)
### RESULTS: Vast majority is centered at 0 with mean of -320 (median: 0) and
# some extreme negative/positive outliers. The 10%-90% percentiles are all zero.
# Has some very few missings.

# [DIAGNOSTIC] Extent of differneces
comma(diag.real9_2a <- sum(test$Receipt_Diff, na.rm=T) )
### RESULTS: -115b difference, though this was expected given the improvements in the 
# the inference procedure itself. Also, the [Receipt0] field was not yet treated for 
# end-of-loan-life events, which treatment may be a large portion of the overall differences.

# [LOOKUPS] Missing case | should be zero-case
Lookup <- subset(test, LoanID == unique(test[is.na(Receipt_Inf), LoanID])[1])
# [LOOKUP] Specific case with vast differences in implied delinquency
Lookup <- subset(test, LoanID == 3000000057511)
# [LOOKUP] Case with non-zero [Receipt_Diff]
Lookup <- subset(test, LoanID == unique(test[Receipt_Diff != 0, LoanID])[1])
# [LOOKUP] Case with non-zero [Receipt_Diff] and redraw_event (broken SAS-logic)
Lookup <- subset(test, LoanID == unique(test[Receipt_Diff != 0 & Redraw_Ind==1, LoanID])[1])
# [LOOKUP] Case with non-zero [Receipt_Diff] that settled | should be zero-set
Lookup <- subset(test, LoanID == unique(test[Receipt_Diff != 0 & EarlySettle_Ind==1, LoanID])[1])
# [LOOKUP] Case with non-zero [Receipt_Diff]  that was written-off
Lookup <- subset(test, LoanID == unique(test[Receipt_Diff != 0 & WOff_Ind==1, LoanID])[1])
dat.raw[ACCT_NO==94466 & Counter == Max_Counter, WriteOff_Amt]
# [LOOKUPS] Case with redraw_event
Lookup <- subset(test, LoanID == unique(test[Redraw_Ind==1, LoanID])[1])
# [LOOKUPS] Case with early settlement
Lookup <- subset(test, LoanID == unique(test[HasSettle==1, LoanID])[1])
# [LOOKUPS] Case with write-off
Lookup <- subset(test, LoanID == unique(test[HasWOff==1, LoanID])[1])

### REFINEMENTS / AREAS OF INVESTIGATION (To be discussed with Amina & Quants):
# 1) Use Principal[1] when creating Balance_prev in first record | Discussed
# 2) There appears to be a problem in SAS-code where Receipt != Instalment at 
#     Redraw_Ind==1 (we expect Receipt == Instalment) | Discussed
# 3) Switch order of behavioural events | Discussed 
# 4) There are cases where Redraw_Ind == 1 & EarlySettle_Ind == 1 at the same time. 
#     | Possible treatment to make Redraw_ind (and Furtherloan_ind) = 0 for terminal cases
# 5) Arrears seem to be interest-bearing | Discussed
# 6) Consider extricating [Receipt] from the SAS-process | Agreed
# 7) Ensure [Receipt]-creation is such that resulting [Balance]=0 for terminal events | Advanced Data Treatment 8


# --- 3. Calculate Receipt for a specific account

# [LOOKUP] Specific case that was written-off with multiple default spells
Lookup <- subset(datCredit_real, ExclusionID==0 & LoanID == 3000012424328)

# Balance | Note: Assign first element to Principal if a new loan
Lookup[ExclusionID==0, Balance_prev := shift(Balance), by=list(LoanID)]
Lookup[ExclusionID==0 & Counter==1, Balance_prev := ifelse(New_Ind, Principal[1], Balance[1]), by=list(LoanID)]
# Interest Rate
Lookup[ExclusionID==0, IntRate_prev := shift(InterestRate_Nom), by=list(LoanID)]
Lookup[ExclusionID==0 & Counter==1, IntRate_prev := InterestRate_Nom[1], by=list(LoanID)]
# Arrears
Lookup[ExclusionID==0, Arrears_prev := shift(Arrears), by=list(LoanID)]
Lookup[ExclusionID==0 & Counter==1, Arrears_prev := 0, by=list(LoanID)]

# - Calculate net cash flow (Balance difference) between two successive points in time
Lookup[ExclusionID==0, CashFlow := (1+IntRate_prev/12)*Balance_prev - Balance, by=list(LoanID)]

# - Arrears difference between two successive points in time
Lookup[ExclusionID==0, Arrears_Diff := (1+IntRate_prev/12)*Arrears_prev - Arrears, by=list(LoanID)]
# - Infer the Receipt R(t) given the following logic with balance B(t), interest i(t):
# 1) For terminal events (settlement, write-off) with event amount E at last time t': 
#     R(t') = (1+i(t-1)/12)*B(t-1) - B(t) + E_s, where R(t) >= 0 and E_s is the settlement amount
#     R(t') = (1+i(t-1)/12)*B(t-1) - B(t) - E_w, where R(t) >= 0 and E_w is the write-off amount
#     This should also redistribute any non-zero balance B(t) (which should rightfully be zero; treated in 
#     Advanced Data Treatment 8) to the receipt amount itself.
# 2) For behavioural events (redraw, further loan) with instalment I(t):
#     R(t) = I(t), where R(t) >= 0. 
#     In order for a redraw/further loan to occur, we must reasonably assume the account to be in good standing,
#     otherwise, why 'disburse' even more funds? As such, we can assume that at least one instalment to be paid.
# 3) For all other events, with calculated cash flow C(t) and Arrears balance A(t):
#     R(t) = C(t) + A(t), where R(t) >= 0
#     Any accrued arrears are offset against the presumably positive cash flow,
#     itself signifying a reducing balance, thereby a payment of sorts. Otherwise, if negative, then no payment.
Lookup[ExclusionID==0, Receipt_Inf2 := round(case_when(
  EarlySettle_Ind==1 ~ pmax((1+IntRate_prev/12)*Balance_prev - Balance + EarlySettle_Amt, 0),
  WOff_Ind==1 ~ pmax((1+IntRate_prev/12)*Balance_prev - Balance - WriteOff_Amt, 0),
  (Redraw_Ind == 1 | FurtherLoan_Ind == 1) ~ pmax(Instalment, 0),
  TRUE ~ pmax(CashFlow + Arrears_Diff, 0)
), digits=2), by=list(LoanID)]
Lookup <- Lookup %>% relocate(CashFlow, Arrears_Diff, Receipt_Inf, Receipt_Inf2, Redraw_Ind,FurtherLoan_Ind,slc_past_due_amt,  .after=Instalment)



# ------ END-RESULT: This experiment led to Advanced Data Treatment 7 wherein we infer/recreate
#   the [Receipt_Inf]-field, itself to be compared with a proper one calculated from the transaction-level later on





# ------ Data Experiment 10: Some basic comparative analyses between fields with similar names

# --- 1. [POST_D7_ACC_STATUS] vs [slc_status_final]

# [DIAGNOSTIC Testing for differences between POST_D7_ACC_STATUS and slc_status_final
# Ignore missings; those relate to the shorter sampling period from which the latter 
# was originally extracted
dat.raw[!is.na(slc_status_final) & slc_status_final != "" & !is.na(POST_D7_ACC_STATUS) &
          POST_D7_ACC_STATUS != "", Acc_Diff := POST_D7_ACC_STATUS == slc_status_final]
(diag.real10_1 <- dat.raw[Counter == 1 & Acc_Diff == F, .N] / dat.raw[Counter == 1, .N] * 100)
(diag.real10_1_rec <- dat.raw[Acc_Diff == F, .N] / dat.raw[!is.na(Acc_Diff), .N] * 100)
### RESULTS: ~ 0.05% of accounts (~1.2% of records) have differing non-missing/non-empty values 
# between [POST_D7_ACC_STATUS] and [slc_status_final]
# [LOOKUP] Example case of disagreeing values
#Lookup <- subset(dat.raw, ACCT_NO == unique(dat.raw[Acc_Diff == F, ACCT_NO])[1] & !is.na(Acc_Diff))
Lookup <- subset(dat.raw, ACCT_NO == 3000000001265 & !is.na(Acc_Diff))


# --- 2. [PRE_D7_ACC_STATUS] vs [slc_status_final_pred7]

# [DIAGNOSTIC Testing for differences between PRE_D7_ACC_STATUS and slc_status_final_pred7
# Ignore missings; those relate to the shorter sampling period from which the latter 
# was originally extracted
dat.raw[!is.na(slc_status_final_pred7) & slc_status_final_pred7 != "" & !is.na(PRE_D7_ACC_STATUS) & 
          PRE_D7_ACC_STATUS != "", Acc_Diff_pre := PRE_D7_ACC_STATUS == slc_status_final_pred7]
(diag.real10_2 <- dat.raw[Counter == 1 & Acc_Diff_pre == F, .N] / dat.raw[Counter == 1, .N] * 100)
(diag.real10_2_rec <- dat.raw[Acc_Diff_pre == F, .N] / dat.raw[!is.na(Acc_Diff_pre), .N] * 100)
### RESULTS: Similarly, ~ 0.05% of accounts (~0.47% of records) have differing non-missing/non-empty values
# between [PRE_D7_ACC_STATUS] and [slc_status_final_pred7].

# - Cleanup
dat.raw[, `:=`(Acc_Diff = NULL, Acc_Diff_pre = NULL)]

# ------ END-RESULT: This experiment led to the following conclusions/insights:
# 1) Using [POST_D7_ACC_STATUS] and [PRE_D7_ACC_STATUS] instead of their alternatives, simply 
#   since the former comes from the {SLC_ACCOUNTS_HISTORY}-table, which is deemed more trustworthy and 'finalised'




# ------ Data Experiment 11: Remedial action for extreme values in realised loss rates [LGD]

# --- 1. Preliminary distributional analysis

# - Distributional analysis of realised loss rate (discounted to default point) | Incomplete Portfolio
testa <- subset(datCredit_real, ExclusionID==0 & DefaultStatus1==1 & TimeInDefSpell == 1 & DefSpellResol_Type_Hist == "WOFF")
describe(testa$LossRate_Real); hist(testa[LossRate_Real >=-0.1,LossRate_Real], breaks=2*testa[,.N]^(1/3))
### RESULTS: High concentration at 100% loss seems suspect. As does the extreme negative outliers
# Mean of -15% (median: 24%)

# [LOOKUP] Generic percentile-based lookups
Lookup <- subset(datCredit_real, ExclusionID==0 & LoanID == 
                   unique(testa[LossRate_Real >= 0.99, LoanID])[2])
Lookup <- subset(datCredit_real, ExclusionID==0 & LoanID == 
                   unique(testa[LossRate_Real <=-300, LoanID])[2])
### RESULTS: High losses seem plausible for a few cases investigated.
# However, the extreme minima ('profitable' loss rates) seem suspect

# [LOOKUP] Specific lookup where ReceiptPV = 1070 and WriteOff_Amt = ~30k coincide
# with a very small default [Balance] = ~3. The current definition yields a loss rate of ~-318% which doesn't make sense
Lookup <- subset(datCredit_real, ExclusionID==0 & LoanID==3000003135629,
                 select=c("LoanID","DefSpell_Key", "Counter", "Age_Adj", "Receipt_Inf", "Instalment", "Balance",
                          "ReceiptPV", "LossRate_Real", "AccountStatus", "g0_Delinq", "DelinqState_g0", "g1_Delinq",
                          "DefSpell_Counter", "TimeInDefSpell", "DefSpellResol_Type_Hist"))



# --- 2. Threshold-based EAD-based de-selection "small-balance EADs" as a remedy
# Question: "How many defaults start off with a very small balance? Does their exclusion help?"
# Experiment set up as follows. Raising the ZAR-threshold will increasingly exclude portion of 
# the loss rate sample (bottom panel). As this happens, the resulting mean of loss rates should react positively (uppr panel)
# by increasing to positive-values (a good thing). Simultaneously, the mean of excluded loss rates (middle panel) should largely stay
# negative as a counter-measure. However, I expect this latter mean to start increasing as well, since an increasing
# portion of "good" loss rates are mixed together with the predominant "bad" loss rates that primarily cause the issue
# (of negative loss rates in aggregate)

# - Distributional analysis on EAD at start of default spells ending in write-off
describe(testa$Balance); hist(testa$Balance, breaks = "FD")
### RESULTS: Large proportion of small balances, with the rest of the density being highly right-skewed
# Mean of 573k while median of 497k.

# - Relationship between EAD and loss rate
ggplot(testa, aes(x=(Balance), y=(LossRate_Real))) + geom_point() + theme_bw()
### RESULTS: No apparent relationship. Exp/Log-transforms did not reveal linearity

# - ZAR-denominated threshold vector for EAD below which default spells will be excluded
bal.v <- c(0, 5, 10, 25, seq(50,1500, by=50), seq(1750,7500, by=250), seq(8000,150000,by=1000))
datResults_Exp11  <- data.table(Threshold=bal.v); 

# - Iterative scanning of loss rate distribution, given each EAD-threshold
for (i in 1:length(bal.v)) {
  
  # - Prevalence of defaults with starting balance below threshold
  tempVar1 <- datCredit_real[ExclusionID==0 & DefSpell_Counter==1 & Balance <=bal.v[i], .N] / 
    datCredit_real[ExclusionID==0 & DefSpell_Counter==1, .N]
  
  # - Sample loss rates having excluded defaults with starting balances < threshold
  samp <- subset(datCredit_real, ExclusionID==0 & TimeInDefSpell == 1 & DefSpellResol_Type_Hist == "WOFF" & 
                   Balance > bal.v[i], select=c("LossRate_Real"))
  
  # - Sample loss rates again, exclusively from the excluded subpopulation
  samp2 <- subset(datCredit_real, ExclusionID==0 & TimeInDefSpell == 1 & DefSpellResol_Type_Hist == "WOFF" & 
                    Balance <= bal.v[i], select=c("LossRate_Real"))
  
  # - Distributional summaries
  tempVar2 <- mean(samp$LossRate_Real, na.rm=T)
  tempVar3 <- median(samp$LossRate_Real, na.rm=T)
  tempVar4 <- quantile(samp$LossRate_Real, prob=0.05, na.rm=T)
  tempVar5 <- quantile(samp$LossRate_Real, prob=0.10, na.rm=T)
  tempVar6 <- mean(samp2$LossRate_Real, na.rm=T)
  tempVar7 <- quantile(samp2$LossRate_Real, prob=0.9, na.rm=T)
  tempVar8 <- quantile(samp2$LossRate_Real, prob=0.95, na.rm=T)
  tempVar9 <- median(samp2$LossRate_Real, na.rm=T)
  
  # - Stamp away results for this iteration
  datResults_Exp11[i, Prevalence := tempVar1]
  datResults_Exp11[i, LossRate_Mean := tempVar2]
  datResults_Exp11[i, LossRate_05_perc := tempVar4]
  datResults_Exp11[i, LossRate_10_perc := tempVar5]
  datResults_Exp11[i, LossRate_50_perc := tempVar3]
  datResults_Exp11[i, LossRate_Excl_Mean := tempVar6]
  datResults_Exp11[i, LossRate_Excl_90_perc := tempVar7]
  datResults_Exp11[i, LossRate_Excl_95_perc := tempVar8]
  datResults_Exp11[i, LossRate_Excl_50_perc := tempVar9]
}

# - Store experimental objects | Memory optimisation
pack.ffdf(paste0(genObjPath,"DatExp11-Results"), datResults_Exp11)

# - Confirm prepared credit data is loaded into memory
if (!exists('datResults_Exp11')) {
  unpack.ffdf(paste0(genObjPath,"DatExp11-Results"), tempPath)
}

# - Prepare graphing object
plot.data <- pivot_longer(datResults_Exp11 , cols=c("Prevalence", "LossRate_Mean", "LossRate_05_perc", 
                                            "LossRate_10_perc", "LossRate_50_perc", "LossRate_Excl_Mean",
                                             "LossRate_Excl_95_perc", "LossRate_Excl_50_perc"),
                          names_to = "Type", values_to = "Value") %>% as.data.table()
plot.data[, FacetGroup := 
            case_when(Type %in% c("LossRate_Mean", "LossRate_05_perc", 
                                  "LossRate_10_perc", "LossRate_50_perc") ~ "a_LossRates",
                      Type %in% c("LossRate_Excl_Mean", "LossRate_Excl_90_perc", "LossRate_Excl_95_perc", "LossRate_Excl_50_perc") ~ "b_LossRates_Excl",
                      TRUE ~ "c_Prevalence")]

# - Graphing options
label.v <- c("LossRate_05_perc"="5th Percentile [Loss Rates]", "LossRate_10_perc" = "10th Percentile [Loss Rates]",
             "LossRate_50_perc"="Median [Loss Rates]", "LossRate_Mean"="Mean [Loss Rates]",
             "Prevalence"="Prevalence [EAD <= Threshold]", "LossRate_Excl_Mean"="Mean of excluded [Loss Rates]",
             "LossRate_Excl_90_perc"="90th Percenttile of excluded [Loss Rates]",
             "LossRate_Excl_95_perc"="95th Percenttile of excluded [Loss Rates]",
             "LossRate_Excl_50_perc"="Median of excluded [Loss Rates]")
facet_names <- c("a_LossRates"="Summary Statistics: included loss rates", "b_LossRates_Excl"="Summary Statistics: excluded loss rates",
                 "c_Prevalence"="Sample rate (of excluded loss rates)")

# - Graph results
thresh_EADLoss <- 6500 # Preliminary threshold choice
(g <- ggplot(plot.data, aes(x=Threshold, y=Value)) + theme_bw() + 
    geom_line(aes(colour=Type, linetype=Type)) + geom_point(aes(colour=Type)) + 
    geom_vline(xintercept=thresh_EADLoss, size=1) + 
    labs(x="Threshold for EAD (ZAR)", y="Value (%)") + 
    facet_grid(FacetGroup ~ ., scales="free", labeller=labeller(FacetGroup=as_labeller(facet_names))) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
    scale_colour_brewer(name="", palette="Dark2", labels=label.v, guide=guide_legend(ncol=3, byrow=T)) + 
    scale_linetype_discrete(name="", labels=label.v) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_continuous(breaks=pretty_breaks(), label=comma))

# - Store graphical results
dpi <- 140
ggsave(g, file=paste0("exp_graphs/Exp11a-LossRate_SampleDynamics.png"),width=1200/dpi, height=1500/dpi,dpi=dpi)

### RESULTS: Cut-off of 6500 seems adequate. However, the 5% and 10% percentiles are still negative,
# which suggests that this threshold-based remedy aren't wholly effective.



# --- 3. Impact analysis of chosen EAD-threshold

# - Sample loss rates having excluded defaults with starting balances < threshold
samp <- subset(datCredit_real, ExclusionID==0 & TimeInDefSpell == 1 & DefSpellResol_Type_Hist == "WOFF" & 
                 Balance > thresh_EADLoss, select=c("LoanID", "LossRate_Real"))

# - Distributional analysis on loss rates, given remedy
describe(samp$LossRate_Real); hist(samp[LossRate_Real >= 0, LossRate_Real], breaks = "FD")
### RESULTS: Although the extreme minima seems to have dampened, there are still significant negatives (up to 10% percentile)
# mean of 21% (median: 22%)

# - Distributional analysis on loss rates, given remedy + exclude negative rates
describe(samp[LossRate_Real >= 0, LossRate_Real]); hist(samp[LossRate_Real >= 0, LossRate_Real], breaks = "FD")
### RESULTS: mean of 39% (median of 33%), significant shifts in central tendencies

# - Sample loss rates again having excluded defaults loss rates < 0
samp2 <- subset(datCredit_real, ExclusionID==0 & TimeInDefSpell == 1 & DefSpellResol_Type_Hist == "WOFF" & 
                  LossRate_Real >= 0, select=c("LoanID", "LossRate_Real"))

# - Distributional analysis on loss rates, given simple exclusion of negative rates
describe(samp2$LossRate_Real); hist(samp2[LossRate_Real >= 0, LossRate_Real], breaks = "FD")
### RESULTS: Similar shape, but significant differences at particular densities
# mean of 44% (median: 36%)

# [LOOKUP] Offending case with a negative loss rate, despite remedy
Lookup <- subset(datCredit_real, ExclusionID==0 & LoanID == unique(samp[LossRate_Real <= -10, LoanID])[1],
                 select=c("LoanID","DefSpell_Key", "Counter", "Age_Adj", "Receipt_Inf", "Instalment", "Balance",
                          "ReceiptPV", "LossRate_Real", "AccountStatus", "g0_Delinq", "DelinqState_g0", "g1_Delinq",
                          "DefSpell_Counter", "TimeInDefSpell", "DefSpellResol_Type_Hist"))

# - Cleanup
rm(plot.data, samp, samp2, tempVar1, tempVar2, tempVar3, tempVar4, tempVar5, tempVar6,tempVar7, tempVar8, tempVar9,
   testa, Lookup, g, bal.v, datResults_Exp11)


# ------ END-RESULT: This experiment led to the following conclusions/insights:
# 1) Experimental setup did not yield a perfect remedy; the use of ZAR-based thresholds by which
#     certain "small/immaterial" default spells (write-offs) did not solve the primary problem
#     of a loss distribution with extreme negative rates.
# 2) The remedy seemingly introduces distributional contamination, relative to the truncated (but otherwise intact)
#     baseline distribution. 
# 3) Remedy rejected as unsafe, likely introducing severe estimation bias when modelling loss severity (LGW).









  