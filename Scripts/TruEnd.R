# ================================== TruEnd-functions ===================================
# Function definitions used in the TruEnd-procedure itself
# ---------------------------------------------------------------------------------------
# SCRIPT AUTHOR(S): Dr Arno Botha
# VERSION: 1.0 (May-2023)
# DESCRIPTION: 
# See Botha, A., Verster, T., Bester, R. (2023). The TruEnd-procedure: Treating trailing
#   zero-valued balances in credit data (forthcoming). Overleaf
#   https://www.overleaf.com/read/dptmgxnvqjjq
# =======================================================================================



# --- Custom function for isolating a TZB-regime for a given vector [vecGive] and threshold [thres], 
# which should be of minimum length [minLenght], while returning a number of specifiable items, 
# based on return type [retType]: 
# 1) "t_z": starting point of TZB-regime, if found
# 2) "t_z-1": true end point, considering the TZB-regime (if found), otherwise the observed account end
# 3) "tzb_len": length of the TZB-regime, if found
# [nonTZB_Val]: optional parameter that is assigned to the TZB-start t_z for non-TZB accounts.
TruEnd_inner <- function(vecGive, thres=0, retType, minLength=2, nonTZB_Val=-1) {
  # - testing conditions
  # vecGive <- matControl[,i]
  # retType <- "epiNum"
  
  # - Find positions in reversed vector where threshold is triggered
  vec.found <- which(rev(vecGive) <= thres)
  
  # - Find complementing position in reversed vector where threshold is first NOT triggered
  lastNotTrig <- which(rev(vecGive) > thres)[1]
  
  # - Logic check: Only execute if any regimes are isolated and if the last of these regimes 
  # exceed a given minimum length, bearing in mind that the "true end" itself would be <= [thres] 
  # and is not technically part of the TZB-regime (hence minus 1)
  if(length(vec.found) > 0 & !is.na(lastNotTrig) &
     (lastNotTrig-vec.found[1]-1) >= minLength) {logiCheck <- T} else {logiCheck <- F}
  
  # - Find start of TZB-regime, i.e., the t_z point
  # Note: 1st parenthesis (A) is the last position where control variable > [thres]
  # Therefore, the "true end" would be (A) + 1, while the TZB-regime starts at (A) + 2
  if (logiCheck) {t_z <- (length(vecGive) - lastNotTrig+1) + 2} else {t_z <- nonTZB_Val}
  # Possible return point: starting point (as index in given vector) of TZB-regime (if found)
  if (retType == "t_z") {return(t_z)}
  
  # - Find "true end" point (as index in given vector) of TZB-regime (if found)
  if (logiCheck) {tru_end <- t_z - 1} else {tru_end <- length(vecGive)}
  # Possible return point
  if (retType == "t_z-1") {return(tru_end)}
  
  # - Length of the isolated TZB-regime, if found
  if (logiCheck) {tzb_len <- length(vecGive)-t_z+1} else {tzb_len <- 0}
  # Possible return point
  if (retType == "tzb_len") {return(tzb_len)}
  
  # Generic return point
  return(list(t_z=t_z, True_End=tru_end, TZB_Length=tzb_len))
}



# --- Multithreadable function for testing the impact of a specified threshold on a given loan portfolio
# Inputs: [matControl]: data matrix to scan against threshold; [thres]: a particular threshold-value;
# [controlVar]: text describing the control variable being evaluated against the given [thres]
# [controlVar2]: text describing the secondary/optional control variable
# [thres2]: threshold corresponding to the secondary control variable
# [controlVar2VoidVal]: threshold value for voiding the secondary control variable (as if its influence is ignored)
# [matBalance]: balance matrix from which M1 and M2 measures are calculated in evaluating thresholds
# [vecMaturity]: vector of observed maturities (total loan ages) up to the 'false'/untreated account termination time points
# [tau]; length of preceding non-TZB period, [it]: current threshold iteration within broader setup;
# [numThresh]: number of thresholds being tested within broader setup.
# [minLength]: passable argument denoting the minimum length of a TZB-regime, by definition
# [reportFlag]: Boolean-value indicating whether to log activities or not
# [logName]: "filename for log"
# [f.obj]: Given objective function that accepts arguments in [args]
# [args]: any external arguments required by [f.obj]; internally-calculated arguments include M1.v, M2.v, t_z.v, true_end.v, tzb_len.v, TruBal.v
TruEnd_outer <- function(matControl, thres=0, controlVar, matBalance, vecMaturity, tau=6, it=1, numThres=1, minLength=1, 
                         matControl2, thres2=0, controlVar2=NA, controlVar2VoidVal=0, reportFlag=F,logName="General",f.obj, args) {
  
  if (reportFlag & (is.na(controlVar2) | thres2==controlVar2VoidVal)) {
    cat(paste0("\n 1)[", it, " of ", numThres, "] Evaluating control variable [", controlVar, 
               "] against threshold ", thres, " .."), file=paste0("assesslog_", logName,".txt"), append=T)
  } else if (reportFlag) {
    cat(paste0("\n 1)[", it, " of ", numThres, "] Evaluating control variable [", controlVar, 
               "] (primary) against threshold ", thres, ", and also secondary control variable [",
               controlVar2, "] against threshold ", thres2, " .."), file=paste0("assesslog_", logName,".txt"), append=T)
  }
  
  # - Abstract the number of accounts from data
  nAcc <- ncol(matControl)
  
  # - Find preliminary t_z points as the start of a TZB-regime, if found
  t_z.v <- sapply(1:nAcc, function(i){
    # testing: i <- 1
    if (is.na(controlVar2) | thres2==controlVar2VoidVal) { # only primary control variable
      TruEnd_inner(vecGive=matControl[,i], thres=thres, retType="t_z", minLength=minLength)
    } else { # 2 control variables given, use minimum to coalesce signals
      min(TruEnd_inner(vecGive=matControl[,i], thres=thres, retType="t_z", minLength=minLength),
          TruEnd_inner(vecGive=matControl2[,i], thres=thres2, retType="t_z", minLength=minLength))
    }
  })
  
  # - Find preliminary "true end" points
  true_end.v <- sapply(1:nAcc, function(i){
    if (is.na(controlVar2) | thres2==controlVar2VoidVal) {  # only primary control variable
      TruEnd_inner(vecGive=matControl[,i], thres=thres, retType="t_z-1", minLength=minLength)  
    } else { # 2 control variables given, use minimum to coalesce signals
      min(TruEnd_inner(vecGive=matControl[,i], thres=thres, retType="t_z-1", minLength=minLength),
          TruEnd_inner(vecGive=matControl2[,i], thres=thres2, retType="t_z-1", minLength=minLength))
    }
  })
  
  # - Calculate length of isolated TZB-regimes, if found
  tzb_len.v <- sapply(1:nAcc, function(i){
    if (is.na(controlVar2) | thres2==controlVar2VoidVal) {  # only primary control variable
      TruEnd_inner(vecGive=matControl[,i], thres=thres, retType="tzb_len", minLength=minLength)
    } else { # 2 control variables given, use maximum to coalesce signals
      max(TruEnd_inner(vecGive=matControl[,i], thres=thres, retType="tzb_len", minLength=minLength),
          TruEnd_inner(vecGive=matControl2[,i], thres=thres2, retType="tzb_len", minLength=minLength))
    }
  })
  
  # - Find balance at the "true end" points, if found
  TruBal.v <- sapply(1:nAcc, function(i, t){
    return(matBalance[t[i],i])}, t=true_end.v)
  
  # - Measure 1: mean TZB-balance, where found
  M1.v <- sapply(1:nAcc, function(i, t_z, m){
    if (t_z[i] > 0) {
      # Isolate entire TZB-regime and calculate mean
      meanBal <- mean(matBalance[t_z[i]:m[i],i],na.rm=T)
    } else {meanBal <- NA}
    return(meanBal)
  }, t_z=t_z.v, m=vecMaturity)
  
  # - Measure 2: non-TZB mean balance
  M2.v <- sapply(1:nAcc, function(i, tru, tau){
    # Isolate the [tau] values in the vector that precede the TZB-regime (if found), and calculate mean
    if (tru[i]-tau > 0) { # sufficient history exists
      mean(matBalance[(tru[i]-tau):tru[i],i], na.rm=T)
    } else { # insufficient history exists, therefore start from element 1
      mean(matBalance[1:tru[i],i], na.rm=T)
    }
  }, tru=true_end.v, tau=tau)
  

  # --- 3. Evaluate objective function, given measures M1 and M2
  # Calculate gien objective function
  obj.v <- do.call(f.obj, args=args)
  
  
  # --- 4. Concatenate results
  datResults.interim <- data.table(ResultSet = logName, Control = controlVar, Threshold = thres, Control2 = controlVar2, Threshold2 = thres2,
                                   FalseEnd_mean = mean(vecMaturity, na.rm=T), FalseEnd_sd = sd(vecMaturity, na.rm=T),
                                   TruEnd_mean = mean(true_end.v, na.rm=T), TruEnd_sd = sd(true_end.v, na.rm=T),
                                   TZB_Length_mean = mean(tzb_len.v, na.rm=T), TZB_Length_sd = sd(tzb_len.v, na.rm=T),
                                   TruBal_mean = mean(TruBal.v, na.rm=T), TruBal_sd = sd(TruBal.v, na.rm=T),
                                   M1_mean = mean(M1.v, na.rm=T), M1_sd = sd(M1.v, na.rm=T),
                                   M2_mean = mean(M2.v, na.rm=T), M2_sd = sd(M2.v, na.rm=T),
                                   TZB_prevalence = sum(t_z.v>=0)/nAcc, Objective = obj.v[1]
  )
  
  return (datResults.interim)
  
}




