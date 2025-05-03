# ============================== CUSTOM FUNCTIONS ==============================
# Defining custom functions used across various projects
# ------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB), Marcel Muller (MM)
# VERSION: 2.0 (Apr-2025)

# DESCRIPTION:
# This script defines various functions that are used elsewhere in this project
# or, indeed, used across other projects. Functions are grouped thematically.
# ==============================================================================




# -------- Ternary functions
# from https://stackoverflow.com/questions/8790143/does-the-ternary-operator-exist-in-r
`%?%` <- function(x, y) list(x = x, y = y)
`%:%` <- function(xy, z) if(xy$x) xy$y else z




# -------- Utility functions

# -- Mode function (R doesn't have a built-int one)
getmode <- function(v) {
  uniqv <- unique(v);
  # discard any missingness
  uniqv <- uniqv[complete.cases(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# -- Memory function using 'gdata' package
getMemUsage <- function(limit=1000){
  require(gdata); require(scales)
  # - Get list of significant object sizes occupied in memory, order ascendingly
  totUsage <- ll()
  memSize <- subset(totUsage, KB >= limit)
  memSize$MB <- memSize$KB/1000
  gc(verbose=F)
  cat("Total memory used: ", comma(sum(totUsage$KB)/1000), "MB\n")
  cat("Big objects size: ", comma(sum(memSize$MB)), "MB\n\n")
  return(  memSize[order(memSize$KB), c(1,3)])
}




# -------- Cleaning functions (Missing/extreme value treatments)

# -- Custom function that curates a main vector [x] to equal the previous/most-recent non-
# missing element in a given vector
imputeLastKnown <- function (x) {
  # -- Testing purposes
  # x <- Lookup$ZeroBal_Remain_Ind; x_lead <- Lookup$ZeroBal_Remain_Ind_lead
  # x <- c(0,0,0,1,1,1,0,1)
  # x <- c(0,0,0,1,1,1,0,NA)
  # x <- c(0,0,0,1,1,1,1,NA)
  # x <- c(0,0,0,1,NA,1,0,NA)
  # x <- c(0,NA)
  
  firstOne <- which(is.na(x))[1]
  if (!is.na(firstOne) & firstOne > 1) {
    x[firstOne] <- x[firstOne-1]
    # call function recursively to fix earlier missing-cases
    return( imputeLastKnown(x))
  } else { # no missing value found, return original vector
    return(x)
  }
}


# -- Custom function that curates a main vector [x] where x[1] is missing.
# This is achieved by finding the first non-missing element and back-filling that value
imputeFirstKnown <- function(x) {
  # -- Testing purposes
  # x <- c(NA, NA, 2,3,4)
  firstOne <- which(!is.na(x))[1]
  if (!is.na(firstOne) & firstOne > 1) {
    x[1:(firstOne-1)] <- x[firstOne]
    return(x)
  } else { # no non-missing value found, return original vector
    return(x)
  }
}


# -- Custom Function by which to adjust for inflation
# Assumes a monthly macroeconomic dataset [macro_data_hist] to exist with [Date_T] and [Inflation] fields
adjInflation <- function(g_start, g_stop) {
  compFact <- macro_data_hist[Date_T >= g_start & Date_T <= g_stop, list(Factor = prod(1 + (Inflation/100)/12))]
  return(compFact)
}


# - Adjusting for inflation (Robust version that accepts a macroeconomic dataset)
# Generating an inflation factor for a given series of yearly inflation growth rates
# Input:  [datMacro]: The dataset containing the yearly inflation growth rate
#         [time]: Name of the time/date variable in [datMacro]
#         [g_start]:  The starting date for the series of inflation growth rates
#         [g_stop]:   The ending date for the series of inflation growth rates
# Output: A factor indicating the cumulative inflation over the period starting at [g_start] and ending [g_stop]
# --- Define custom function for computing inflation/deflation factors
adjInflation_MV <- function(datMacro, time, Inflation_Growth, g_start, g_stop) {
  # datMacro=datMV; time="Date"; g_start<-date("2015-02-28"); g_stop<-date("2022-12-31"); Inflation_Growth<-"M_Inflation_Growth"
  compFact <- as.numeric(datMacro[get(time) >= g_start & get(time) <= g_stop, list(Factor = prod(1 + (get(Inflation_Growth))/12))])
  return(compFact)
  # rm(datMacro, time, g_start, g_stop, Inflation_Growth); gc()
}
# - Unit test
# if (!exists('datMV')) unpack.ffdf(paste0(genPath,"datMV"), tempPath)
# (test <- adjInflation(datMacro=datMV, time="Date", g_start=date("2015-02-28"), g_stop=date("2022-12-31"), Inflation_Growth="M_Inflation_Growth"))
# rm(datMV, test); gc()


# -- Function to convert NaN-values or infinite values within a vector to the given value 
Treat_NaN <- function(vec, replaceVal=0) {
  vec[is.nan(vec) | is.infinite(vec)] <- replaceVal
  return (vec)
}


# -- Function for applying Winsorisation for dealing with outliers
# All values above (and or below) a certain quantile(s) are assigned the value(s) of that specific quantile(s).
# Input: [x]: A real valued vector
# Output: A vector containing values betweem the specified lower and upper quantile of the input vector.
winsorise <- function(x,lower_quant=0.025,upper_quant=0.975){
  # --- Testing purposes
  #x <- 1560*rbeta(10000, shape1=1, shape2=20); hist(x)
  #lower_quant <- 0.05; upper_quant <- 0.90
  # - Obtainnig the upper and lower quantiles of the distribution and creating indicator variables for applying the winsorisation
  wins_ind_low <- as.numeric(x<quantile(x,lower_quant))
  wins_ind_up <- as.numeric(x>quantile(x,upper_quant))
  wins_ind_bet <- abs(wins_ind_low + wins_ind_up-1)
  
  # - Applying the winsorisation
  z <- wins_ind_bet*x + wins_ind_low*quantile(x,lower_quant) + wins_ind_up*quantile(x,upper_quant)
  
  # - Returning the winsorised vector
  return(z)
  
  #hist(z); rm(x,z,lower_quant,upper_quant,wins_ind_low,wins_ind_up,wins_ind_bet)
}




# -------- Interleaving & Interpolation functions

# -- Coalescing function to facilitate data fusion between two given vectors
# Input: two scalar values (x & y) that may have selective missingness in either side (left: x; right: y)
# Output: Returns the non-missing side. If both are non-missing, then returns the (given) preference.
interleave <- function(x,y, na.value = as.integer(NA), pref='X') {
  # ensure require(dplyr)
  case_when(!is.na(x) & is.na(y) ~ x,
            is.na(x) & !is.na(y) ~ y,
            is.na(x) & is.na(y) ~ na.value,
            x == y ~ x,
            x != y & pref=='X' ~ x,
            x != y & pref=='Y' ~ y,
  )
}


# -- Missing value Treatment: Interpolate the values between two known non-missing points
# Assumes all missingness are 'encased' between two known points.
# Input: [given]: a time series possibly with some missing values for which we like to interpolate;
#     [shouldRollBackward]: If the first element is missing, should we try to fix this by 'back-interpolating'
#       from the first non-missing point found?;
#     [SilenceWarnings]: Self-explanatory;
#     [shouldRollForward]: When there is only a single non-missing element, should we simply copy that value forward?
# Output: Linearly interpolated vector
interPol <- function(given, shouldRollForward=T, shouldRollBackward=T, SilenceWarnings=T) {
  
  # -- Testing conditions
  #given <- macro_data_hist$Inflation # for testing
  #given <- as.vector(subset(macro_data, Scenario=="Historic")[order(Date_T), RealGDP_Growth_yoy])
  #unique(macro_data$Scenario)
  #given <- as.vector(subset(macro_data, Scenario=="Baseline")[order(Date_T), rbqn_rb5339q])
  # (given <- as.vector(subset(macro_data, Scenario=="SevereStress")[order(Date_T), Consumption_Level_1q]))
  
  
  # first, check if there are any non-missing element
  if (all(is.na(given))) {
    # yes, there is, so just return the same input values and throw a warning (if allowed)
    if (SilenceWarnings==F) {
      warning("All data is missing, returning NA throughout..") 
    }
    return(given)
  }
  
  # second, check if there is any missing value; if so, then exit the function
  if (all(!is.na(given))) {
    return(given)
  }
  
  # third, check if first value is missing, which can hamper our interpolation procedure
  if (is.na(given[1])) {
    # yup, so should we try to fix this by 'back-interpolating' based on the first set of 2 non-missing values in the series?
    if (shouldRollBackward == T) {
      
      start.point <- 1 # starting point for filling in interpolated vaues at the end of this procedure
      
      # find first non-missing value in the series, which will be our 'ending value' for interpolating backwards
      end.point <- which(!is.na(given))[1]-1 # position before first non-missing element
      end.val <- given[end.point+1] # first non-missing element
      
      # we need to find second non-missing value and perform an 'interim' interpolation so that we have a one-period value
      # by which to change [end.val] backwards to [start.point] at the 'same speed' (as an assumption)
      start.point2 <- which(!is.na(given))[1]+1 # position after first non-missing element
      start.val2 <- given[start.point2-1] # first non-missing element
      end.point2 <- which(!is.na(given))[2]-1 # position before second non-missing element
      end.val2 <- given[end.point2+1] # second non-missing element
      
      # interpolate across this range, including the two values as outer bounds (therefore add 2 to the interpolation length)
      # note that (end.point - start.point + 1) denotes the length of this missingness-episode
      inter.vals <- seq(from=start.val2, to=end.val2, length.out = end.point2 - start.point2 + 1 + 2)
      
      # - might as well linearly interpolate here (saving a computing cycle of the while loop later on) ..
      # delete the first and last observation (they are the outer values outside of the missingness range)
      # and assign these interpolated values to the given vector
      given[start.point2:end.point2] <- inter.vals[2:(end.point2 - start.point2 + 2)]
      
      # check if we have non-zero elements at both sides
      if (start.val2 == 0 & end.val2 == 0) {
        # yes, so by-pass this treatment and just fill with 0s
        given[start.point:end.point] <- rep(0, end.point-start.point + 1)
      } else {
        
        # get interpolation 'speed'
        speed <- diff(given[start.point2:(start.point2+1)]) / given[start.point2]
        # given[start.point2]*(1+speed) # test
        
        # 'discount' the value backwards from the first non-missing value, using the previously calculated speed as the 'discount rate'
        for (i in end.point:start.point ) {
          given[i] <- given[i+1]*(1+speed)^(-1)
        } 
      }
      
    } else {
      # no we cannot. So throw error and exit
      stop("Error: Base assumption violated - First observation is missing, cannot interpolate. Exiting ..") 
    }
  }
  
  # repeat until no more missingness in given vector
  while ( any(is.na(given)) ) {
    
    # -- testing conditions
    #given <- c(2,NA,NA,5,NA,NA,8) # works
    #given <- c(2,NA,NA,5,6, NA,NA, 9) # works
    #given <- c(2,NA,NA,5,6, NA, NA, NA)
    
    # find the indices of all missing observations
    miss.ind <- which(is.na(given))
    
    # find "episodes" of missingness in these indices, since there may be more than 1 episode in the general case,
    # for which we need to repeat this procedure.
    # 1. Do this by first isolating the cases where the lagged differences are greater than 1
    # 2. Add 1 to these found positions to move to the "initial starting points" of the next episode in succession
    # 3. Pre-fix this vector with '1' to re-include the first 'episode' that was deselected previously
    # 4. Given this vector of indices (of indices), return starting positions again
    episode.starting.times <- miss.ind[c(1, which(diff(miss.ind) > 1) + 1)]
    
    # - check if we have data points outside of the first episode from which to interpolate
    # get staring point of first episode of missingness
    start.point <- episode.starting.times[1]
    # get ending point of first episode (got to test first if we have multiple episodes and diverge logic from there)
    if (length(episode.starting.times) > 1) {
      # we have multiple episodes. Therefore, scan the series from missingness's start up to the first non-missing element, then minus 1
      # add this to the starting point, minus 1 to exclude the first missing value (otherwise we are double-counting it when adding this range)
      end.point <- start.point + (Position(function(x) {!is.na(x)}, x=given[start.point:(episode.starting.times[2]-1)] ) - 1) - 1
    } else {
      # we don't have multiple episodes. Therefore, take last known missingness index
      end.point <- miss.ind[length(miss.ind)]
    }
    
    # given the starting and ending points for the actual interpolation, test for non-missing data outside of this range from 
    # which we need to interpolate
    if (!is.na(given[start.point-1]) & !is.na(given[end.point+1])) {# returns true if we can interpolate (no missingness outside of range)
      start.val <- given[start.point-1]
      end.val <- given[end.point+1]
      # interpolate across this range, including the two values as outer bounds (therefore add 2 to the interpolation length)
      # note that (end.point - start.point + 1) denotes the length of this missingness episode
      inter.vals <- seq(from=start.val, to=end.val, length.out = (end.point - start.point + 1) + 2)
      # delete the first and last observation (they are the outer values outside of the missingness range)
      # and assign these interpolated values to the given vector
      given[start.point:end.point] <- inter.vals[2:(end.point - start.point + 2)]
      
    } else {
      # assumption violated or episode's length = 1. Check if we can simply replace NAs with last known value in either case?
      if (shouldRollForward == T){
        if (SilenceWarnings==F) {
          warning("Base assumption violated - no available data outside of missingness range from which to interpolate. Rolling values forward instead ..")
        }
        # by definition, we must have a non-missing first element (should start.point >= 2)
        start.val <- given[start.point-1]
        given[start.point:end.point] <- rep(start.val, (end.point - start.point + 1)) # just repeat for the length of the missingness episode
        
      } else {
        # no we cannot. So throw error and exit
        stop("Error: Base assumption violated - no available data outside of missingness range from which to interpolate. Exiting ..") 
      }
    }
    
  }
  
  return(given)
}





# -------- Scaling functions

# -- A few scaling functions to standardize given vectors unto a uniform scale
# Input: [given]: a real-valued vector
# Output: standardized vector
### NOTE: the shifting parameter can be useful when the given vector contains excessive zero-values.

# 1) Range-based scaler | vectors will have equal ranges (min-max); 
scaler <- function(given, shift=TRUE){
  if (shift==T){
    output <- (given - min(given,na.rm=T)) / (max(given, na.rm=T) - min(given, na.rm=T))
  } else {
    output <- (given) / (max(given, na.rm=T) - min(given, na.rm=T))
  }
  return(output)
}
# 2) Z-score/normalized scaler | vectors should roughly be N(0,1) distributed
scaler.norm <- function(given, shift=TRUE){
  # (given <- as.vector(subset(macro_data_hist1, Scenario=="Baseline")$DebtToIncome_Rate)) # for testing
  if (shift==T){
    output <- (given - mean(given,na.rm=T)) / (sqrt(var(given,na.rm=T)))
  } else {
    output <- (given) / (sqrt(var(given,na.rm=T)))
  }
  # check for NaN values (which can result if there is 0 variance)
  if (all(is.na(output))) {
    # just assign the central value, in this case, 0
    output <- rep(0, length(output))
  }
  return(output)
}





# -------- Transformation functions

# -- Yeo-Johnson transformation function
# A function for applying Yeo-Johnson transformation to a given vector. The optimal transformation is selected based on either a normal log-likelihood
# function of the transformed vector. The optimal transformation is thus chosen based on the best approximation to normality.
# Input: [x]: a real-valued vector
# Output:vector transformed with an optimal power transformation
transform_yj <- function(x, bound_lower=-2, bound_upper=2, lambda_inc=0.5, verbose=FALSE, plotopt=FALSE, plotqq=FALSE, norm_test=FALSE, loss_func="loglik"){
  # --- Unit test
  # x <- 1560*rbeta(10000, shape1=1, shape2=20); hist(x)
  # x <- rgamma(10000, shape=2, rate=5); hist(x) 
  # bound_lower<--5; bound_upper<-5; lambda_inc<-0.5; verbose<-FALSE; plotopt<-TRUE; plotqq<-TRUE; norm_test=TRUE; loss_func<-"loglik"
  
  # - Preliminaries
  require(MASS) # Ensure the requried pacakage in loaded of the Box-Cox function
  lambda_search <- seq(from=bound_lower, to=bound_upper, by=lambda_inc) # Check if lambda_search exists and if not, assign a value: This parameter specifies the search space to obatin the optimal power transformation in th boxcox function
  
  # - Ensuring the plotting area is ready for possible plots 
  par(mfcol=c(1,1))
  
  # - Selecting the optimal lambda1 parameter based on the choice of the loss function
  lambda <- boxcox((x-min(x)+0.000001)~1,lambda=seq(from=bound_lower,to=bound_upper,by=lambda_inc), plotit=FALSE)
  lambda_opt <- lambda$x[which.max(lambda$y)]
  lambda_yj <- lambda_opt
  
  # - Applying the Yeo-Johnson transformation with the optimal lambda parameter
  con1 <- as.integer(x>=0)*(lambda_yj!=0)
  con2 <- as.integer(x>=0)*(lambda_yj==0) 
  con3 <- as.integer(x<0)*(lambda_yj!=2)
  con4 <- as.integer(x<0)*(lambda_yj==2)
  
  y <- ((con1*x+1)^lambda_yj-1)/(ifelse(lambda_yj!=0,lambda_yj,1)) + log(con2*x+1) - ((-con3*x+1)^(2-lambda_yj)-1)/(2-ifelse(lambda_yj!=2,lambda_yj,1)) - log(-con4*x+1)
  
  
  # - Reporting the optimal lambda parameter, as well as the corresponding log-likelihood
  cat('\nNOTE:\tThe optimal power-transformation is lambda1 = ', lambda_yj)
  cat('\n \tThe optimal transformation has a log-likelihood =', round(max(lambda$y[which(lambda$x==lambda_yj)])), '\n')
  
  # - Plotting two qq-plots to show the normality of the given vector (x) before and after the optimal transformation 
  if(plotqq==TRUE){
    par(mfcol= c(1,2)); qqnorm(x); qqline(x, distribution = qnorm); qqnorm(y); qqline(y, distribution = qnorm);
  }
  
  # - Conducting a KS test for normality
  if(norm_test){
    ks_test_x <- ks.test(x, "pnorm")$p.value
    ks_test_y <- ks.test(y, "pnorm")$p.value
    
    cat('\nNOTE:\tThe KS-test for normality on the un-transformed data yields a p-value of ', ks_test_x)
    cat('\n \tThe KS-test for normality on the transformed data yields a p-value of ', ks_test_y)
  }
  
  # - Return the transformed vector
  cat('\n \n')
  return(y)
  #rm(x,y,bound_lower,bound_upper, verbose, plotopt, plotqq, lambda1, lambda2, lambda_search, norm_test, loss_func)
}
# Some more testing conditions
# transform_yj(x=1560*rbeta(10000, shape1=1, shape2=20), bound_lower=4, plotopt=TRUE)
# transform_yj(x=1560*rbeta(10000, shape1=1, shape2=20),bound_lower=-2,bound_upper=2,lambda_inc=0.1, plotopt=TRUE, plotqq=TRUE, norm_test=TRUE)


# -- Box-Cox transformation function
# A function for applying Box-Cox transformation to a given vector. The optimal transformation is selected based on either a normal log-likelihood
# function or the skewness of the transformed vector. The optimal transformation is thus chosen based on the best approximation to normality.
# Input: [x]: a real-valued vector
# Output:vector transformed with an optimal power transformation
transform_bc <- function(x, bound_lower=-2, bound_upper=2, lambda_inc=0.5 , anchor1=FALSE, verbose=FALSE, plotopt=FALSE, plotqq=FALSE, norm_test=FALSE, loss_func="loglik"){
  # --- Testing purposes
  #x <- 1560*rbeta(10000, shape1=1, shape2=20); hist(x)
  #x <- rgamma(10000, shape=2, rate=5); hist(x) 
  #anchor1 <- FALSE; bound_lower<--5; bound_upper<-5; lambda_inc<-0.5; anchor1<-FALSE; verbose<-FALSE; plotopt<-TRUE; plotqq<-TRUE; norm_test=TRUE; loss_func<-"skewness";
  
  # - Preliminaries
  require(MASS) # Ensure the requried pacakage in loaded of the Box-Cox function
  lambda_search <- seq(from=bound_lower, to=bound_upper, by=lambda_inc) # Check if lambda_search exists and if not, assign a value: This parameter specifies the search space to obatin the optimal power transformation in th boxcox function
  
  # - Selecting lambda2: parameter for ensuring that all values of x > 0 since Box-Cox can't handle values x>=0
  if (anchor1==TRUE){ # Anchor the minimum of the distribution at 1 (i.e., shift entire distribution such that the minimum is one)
    lambda2 <- 1 - min(x)
  } else if (anchor1==FALSE) {# Do not anchor the minimum of the distribution at 1, but ensure that values are still > 0.0000001 to enable the Box-Cox transformations to be successfully executed with in boxcox()
    if (min(x)==0){
      lambda2 <- 0.0000001
      if(verbose) cat("\nNOTE: Zero values detected in x, minimum of distribution shifted to 0.0000001")
    } else {
      lambda2 <- ifelse(min(x)<0,-min(x) + 0.0000001,0)
    }
  } 
  
  # Ensuring the plotting area is ready for possible plots 
  par(mfcol=c(1,1))
  
  # - Selecting the optimal lambda1 parameter based on the choice of the loss function
  if(loss_func=="loglik"){ # Optimal transformation chosen based on log-likelihood function (maximum)
    opt_lambda <- boxcox(x+lambda2 ~ 1, lambda = lambda_search, plotit = plotopt)
    lambda1 <- as.numeric(sprintf("%.1f",opt_lambda$x[which.max(opt_lambda$y)])) # Rounding to the nearest 0.1 for easier interpretation.
  } else if (loss_func=="skewness"){ # Optimal transformation chosen based on skewness of distribution (absolute minimum)
    skew <- rep(0,length(lambda_search))
    for (i in 1:length(lambda_search)){
      transformed_x <- if(lambda_search[i]!=0) ((x+lambda2)^lambda_search[i]-1)/lambda_search[i] else log(x+lambda2)
      skew[i] <- (skewness(transformed_x))
    }
    if (plotopt) {plot(x=lambda_search,y=skew, type="l", xlab=bquote(lambda), ylab="Skewness"); abline(v=lambda_search[which.min(abs(skew))], lty=3)}
    lambda1 <- as.numeric(sprintf("%.1f",lambda_search[which.min(skew)])) 
  }
  
  # - Applying the optimal power transformation using the Box-Cox transformation as developed by Box & Cox (1964); URL = https://www.semanticscholar.org/paper/An-Analysis-of-Transformations-Box-Cox/6e820cf11712b9041bb625634612a535476f0960
  if (lambda1 != 0){
    y <- ((x+lambda2)^lambda1-1)/lambda1
  } else if (lambda1 == 0){
    y <- log(x+lambda2)
  }

  # - Reporting the optimal lambda1 and lambda2 parameters, as well as the corresponding log-likelihood
  cat('\nNOTE:\tThe optimal power-transformation is lambda1 = ', lambda1)
  cat('\n \tThe distribution was shifted with lambda2 =', sprintf("%.0000001f",lambda2))
  if (loss_func=="log_lik") cat('\n \tThe optimal transformation has a log-likelihood =', round(max(opt_lambda$y)), '\n')
  if (loss_func=="skewness") cat('\n \tThe optimal transformation has a skewness =', min(abs(skew)), '\n')
  
  # - Plotting two qq-plots to show the normality of the given vector (x) before and after the optimal transformation 
  if(plotqq==TRUE){
    par(mfcol= c(1,2)); qqnorm(x); qqline(x, distribution = qnorm); qqnorm(y); qqline(y, distribution = qnorm);
  }
  
  # - Conducting a KS test for normality
  if(norm_test){
    ks_test_x <- ks.test(x, "pnorm")$p.value
    ks_test_y <- ks.test(y, "pnorm")$p.value
    
    cat('\nNOTE:\tThe KS-test for normality on the un-transformed data yields a p-value of ', ks_test_x)
    cat('\n \tThe KS-test for normality on the transformed data yields a p-value of ', ks_test_y)
  }
  
  # - Return the transformed vector
  cat('\n \n')
  return(y)
  #rm(x,y,bound_lower,bound_upper, verbose, plotopt, plotqq, lambda1, lambda2, lambda_search, norm_test, loss_func)
}

# Some more testing conditions
#bc_transform(x=1560*rbeta(10000, shape1=1, shape2=20), anchor1 = FALSE, bound_lower=4, plotopt=TRUE)
#bc_transform(x=1560*rbeta(10000, shape1=1, shape2=20),bound_lower=-2,bound_upper=2,lambda_inc=0.1,  anchor1=FALSE, plotopt=TRUE, plotqq=TRUE, norm_test=TRUE)




# -------- Performance measurement functions for fitted models

# - function to return the VIF of the variables
# Input: model formula
# Output: VIF-value
multicollinearity <- function(model){
  vif_model <- VIF(model)
  return(vif_model)
}

# - function to compute the accuracy, precision, and recall
# Input:
# 3 parameters relating to the model that must be analysed, the data set that
# contains the target variable, and the cut-off imposed that flags accounts as 
# SICR, i.e., the cut-off imposed on the posterior probability output
# Output:
# A data table containing the accuracy, precision, and recall of the model under consideration
accuracy_prec_rec <- function(model, dat_set, cut_off){
  
  predict_model <- predict(model, dat_set, type = 'response')
  
  # Confusion matrix
  confusion_mat_model <- table(dat_set$DefaultStatus1_LS, predict_model > cut_off)
  
  # Accuracy
  accuracy_model <- sum(diag(confusion_mat_model)) / sum(confusion_mat_model)
  
  # Precision
  tp <- confusion_mat_model[2, 2]
  fp <- confusion_mat_model[1, 2]
  
  Prec <- tp / (tp + fp)
  
  # Recall
  tp <- confusion_mat_model[2, 2]
  fn <- confusion_mat_model[2, 1]
  
  Rec <- tp / (tp + fn)
  
  # Create a dataset to return the measures
  Performance_model <- data.table(Measures = c("Accuracy", "Precision", "Recall"),
                                  Values = c(accuracy_model, Prec, Rec))
  
  return(Performance_model)
  
}






