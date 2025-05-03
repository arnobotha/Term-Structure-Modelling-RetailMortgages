# Cumulative cases/ dynamic controls vs Incident Cases vs dynamic controls
#-------------------------

# ------- 1, Preliminaries
# --- Confirm prepared training and validation dataset is loaded into memory
if (!exists('dat_train')) unpack.ffdf(paste0(genPath,"dat_train"), tempPath)
if (!exists('dat_valid')) unpack.ffdf(paste0(genPath,"dat_valid"), tempPath)

# --- Remove some unessary variables to improve performance
dat_train <- dat_train[1:1000000,-(63:200)]
dat_valid <- dat_valid[1:1000000,-(63:200)]


# --- Transforming the Default status into a numeric variable to ensure compatibility with time-dependent AUC/ ROC functions
dat_train[,DefaultStatus1 := as.numeric(DefaultStatus1)]; dat_train[,DefaultStatus1 := ifelse(DefaultStatus1==1,0,1)]

dat_valid[,DefaultStatus1 := as.numeric(DefaultStatus1)]; dat_valid[,DefaultStatus1 := ifelse(DefaultStatus1==1,0,1)]
# unique(dat_valid[PerfSpell_Age>=36,PerfSpell_Key])
# dat_sample <- dat_valid[PerfSpell_Key %in% c("3000000157761_1", "3000000199847_1", "3000000231791_1", "3000000380002_1", "3000000550374_1"),]
# --- Grouping the data according to PerfSpell_Key and TimeInPerfSpell
dat_train <- dat_train %>% group_by(PerfSpell_Key, TimeInPerfSpell)
dat_valid <- dat_valid %>% group_by(PerfSpell_Key, TimeInPerfSpell)
# dat_sample <- dat_sample %>% group_by(PerfSpell_Key, TimeInPerfSpell)
# --- Load model into memory
if (!exists('cph_def_com_scaled3')) unpack.ffdf(paste0(genPath,"cph_def_fin"), tempPath)

# --- Obtaining predictions using the final default Cox model on the validation dataset.
dat_valid$lp <- predict(cph_def_com_scaled3, dat_valid, type = "lp")
# dat_sample$lp <- predict(cph_def_com_scaled3, dat_sample, type = "lp")

# ------- 2. Incident cases/ dynamic controls
# Incident sensitivity is the probability that an individual has a marker (predicted) value 
#   greater than c among those individuals who experienced the event AT TIME t. 
# Dynamic specificity is the probability that an individual has a marker value less than or 
#   equal to among the individuals that remain event-free AT TIME t. 
# I/D is appropriate to apply when the exact event time is known and we want to discriminate 
#   between individuals experiencing the event and those event-free AT A GIVEN EVENT TIME t.

# - Assessing the predictive power
# Define a helper function to evaluate at various t
risksetROC_helper <- function(dat_explain, Input_Names = c("TimeInPerfSpell", "DefaultStatus1"), marker = cph_lp, predict.time=12) {
  risksetROC::risksetROC(Stime        = dat_explain[[Input_Names[1]]],  
                         status       = dat_explain[[Input_Names[2]]],                 
                         marker       = marker,                             
                         entry        = dat_explain[[Input_Names[1]]]-1,                              
                         predict.time = predict.time,
                         method       = "Cox",
                         plot = FALSE)
}

predict.time <- 3

risksetROC_data <- data.frame(Predict_Time = NULL, AUC = NULL, FP = NULL, TP = NULL, Marker = NULL)
start_time_ID <- Sys.time()
for (i in 1:length(predict.time)){
  temp <- risksetROC_helper(dat_explain = dat_valid, Input_Names = c("TimeInPerfSpell", "DefaultStatus1"), marker = dat_valid$lp, predict.time = predict.time[i])
  risksetROC_data <- rbind(risksetROC_data, data.frame(Predict_Time = predict.time[i], AUC = temp$AUC, FP = temp$FP, TP = temp$TP))
}
end_time_ID <- Sys.time()
(RunTime_ID <- end_time_ID - start_time_ID)
rm(temp)
risksetROC_data <- risksetROC_data %>% arrange(Predict_Time, FP, TP) # Arranging the dataset according to the given vector of prediction times (and then by FPs and TPs)

output <- data.frame(predict_times = predict.time,
                     auc = unique(risksetROC_data$AUC))
### RESULTS:~ AUC-value at 3-months = 0.9850614
###           Runtime for one iteration = 1.255321 minutes

# --- Graphs
# - Side-by-side graphs
g_auc <- ggplot(data = risksetROC_data, mapping = aes(x = FP, y = TP)) +
         geom_line(col = 'blue') +
         geom_abline(col='red') +
         geom_label(data = risksetROC_data %>% dplyr::select(Predict_Time,AUC) %>% unique,
                    mapping = aes(label = sprintf("%.3f", AUC)), x = 0.5, y = 0.5) +
         facet_wrap(vars(Predict_Time)) +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
               legend.key = element_blank(),
               plot.title = element_text(hjust = 0.5),
               strip.background = element_blank()) +
         xlab("FP") + ggtitle('ROC Curve(s) for the given CPH Model')




# ------ 3. Cumulative sensitivity and dynamic specificity (C/D)
# Cumulative sensitivity is the probability that an individual has a marker (predicted) value 
#   greater than c among those individuals who experienced the event DURING OR BEFORE TIME t. 
# Dynamic specificity is the probability that an individual has a marker value less than or 
#   equal to among the individuals that remain event-free AT TIME t. 
# C/D is appropriate to apply when the exact event time is known and we want to discriminate 
#   between individuals experiencing the event and those event-free AT A GIVEN EVENT TIME t....???

### NOTE:~  Test by reducing the sample size (may be necessary to adjust and consider fixed sampling sizes)
# test <- dat_sample
# test$lp <- predict(cph_def_com_scaled3, test, type = "lp")
# --- KM Estimator
start_time_CD <- Sys.time()
survivalROC(Stime = test$TimeInPerfSpell,
            status = test$DefaultStatus1,
            marker = test$lp,
            entry = test$TimeInPerfSpell-1,
            predict.time = 1,
            method = 'KM')
end_time_CD <- Sys.time()
(RunTime_CD <- end_time_CD - start_time_CD)
### NaN produced for all TPs. Investigate further...


# --- NN Estimator
test <- dat_valid[100000:101000,]
nobs <- nrow(test)

start_time_CD <- Sys.time()
survivalROC(Stime = test$TimeInPerfSpell,
            status = test$DefaultStatus1,
            marker = test$lp, ### NOTE:~ Investigate marker (lp vs risk score); Look at which one the model was built for.
            entry = test$TimeInPerfSpell-1,
            predict.time = 1,
            method = 'NNE',
            lambda = 0.25*nobs^(-0.20))
end_time_CD <- Sys.time()
(RunTime_CD <- end_time_CD - start_time_CD)
# Error in if (n > 0) s0 <- s0 * (1 - d/n) : missing value where TRUE/FALSE needed
### NOTE:~ It might be that DefaultStatus1 needs to be a factor and not numerical => Looks like it can take both factors and numerical values.
# Error in order(cut.values) : argument 1 is not a vector.
### Note:~ NaN still produced using a small sample of 5 PerformanceSpell_Keys.



survivalROC_helper <- function(t) {
  survivalROC(Stime = dat_sample$TimeInPerfSpell,
              status = dat_sample$DefaultStatus1,
              marker = dat_sample$lp,
              entry = dat_sample$TimeInPerfSpell-1,
              predict.time = t,
              method       = "NNE",
              span = 0.25 * nrow(ovarian)^(-0.20))
}
## Evaluate every 180 days
survivalROC_data <- tibble(t = 12 * c(1,2)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           tibble(obj[c("cut.values","TP","FP")])
         }))





library(tidyverse)
## https://github.com/tidyverse/tibble/issues/395
options(crayon.enabled = FALSE)
## Used for the dataset.
library(survival)
## Used for visualizaiton.
library(survminer)
## Load the Ovarian Cancer Survival Data
data(ovarian)
## Turn into a data_frame
ovarian <- tibble(ovarian)
## Plot
ggsurvplot(survfit(Surv(futime, fustat) ~ 1,
                   data = ovarian),
           risk.table = TRUE,
           break.time.by = 180)

## Fit a Cox model
coxph1 <- coxph(formula = Surv(futime, fustat) ~ pspline(age, df = 4) + factor(resid.ds) +
                  factor(rx) + factor(ecog.ps),
                data    = ovarian)
## Obtain the linear predictor
ovarian$lp <- predict(coxph1, type = "lp")
ovarian

library(survivalROC)
## Define a helper functio nto evaluate at various t
survivalROC_helper <- function(t) {
  survivalROC(Stime        = ovarian$futime,
              status       = ovarian$fustat,
              marker       = ovarian$lp,
              entry = 0,
              predict.time = t,
              method       = "NNE",
              span = 0.25 * nrow(ovarian)^(-0.20))
}
## Evaluate every 180 days
survivalROC_data <- tibble(t = 180 * c(1,2,3,4,5,6)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           tibble(obj[c("cut.values","TP","FP")])
         }))

