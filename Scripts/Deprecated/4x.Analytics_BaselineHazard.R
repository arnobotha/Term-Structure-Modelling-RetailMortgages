# ================================== Analytics : Baseline Hazards ==================================
# Ancillary script for exploring different baseline hazard estimators.
# --------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller, Dr Arno Botha (AB)
# --------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R
#
# -- INPUT:
#   1) creditdata_final3  | Prepared credit dataset from scirpt 2d.
#
# -- OUTPUT:
#   1) Some graphs used for comparison purposes and statistics.
# ==================================================================================================


### AB [2025-05-05]: Marked as deprecated due to fundamental errors. Useful parts mined and refactored elsewhere


# ------ 1. Preliminaries
# - Confirm prepared credit data is loaded into memory
if (!exists('datCredit_train_PWPST')) {  unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath)}

# - Initialize variables
vars2 <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Build model based on variables
cph_1 <- coxph(as.formula(paste0("Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1) ~ ", 
                                 paste(vars2,collapse=" + "), " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, datCredit_train_PWPST, ties="efron")
summary(cph_1); concordance(cph_1)



# ------ 2. Cumulative hazard using the Nelson-Aalen estimator
# --- Creating a counting process style dataset
dat_count1 <- datCredit_train_PWPST[,.N, by=TimeInPerfSpell]; names(dat_count1) <- c("Time", "AtRisk") # This includes defaulted observations
dat_count2 <- datCredit_train_PWPST[DefaultStatus1==1, .N, by=TimeInPerfSpell]; names(dat_count2) <- c("Time", "Defaulted")
NelAal <- merge(dat_count1, dat_count2, by=c("Time")); rm(dat_count1, dat_count2)


# --- Computing the hazard and cumulative hazard
NelAal[, Haz := Defaulted/ AtRisk]
NelAal[, CumHaz := cumsum(Haz)]
plot(NelAal$CumHaz)





### AB [2025-05-05]: Times ought to be unique failure times. Here is another more correct algorithm
### AB [2025-05-05]: Start of algorithm
dat <- datCredit_train_PWPST[!is.na(PerfSpell_Num), list(PerfSpell_Key, SpellNum=PerfSpell_Num, 
                                                         Start=TimeInPerfSpell-1, Stop=TimeInPerfSpell, Event=DefaultStatus1,
                                                         SpellAge=PerfSpell_Age,ResolType=PerfSpellResol_Type_Hist)]
setkey(dat, PerfSpell_Key, Stop)
vFailTimes <- sort(unique(dat$Stop[dat$Event==1]))
vCHaz <- numeric(length(vFailTimes)); CHaz <- 0
ind <- 1
for (t in vFailTimes) {
  # number of defaults
  d_t <- sum(dat$Stop==t & dat$Event==1)
  # number at risk
  n_t <- sum(dat$Start <= t & dat$Stop > t)
  # update cumulative hazard
  CHaz <- CHaz + d_t / n_t
  vCHaz[ind] <- CHaz
  ind <- ind + 1
}
plot(vFailTimes, vCHaz, xlim=c(0,600))
lines(NelAal$CumHaz, col="red")

# - Compare with cumulative hazard from KM-analysis
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_All <- survfit(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=PerfSpell_Event==1, type="counting") ~ 1, 
                  id=PerfSpell_Key, data=datCredit_train_PWPST)
CHaz_KM <- -log(km_All$surv)

# Plot differences
plot(km_All$time, CHaz_KM, type = "s", col = "blue", xlim = c(0, 600),
     ylab = "Cumulative Hazard", xlab = "Time", main = "Nelson-Aalen vs KM-derived Hazard")
lines(vFailTimes, vCHaz, type = "s", col = "red", lty = 2)
legend("topleft", legend = c("KM (from survfit)", "Nelson-Aalen"), col = c("blue", "red"), lty = 1:2)
### RESULTS: Basically equal

### CONCLUSION: KM-derived cumulative hazard is a smoothed version of the Nelson-Aalen estimator
### AB [2025-05-05]: End of algorithm






# ------ 3. Cumulative hazard using the Breslow estimator
# --- Computing the risk score of each observation
datCredit_train_PWPST[, RiskScore := predict(cph_1, type = "risk", newdata=datCredit_train_PWPST)]

# --- Computing the sum of the risk scores at each period
dat_count3 <- datCredit_train_PWPST[, sum(RiskScore, rm.na=T), by=list(TimeInPerfSpell)]
colnames(dat_count3) <- c("Time", "Sum_RiskScore")

# --- Computing the number of defaults at each period
dat_count4 <- datCredit_train_PWPST[DefaultStatus1==1, .N, by=list(TimeInPerfSpell)]
colnames(dat_count4) <- c("Time", "Defaults")

# --- Computing the Breslow estimator for the cumulative hazard
Bres <- merge(dat_count3, dat_count4, by=c("Time"), all.x=T); gc() # Merging the two datasets
Bres <- Bres[!is.na(Defaults) & !is.na(Time),] # Removing all periods where there were no defaults
Bres[, Cum_Defaults := cumsum(Defaults)] # Getting the cumulative sums of default at each default period
Bres[, Haz := Defaults/Sum_RiskScore] # Computing the hazard
Bres[, CumHaz := cumsum(Haz)] # Computing the cumulative hazard
plot(Bres$Time, Bres$CumHaz)




### AB [2025-05-05]: Times ought to be unique failure times for Breslow-estimator
### AB [2025-05-05]: Start of algorithm
datCredit_train_PWPST[, RiskScore2 := predict(cph_1, type = "exp", newdata=datCredit_train_PWPST)]
dat <- datCredit_train_PWPST[!is.na(PerfSpell_Num), list(PerfSpell_Key, SpellNum=PerfSpell_Num, 
                                                         Start=TimeInPerfSpell-1, Stop=TimeInPerfSpell, Event=DefaultStatus1,
                                                         SpellAge=PerfSpell_Age,ResolType=PerfSpellResol_Type_Hist,
                                                         RiskScore=RiskScore2,
                                                         Stratum=PerfSpell_Num_binned)]
setkey(dat, PerfSpell_Key, Stop)

vFailTimes <- sort(unique(dat$Stop[dat$Event==1]))
vCHaz0 <- numeric(length(vFailTimes))

for (i in seq_along(vFailTimes)) {
  t <- vFailTimes[i]
  # number of defaults
  d_t <- sum(dat$Stop==t & dat$Event==1)
  # Risk set
  r_t <- sum(dat$RiskScore[dat$Start <= t & dat$Stop > t])
  # Update cumulative baseline hazard \delta(t) = d_t / r_t | Breslow's estimator
  vCHaz0[i] <- (ifelse(i > 1, vCHaz0[i-1], 0)) + (d_t / r_t)
}
plot(vFailTimes, vCHaz0)
### RESULTS: Flawed. ChatGPT is of no use in debugging code. Effort abandoned
# Needs to estimate H0 per stratum, but abandoned.

# - Compare with survfit()-output
objSurvfit <- survfit(cph_1)
plot(objSurvfit$time, objSurvfit$cumhaz)

### AB [2025-05-05]: End of algorithm





# ------ 4. Comparing the Nelson-Aalen and Breslow Estimator

### AB [2025-05-05]: Flawed estimates of either NA and Breslow estimator. Deprecated.

# --- Restricting the comparison to only include defaults occuring befor TimeInPerfSpell (Time) 300
NelAal <- NelAal[Time<=300, list(Time, CumHaz, Haz)]; NelAal <- NelAal[Time!=1,]
Bres <- Bres[Time<=300, list(Time, CumHaz, Haz)]

# --- Testing how similar the two estimators are
all.equal(NelAal$CumHaz, Bres$CumHaz)
### RESULTS:~ Not all equal, with mean relative difference of 0.189125

# --- Getting the data in a long format for using ggplot
dat_plot <- rbind(NelAal[, Est := "a_NelAal"],
                  Bres[, Est:= "b_Bres"])

# --- Plotting the differences
# - MAE between NA and Breslow estimator
mae_CumHaz_NA_Bres <- sum(abs(dat_plot[Est=="a_NelAal", CumHaz] - dat_plot[Est=="b_Bres", CumHaz]))/dat_plot[Est=="a_NelAal",.N]
mae_Haz_NA_Bres <- sum(abs(dat_plot[Est=="a_NelAal", Haz] - dat_plot[Est=="b_Bres", Haz]))/dat_plot[Est=="a_NelAal",.N]

# - Graphing parameters
col.v <- brewer.pal(9, "Spectral")[c(3,9)]
label.v <- c("a_NelAal"="NA-Estimator",
             "b_Bres"="BR-Estimator")
shape.v <- c(16,4); linetype.v <- c("solid", "dotted")
chosenFont <- "Cambria"; dpi <- 180

# - Cumulative Hazard Comparison
(g_comp1 <- ggplot(dat_plot[Est %in% c("a_NelAal", "b_Bres")], aes(x=Time, y=CumHaz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(H(t))), family=chosenFont) +  
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    #annotations
    annotate(geom="text", x=120, y=0.9,
             label=paste("MAE between NA- and BR-estimator = ", sprintf("%.4f",mae_CumHaz_NA_Bres)),
             family=chosenFont, size=4) +
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v, labels=label.v) + 
    scale_shape_manual(name="Estimator", values=shape.v, labels=label.v) + 
    scale_linetype_manual(name="Estimator", values=linetype.v, labels=label.v) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v))

# Saving the plot
ggsave(file = "NA_BR_CumHaz.png", plot = g_comp1, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)

# - Hazard Comparison
(g_comp2 <- ggplot(dat_plot, aes(x=Time, y=Haz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(h(t))), family=chosenFont) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    #annotations
    annotate(geom="text", x=125, y=0.0135,
             label=paste("MAE between NA- and BR-estimator = ", sprintf("%.4f",mae_Haz_NA_Bres)),
             family=chosenFont, size=4) +
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v, labels=label.v) + 
    scale_shape_manual(name="Estimator", values=shape.v, labels=label.v) + 
    scale_linetype_manual(name="Estimator", values=linetype.v, labels=label.v) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v))

# Saving the plot
ggsave(file = "NA_BR_Haz.png", plot = g_comp2, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)




# ------ 5. Cumulative Hazard using the Kaplan-Meier Estimator (actually known as the Peterson estimator from Peterson (1977) => H(t) = -log(S(t)))

### AB [2025-05-05]: Deprecated since comparison already performed earlier.

# --- Compute Kaplan-Meier survival estimates (product-limit) for default-event | Spell-level with right-censoring
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as performing spell key, with no stratification
# --- Creating a dataset containing the at-risk popoulation and defaults
dat_count5 <- datCredit_train_PWPST[,.N, by=TimeInPerfSpell]; names(dat_count5) <- c("Time", "AtRisk") # This includes defaulted accounts
dat_count6 <- datCredit_train_PWPST[DefaultStatus1==1, .N, by=TimeInPerfSpell]; names(dat_count6) <- c("Time", "Defaulted")
KapMei <- merge(dat_count5, dat_count6, all.x=T, by=c("Time")); rm(dat_count5, dat_count6)
KapMei[is.na(Defaulted), Defaulted := 0] # Correcting for missing values

# --- Computing the survival estimates
KapMei[, Alpha := (AtRisk - Defaulted)/AtRisk] # Proportion of accounts surviving in each period
KapMei[, KM_Surv := cumprod(Alpha)] # Proportion of accounts surviving in each consecutive period

# --- Computing the cumulative hazard and hazard
KapMei[, CumHaz := -log(KM_Surv)]
KapMei[, Haz := CumHaz - lag(CumHaz)]

# --- Comparing the NA- and KM-estimator
all.equal(dat_plot[Est=="a_NelAal",CumHaz], KapMei[Time <= 300 & Time!=1, ]$CumHaz)
### RESULTS:~ Mean relative difference: 0.002010412
plot(dat_plot[Est=="a_NelAal",CumHaz] - KapMei[Time <= 300 & Time!=1, ]$CumHaz)
### RESULTS:~ Clearly a difference between the estimators

# --- Plotting the differences
dat_plot <- rbind(dat_plot , KapMei[, Est := "c_KapMei"][Time <= 300 & Time != 1, list(Time, CumHaz, Haz, Est)])

# - MAE between NA- and KM estimator
mae_CumHaz_NA_KM <- sum(abs(dat_plot[Est=="a_NelAal", CumHaz] - dat_plot[Est=="c_KapMei", CumHaz]))/dat_plot[Est=="a_NelAal",.N]
mae_Haz_NA_KM <- sum(abs(dat_plot[Est=="a_NelAal", Haz] - dat_plot[Est=="c_KapMei", Haz]))/dat_plot[Est=="a_NelAal",.N]

# - Graphing parameters
label.v2 <- c("a_NelAal"="NA",
              "c_KapMei"="KM")

# - Cumulative Hazard Comparison
(g_comp3 <- ggplot(dat_plot[Est %in% c("a_NelAal", "c_KapMei")], aes(x=Time, y=CumHaz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(H(t))), family=chosenFont) +  
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    #annotations
    annotate(geom="text", x=120, y=0.9,
             label=paste("MAE between NA- and KM-estimator = ", sprintf("%.4f",mae_CumHaz_NA_KM)),
             family=chosenFont, size=4) +
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v, labels=label.v2) + 
    scale_shape_manual(name="Estimator", values=shape.v, labels=label.v2) + 
    scale_linetype_manual(name="Estimator", values=linetype.v, labels=label.v2) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v2))

# Saving the plot
ggsave(file = "NA_KM_CumHaz.png", plot = g_comp3, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)

# - Hazard Comparison
(g_comp4 <- ggplot(dat_plot[Est %in% c("a_NelAal", "c_KapMei")], aes(x=Time, y=Haz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(h(t))), family=chosenFont) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    #annotations
    annotate(geom="text", x=125, y=0.0135,
             label=paste("MAE between NA- and KM-estimator = ", sprintf("%.4f",mae_Haz_NA_KM)),
             family=chosenFont, size=4) +
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v, labels=label.v2) + 
    scale_shape_manual(name="Estimator", values=shape.v, labels=label.v2) + 
    scale_linetype_manual(name="Estimator", values=linetype.v, labels=label.v2) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v2))

# Saving the plot
ggsave(file = "NA_KM_Haz.png", plot = g_comp3, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)




# ------ 6. Kalbfleisch-Prentice (KP) Estimator

### AB [2025-05-05]: Flawed since times ought to be unique failure times. Deprecated.


# --- Computing the KP estimator
dat_count7 <- datCredit_train_PWPST[, sum(RiskScore, rm.na=T), by=list(TimeInPerfSpell)]; colnames(dat_count7) <- c("Time", "Sum_RiskScore_All")
dat_count8 <- datCredit_train_PWPST[DefaultStatus1==1, sum(RiskScore, rm.na=T), by=list(TimeInPerfSpell)]; colnames(dat_count8) <- c("Time", "Sum_RiskScore_Default")
KalPre <- merge(KapMei[, c("Time", "AtRisk", "Defaulted")], dat_count7, all.x=T, by="Time"); KalPre <- merge(KalPre, dat_count8, all.x = T, by="Time")

# --- Computing the survival estimates
# - Default times where there is more than one default
Times <- KalPre[Defaulted>0 & Time!=1, Time]

# - Looping over each default time (with more than one event)
for (i in Times){
  FindRoot = function(alpha){
    a <- KalPre[Time==Times[i], Sum_RiskScore_All] # The first quantity to estimate: Total sum of the risk scores for the associated default time
    b <- sum(datCredit_train_PWPST[TimeInPerfSpell==Times[i] & DefaultStatus1==1, RiskScore]/(1 - alpha^datCredit_train_PWPST[TimeInPerfSpell==Times[i] & DefaultStatus1==1, RiskScore])) # The second quantity to estimate: 
    return(a - b)
  }  
  
  KalPre[Time==Times[i]-1, Alpha := uniroot(f = FindRoot, interval = c(0, 0.9999), tol = 0.00001)$root]
  
}

# - Correcting for the first survival quantity (should be one)
KalPre[Time==1, Alpha := 1]

# - Correcting for a failed iteration
KalPre[36,]$Alpha <- KalPre[35,Alpha] + (KalPre[37,Alpha] - KalPre[35,Alpha])/2

# - Computing the KP-estimator
KalPre[, KP_Surv := cumprod(Alpha)]

# --- Computing the cumulative hazard and hazard
KalPre[, CumHaz := -log(KP_Surv)]
KalPre[, Haz := CumHaz - lag(CumHaz)]

# --- Comparing the KM- and KP-estimator
all.equal(dat_plot[Est=="c_KapMei",CumHaz], KalPre[Time <= 300 & Time!=1,]$CumHaz)
### RESULTS:~ Mean relative difference: 0.1902455
plot(dat_plot[Est=="c_KapMei",CumHaz] - KalPre[Time <= 300 & Time !=1,]$CumHaz)
### RESULTS:~ Clearly a difference between the estimators

# --- Plotting the differences
dat_plot <- rbind(dat_plot, KalPre[, Est:="d_KalPre"][Time <= 300 & Time !=1, c("Time", "CumHaz", "Haz", "Est")])

# - MAE between NA and Breslow estimator
mae_CumHaz_KM_KP <- sum(abs(dat_plot[Est=="c_KapMei", CumHaz] - dat_plot[Est=="d_KalPre", CumHaz]))/dat_plot[Est=="c_KapMei",.N]
mae_Haz_KM_KP <- sum(abs(dat_plot[Est=="c_KapMei", Haz] - dat_plot[Est=="d_KalPre", Haz]))/dat_plot[Est=="c_KapMei",.N]

# - Graphing parameters
label.v3 <- c("c_KapMei"="KM",
              "d_KalPre"="KP")

# - Cumulative Hazard Comparison
(g_comp5 <- ggplot(dat_plot[Est %in% c("c_KapMei", "d_KalPre")], aes(x=Time, y=CumHaz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(H(t))), family=chosenFont) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    #annotations
    annotate(geom="text", x=120, y=0.9,
             label=paste("MAE between KM- and KP-estimator = ", sprintf("%.4f",mae_CumHaz_KM_KP)),
             family=chosenFont, size=4) +
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v, labels=label.v3) + 
    scale_shape_manual(name="Estimator", values=shape.v, labels=label.v3) + 
    scale_linetype_manual(name="Estimator", values=linetype.v, labels=label.v3) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v3))

# Saving the plot
ggsave(file = "KM_KP_CumHaz.png", plot = g_comp5, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)

# - Hazard Comparison
(g_comp6 <- ggplot(dat_plot[Est %in% c("c_KapMei", "d_KalPre")], aes(x=Time, y=Haz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(h(t))), family=chosenFont) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    #annotations
    annotate(geom="text", x=125, y=0.0135,
             label=paste("MAE between KM- and KP-estimator = ", sprintf("%.4f",mae_Haz_KM_KP)),
             family=chosenFont, size=4) +
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v, labels=label.v3) + 
    scale_shape_manual(name="Estimator", values=shape.v, labels=label.v3) + 
    scale_linetype_manual(name="Estimator", values=linetype.v, labels=label.v3) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v3))

# Saving the plot
ggsave(file = "KM_KP_Haz.png", plot = g_comp6, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)




# ------ 7. NA- vs BR- vs KM- vs KP
# - Graphing parameters
col.v4 <- brewer.pal(9, "Spectral")[c(1,3,9,11)]
shape.v4 <- c(16,4,1,18); linetype.v4 <- c("solid", "dotted", "dashed", "dotdash")
label.v4 <- c("a_NelAal"="NA",
              "b_Bres"="BR",
              "c_KapMei"="KM",
              "d_KalPre"="KP")

# - Cumulative Hazard Comparison
(g_comp7 <- ggplot(dat_plot, aes(x=Time, y=CumHaz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(H(t))), family=chosenFont) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v4, labels=label.v4) + 
    scale_shape_manual(name="Estimator", values=shape.v4, labels=label.v4) + 
    scale_linetype_manual(name="Estimator", values=linetype.v4, labels=label.v4) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v4))

# Saving the plot
ggsave(file = "NA_BR_KM_KP_CumHaz.png", plot = g_comp7, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)

# - Hazard Comparison
(g_comp8 <- ggplot(dat_plot, aes(x=Time, y=Haz)) + 
    theme_minimal() +
    labs(x=bquote(italic(t)), y=bquote(italic(h(t))), family=chosenFont) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Est, linetype=Est)) + 
    geom_point(aes(colour=Est, shape=Est)) + 
    # facets & scale options
    scale_colour_manual(name="Estimator", values=col.v4, labels=label.v4) + 
    scale_shape_manual(name="Estimator", values=shape.v4, labels=label.v4) + 
    scale_linetype_manual(name="Estimator", values=linetype.v4, labels=label.v4) + 
    scale_size_manual(name="Estimator", values=c(0.5,0.5), labels=label.v4))

# Saving the plot
ggsave(file = "NA_BR_KM_KP_Haz.png", plot = g_comp8, device = 'png', path = genFigPath_CumHaz, height = 1600/dpi, width = 2000/dpi, dpi=dpi)


