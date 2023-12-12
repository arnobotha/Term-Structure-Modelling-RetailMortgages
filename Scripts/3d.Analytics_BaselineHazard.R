# ================================== Analytics : Baseline Hazards ==================================
# Analytical script for exploring different baseline hazard estimators.
# --------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller
# --------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 2d.Data_Fusion.R
#
# -- INPUT:
#   1) creditdata_final3  | Prepared credit dataset from scirpt 2d.
#
# -- OUTPUT:
#   1) Some graphs used for comparison purposes and statistics.
# ==================================================================================================

# ------ 1. Preliminaries
# - Confirm prepared credit data is loaded into memory
if (!exists('datCredit_real')) {  unpack.ffdf(paste0(genPath,"creditdata_analytics"), tempPath)}

# --- Remove some unnecessary variables to improve performance
colnames_datCredit_real <- colnames(datCredit_real)
datCredit_real <- datCredit_real %>% subset(select = -c(which(colnames_datCredit_real == "M_Inflation_Growth_1"):which(colnames_datCredit_real == "M_Repo_Rate_Vol_12")))

# --- Fitting a Cox model
# - Grouping the table according to PerfSpell_Key and TimeInPerfSpell
datCredit_real <- datCredit_real %>% group_by(PerfSpell_Key, TimeInPerfSpell)

# - Fitting the model
cph_1 <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1) ~
                 Principal + Balance + PerfSpell_Num
               , id= PerfSpell_Key, data = datCredit_real)
# Principal     :  coef = -0.00000260788; exp(coef) = 0.99999739212; se(coef) = 0.00000006218; z = -41.94; p < 0.0000000000000002;
# Balance       :  coef = 0.00000238474; exp(coef) =  1.00000238475; se(coef) = 0.00000006461; z = 36.91; p < 0.0000000000000002;
# PerfSpell_Num :  coef = 0.78524638625; exp(coef) =  2.19294718621; se(coef) = 0.00683701232; z = 114.85; p < 0.0000000000000002;

# --- Convert dataset back to datatable
setDT(datCredit_real)

# --- Creating an datatable, not in counting format
datSurv <- subset(datCredit_real, PerfSpell_Counter==1,
                  select=c("LoanID", "PerfSpell_Key", "PerfSpell_Num", "PerfSpell_Counter", "TimeInPerfSpell", "PerfSpell_Age",
                           "PerfSpellResol_Type_Hist")); gc()

# --- Creating an survival object
surv_obj <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, event=PerfSpellResol_Type_Hist=="Defaulted",
                         type="counting") ~ 1, id=PerfSpell_Key, data=datSurv)

# --- Graphing parameters
chosenFont <- "Cambria"
dpi <- 180




# ------ 2. Cumulative hazard using the Nelson-Aalen estimator
# --- Creating a counting process style dataset
dat_count1 <- datCredit_real[,.N, by=TimeInPerfSpell]; names(dat_count1) <- c("Time", "AtRisk") # This includes defaulted observations
dat_count2 <- datCredit_real[DefaultStatus1==1, .N, by=TimeInPerfSpell]; names(dat_count2) <- c("Time", "Defaulted")
NelAal <- merge(dat_count1, dat_count2, by=c("Time")); rm(dat_count1, dat_count2)


# --- Computing the hazard and cumulative hazard
NelAal[, Haz := Defaulted/ AtRisk]
NelAal[, CumHaz := cumsum(Haz)]




# ------ 3. Cumulative hazard using the Breslow estimator
# --- Computing the risk score of each observation
datCredit_real[, RiskScore := predict(cph_1, type = "risk")]

# --- Computing the sum of the risk scores at each period
dat_count3 <- datCredit_real[, sum(RiskScore, rm.na=T), by=list(TimeInPerfSpell)]
colnames(dat_count3) <- c("Time", "Sum_RiskScore")

# --- Computing the number of defaults at each period
dat_count4 <- datCredit_real[DefaultStatus1==1, .N, by=list(TimeInPerfSpell)]
colnames(dat_count4) <- c("Time", "Defaults")

# --- Computing the Breslow estimator for the cumulative hazard
Bres <- merge(dat_count3, dat_count4, by=c("Time"), all.x=T); gc() # Merging the two datasets
Bres <- Bres[!is.na(Defaults),] # Removing all periods where there were no defaults
Bres[, Cum_Defaults := cumsum(Defaults)] # Getting the cumulative sums of default at each default period
Bres[, Haz := Defaults/Sum_RiskScore] # Computing the hazard
Bres[, CumHaz := cumsum(Haz)] # Computing the cumulative hazard




# ------ 4. Comparing the Nelson-Aalen and Breslow Estimator
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
# --- Compute Kaplan-Meier survival estimates (product-limit) for default-event | Spell-level with right-censoring
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as performing spell key, with no stratification
# --- Creating a dataset containing the at-risk popoulation and defaults
dat_count5 <- datCredit_real[,.N, by=TimeInPerfSpell]; names(dat_count5) <- c("Time", "AtRisk") # This includes defaulted accounts
dat_count6 <- datCredit_real[DefaultStatus1==1, .N, by=TimeInPerfSpell]; names(dat_count6) <- c("Time", "Defaulted")
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
# --- Computing the KP estimator
dat_count7 <- datCredit_real[, sum(RiskScore, rm.na=T), by=list(TimeInPerfSpell)]; colnames(dat_count7) <- c("Time", "Sum_RiskScore_All")
dat_count8 <- datCredit_real[DefaultStatus1==1, sum(RiskScore, rm.na=T), by=list(TimeInPerfSpell)]; colnames(dat_count8) <- c("Time", "Sum_RiskScore_Default")
KalPre <- merge(KapMei[, c("Time", "AtRisk", "Defaulted")], dat_count7, all.x=T, by="Time"); KalPre <- merge(KalPre, dat_count8, all.x = T, by="Time")

# --- Computing the survival estimates
# - Default times where there is more than one default
Times <- KalPre[Defaulted>0 & Time!=1, Time]

# - Looping over each default time (with more than one event)
for (i in Times){
  FindRoot = function(alpha){
    a <- KalPre[Time==Times[i], Sum_RiskScore_All] # The first quantity to estimate: Total sum of the risk scores for the associated default time
    b <- sum(datCredit_real[TimeInPerfSpell==Times[i] & DefaultStatus1==1, RiskScore]/(1 - alpha^datCredit_real[TimeInPerfSpell==Times[i] & DefaultStatus1==1, RiskScore])) # The second quantity to estimate: 
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


