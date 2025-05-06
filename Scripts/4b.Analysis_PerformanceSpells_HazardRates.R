# ================================= HAZARD RATE ANALYSIS =========================================
# Analysis performed on the hazard rate of different performance spells. The 
# analysis aids in the decision making about a suitable cut-off point beyond 
# which all performance spells are grouped together.
# # Refer to LifetimePD-TermStructure-RecurrentEvents for varieties of KM-analysis across 
# including/excluding left-truncation and first-spell vs multi-spell.
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB), Bernard Scheepers (BS)
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R

# -- Inputs:
#   - datCredit_real | Prepared from script 2f.
#
# -- Outputs:
#   - Cumulative baseline hazard rates by performance spells (2 different groupings)
#   - Event probabilities
#   - Hazard rates
#   - Kaplan-Meyer analysis by performance spells
#   - datSurv objects | Respective to each setting, containiing survival, cumulative hezard, & event probabilities
# ------------------------------------------------------------------------------------------------------




# -------- 0 Preliminarnies

# - Generic parameters
mainEventName <- "Default"

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)





# -------- 1. Kaplan-Meier analysis across all performance spells | Left-truncation included

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_All <- survfit(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=PerfSpell_Event==1, type="counting") ~ 1, 
                  id=PerfSpell_Key, data=datCredit_real)
summary(km_All)$table # overall summary statistics
plot(km_All, conf.int = T) # survival curve

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_All$time[km_All$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_All$surv*diff(c(0,km_All$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_All$surv[max(which(km_All$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_All$time), "): \n\t",
           percent(km_All$surv[max(which(km_All$time <= max(km_All$time)))], accuracy=0.01)))



# --- Cumulative Lifetime Distribution F(t)=1-S(t)

# - Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_All)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
maxTime <- max(km_All$time)
chosenFont <- "Cambria"

# - Graphing logic: Entire domain of time to event
gsurv_Ft <- ggsurvplot(km_All, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                               break.time.by=round(maxTime)/8, palette=vCol,
                               xlab = bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell"),
                               ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                               xlim=c(0, maxTime+1), censor=F, 
                               ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                               tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                               cumevents=T, cumevents.title="Cumulative number of events", 
                               cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                               risk.table.title = "Number in (% of) sample at risk of main event",
                               font.family=chosenFont, fontsize=2.5,
                               surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  annotate("text", x = medLife, y = 0.5,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft$plot <- gsurv_Ft$plot  +    
    geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
    geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
    annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
             vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
### RESULTS: Based on survival analysis, about 100% of the dataset was consumed at t=300; We shall henceforth remove outliers purely
# for graphing purposes
maxPeriod <- 300

# - Save graph
dpi <- 190
ggsave(gsurv_Ft$plot, file=paste0(genFigPath, "Kaplan-Meier/CumulLifetimeProb_", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
dpi <- 185 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft,newpage=F), file=paste0(genFigPath,"Kaplan-Meier/CumulLifetimeProb_", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc_RiskTable.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")


# - Graphing logic: Limited domain of time to event [maxPeriod]
gsurv_Ft2 <- ggsurvplot(km_All, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                       break.time.by=round(maxPeriod)/8, palette=vCol,
                       xlab = bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell"),
                       ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                       xlim=c(0, maxPeriod), censor=F, 
                       ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                       tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                       cumevents=T, cumevents.title="Cumulative number of events", 
                       cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                       risk.table.title = "Number in (% of) sample at risk of main event",
                       font.family=chosenFont, fontsize=2.5,
                       surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft2$plot <- gsurv_Ft2$plot  +    
  annotate("text", x = medLife, y = 0.51,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=1, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft2$plot <- gsurv_Ft2$plot  +    
  geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  annotate("text", x = percTimepoint, y = 1 - survPercentile+0.01, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
           vjust=-0.5, hjust=1, color = "black", size = 2.5, family = chosenFont)

# - Save graph
ggsave(print(gsurv_Ft2,newpage=F), file=paste0(genFigPath,"Kaplan-Meier/CumulLifetimeProb_", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc_RiskTable_LimPeriod.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")



# --- Create data and graphing objects
# - Create survival table using surv_summary()
(datSurv <- surv_summary(km_All)) # Survival table
datSurv <- datSurv %>% rename(Time=time, AtRisk_n=`n.risk`, Event_n=`n.event`, Censored_n=`n.censor`, SurvivalProb_KM=`surv`) %>%
  mutate(Hazard_Actual = Event_n/AtRisk_n, Hazard_Actual2 = 1 - SurvivalProb_KM/shift(SurvivalProb_KM,n=1,fill=1)) %>% 
  mutate(Hazard_Actual2 = ifelse(is.na(Hazard_Actual2), 0, Hazard_Actual2)) %>% # Handle NaN-values
  mutate(CHaz = cumsum(Hazard_Actual),
         CHaz2 = -log(SurvivalProb_KM),# Created as a sanity check
         SurvivalProb_KM_Disc = cumprod(1-Hazard_Actual)) %>% # Created as a sanity check
  mutate(EventRate = Hazard_Actual*shift(SurvivalProb_KM, n=1, fill=1)) %>%  # probability mass function f(t)=h(t).S(t-1)
  filter(Event_n > 0 | Censored_n >0) %>% as.data.table()

# - conduct sanity checks
all.equal(datSurv$Hazard_Actual, datSurv$Hazard_Actual2) # Should be TRUE
all.equal(datSurv$SurvivalProb_KM, datSurv$SurvivalProb_KM_Disc) # Should be TRUE
all.equal(km_All$cumhaz, datSurv$CHaz) # Should be TRUE
all.equal(datSurv$CHaz, datSurv$CHaz2)
plot(datSurv$CHaz - datSurv$CHaz2, type="b")
### RESULTS: CHaz2 is very similar to CHaz, derived fundamentally from KM-estimate of S(t), though kept for comparative purposes

# - Remove sanity checks
datSurv[, c("SurvivalProb_KM_Disc","CHaz2", "Hazard_Actual2") := NULL]

# - Distributional analyses
# Hazard rate
describe(datSurv$Hazard_Actual); hist(datSurv$Hazard_Actual, breaks="FD")
plot(datSurv$Hazard_Actual, type="b")
### RESULTS: No apparent shape though increasingly affected by outliers as time increases,
# so we should put less stock in those right-most results on the x-axis.

# Event rate
describe(datSurv$EventRate); hist(datSurv$EventRate, breaks="FD")
plot(datSurv$EventRate, type="b")
### RESULTS: Shows a U-shaped distribution, as expected



# --- Graphing the hazard rate h(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.2
smthHazard_Act <- loess(datSurv$Hazard_Actual ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthHazard_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthHazard_Act, newdata=datSurv, se=T)

# - Add smoothed hazard to graphing object
datSurv[, Hazard_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, Hazard_spline_upper := Hazard_spline + critVal*vPredSmth$se.fit]
datSurv[, Hazard_spline_lower := Hazard_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] hazard'*~italic(h(t))*' with 95% CI'))
datSurv[, Group := "1"] # only necessary for aesthetics when using geom_smooth() during plotting

# - Create main graph 
(gsurv_ht <- ggplot(datSurv[Time<= maxPeriod, ], aes(x=Time, y=Hazard_Actual)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=Hazard_Actual), colour=vCol2[2]) + 
    geom_line(aes(y=Hazard_Actual), colour=vCol2[2], linetype="solid") + 
    #geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    # Smoothed quantity
    geom_line(aes(y=Hazard_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=Hazard_spline_lower, ymax=Hazard_spline_upper), fill=vCol2[1], alpha=0.25)+
    # Scales and options
    labs(y=bquote(plain(Discrete~hazard~~italic(h(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be U-shaped, as expected

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ht, file=paste0(genFigPath, "Kaplan-Meier/DiscreteHazard_", mainEventName,"_SpellLevel_MultiSpell-LatentComp_InclLeftTrunc-AG.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# --- Graphing the event density / probability mass function f(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.3
smthEventRate_Act <- loess(datSurv$EventRate ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthEventRate_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthEventRate_Act, newdata=datSurv, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, EventRate_spline_upper := EventRate_spline + critVal*vPredSmth$se.fit]
datSurv[, EventRate_spline_lower := EventRate_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol3 <- brewer.pal(10, "Paired")[c(4,3)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] estimate'*~italic(f(t))*' with 95% CI'))

# - Create main graph 
(gsurv_ft <- ggplot(datSurv[Time <= maxPeriod,], aes(x=Time, y=EventRate)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=EventRate), colour=vCol3[2]) + 
    geom_line(aes(y=EventRate), colour=vCol3[2], linetype="solid") + 
    # Smoothed quantity
    geom_line(aes(y=EventRate_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=EventRate_spline_lower, ymax=EventRate_spline_upper), fill=vCol3[1], alpha=0.15)+
    # Scales and options
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol3, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol3, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(colour = "none"))

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "Kaplan-Meier/EventProb-", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), datSurv)

# - Housekeeping
rm(gsurv_Ft, gsurv_Ft2, gsurv_ht, gsurv_ft, km_All, median_survival, datSurv, vlabel, 
   smthEventRate_Act, smthHazard_Act, vPredSmth)





# -------- 2. Kaplan-Meier analysis on first performance spell | Left-truncation included

# - Initialize a dataset with only the first performance spell
dat <- datCredit_real[PerfSpell_Num==1,]

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_first <- survfit(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=PerfSpell_Event==1,type="counting") ~ 1, 
                  id=PerfSpell_Key, data=dat)
summary(km_first)$table # overall summary statistics
plot(km_first, conf.int = T) # survival curve

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_first$time[km_first$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_first$surv*diff(c(0,km_first$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_first$surv[max(which(km_first$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_first$time), "): \n\t",
           percent(km_first$surv[max(which(km_first$time <= max(km_first$time)))], accuracy=0.01)))



# --- Cumulative Lifetime Distribution F(t)=1-S(t)

# - Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_first)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
maxTime <- max(km_first$time)
chosenFont <- "Cambria"
dpi <- 190

# - Graphing logic: Entire domain of time to event
gsurv_Ft <- ggsurvplot(km_first, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                       break.time.by=round(maxTime)/8, palette=vCol,
                       xlab = bquote(Discrete~time~italic(t)*" (months) in spell: First-spell"),
                       ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                       xlim=c(0, maxTime+1), censor=F, 
                       ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                       tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                       cumevents=T, cumevents.title="Cumulative number of events", 
                       cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                       risk.table.title = "Number in (% of) sample at risk of main event",
                       font.family=chosenFont, fontsize=2.5,
                       surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  annotate("text", x = medLife, y = 0.5,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
### RESULTS: Based on survival analysis, about 100% of the dataset was consumed at t=300; We shall henceforth remove outliers purely
# for graphing purposes
maxPeriod <- 300

# - Save graph
ggsave(gsurv_Ft$plot, file=paste0(genFigPath, "Kaplan-Meier/CumulLifetimeProb_", mainEventName,"_SpellLevel_FirstSpell_LatentComp_InclLeftTrunc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
dpi <- 185 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft,newpage=F), file=paste0(genFigPath,"Kaplan-Meier/CumulLifetimeProb_", mainEventName,"_SpellLevel_FirstSpell_LatentComp_InclLeftTrunc_RiskTable.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")


# --- Create data and graphing objects
# - Create survival table using surv_summary()
(datSurv <- surv_summary(km_first)) # Survival table
datSurv <- datSurv %>% rename(Time=time, AtRisk_n=`n.risk`, Event_n=`n.event`, Censored_n=`n.censor`, SurvivalProb_KM=`surv`) %>%
  mutate(Hazard_Actual = Event_n/AtRisk_n, Hazard_Actual2 = 1 - SurvivalProb_KM/shift(SurvivalProb_KM,n=1,fill=1)) %>% 
  mutate(Hazard_Actual2 = ifelse(is.na(Hazard_Actual2), 0, Hazard_Actual2)) %>% # Handle NaN-values
  mutate(CHaz = cumsum(Hazard_Actual),
         CHaz2 = -log(SurvivalProb_KM),# Created as a sanity check
         SurvivalProb_KM_Disc = cumprod(1-Hazard_Actual)) %>% # Created as a sanity check
  mutate(EventRate = Hazard_Actual*shift(SurvivalProb_KM, n=1, fill=1)) %>%  # probability mass function f(t)=h(t).S(t-1)
  filter(Event_n > 0 | Censored_n >0) %>% as.data.table()

# - conduct sanity checks
all.equal(datSurv$Hazard_Actual, datSurv$Hazard_Actual2) # Should be TRUE
all.equal(datSurv$SurvivalProb_KM, datSurv$SurvivalProb_KM_Disc) # Should be TRUE
all.equal(km_first$cumhaz, datSurv$CHaz) # Should be TRUE
all.equal(datSurv$CHaz, datSurv$CHaz2)
plot(datSurv$CHaz - datSurv$CHaz2, type="b")
### RESULTS: CHaz2 is very similar to CHaz, derived fundamentally from KM-estimate of S(t), though kept for comparative purposes

# - Remove sanity checks
datSurv[, c("SurvivalProb_KM_Disc","CHaz2", "Hazard_Actual2") := NULL]

# - Distributional analyses
# Hazard rate
describe(datSurv$Hazard_Actual); hist(datSurv$Hazard_Actual, breaks="FD")
plot(datSurv$Hazard_Actual, type="b")
### RESULTS: No apparent shape though increasingly affected by outliers as time increases,
# so we should put less stock in those right-most results on the x-axis.

# Event rate
describe(datSurv$EventRate); hist(datSurv$EventRate, breaks="FD")
plot(datSurv$EventRate, type="b")
### RESULTS: Shows a U-shaped distribution, as expected



# --- Graphing the hazard rate h(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.2
smthHazard_Act <- loess(datSurv$Hazard_Actual ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthHazard_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthHazard_Act, newdata=datSurv, se=T)

# - Add smoothed hazard to graphing object
datSurv[, Hazard_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, Hazard_spline_upper := Hazard_spline + critVal*vPredSmth$se.fit]
datSurv[, Hazard_spline_lower := Hazard_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] hazard'*~italic(h(t))*' with 95% CI'))
datSurv[, Group := "1"] # only necessary for aesthetics when using geom_smooth() during plotting

# - Create main graph 
(gsurv_ht <- ggplot(datSurv[Time<= maxPeriod, ], aes(x=Time, y=Hazard_Actual)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=Hazard_Actual), colour=vCol2[2]) + 
    geom_line(aes(y=Hazard_Actual), colour=vCol2[2], linetype="solid") + 
    #geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    # Smoothed quantity
    geom_line(aes(y=Hazard_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=Hazard_spline_lower, ymax=Hazard_spline_upper), fill=vCol2[1], alpha=0.25)+
    # Scales and options
    labs(y=bquote(plain(Discrete~hazard~~italic(h(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be U-shaped, as expected

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ht, file=paste0(genFigPath, "Kaplan-Meier/DiscreteHazard_", mainEventName,"_SpellLevel_FirstSpell-LatentComp_InclLeftTrunc-AG.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# --- Graphing the event density / probability mass function f(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.3
smthEventRate_Act <- loess(datSurv$EventRate ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthEventRate_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthEventRate_Act, newdata=datSurv, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, EventRate_spline_upper := EventRate_spline + critVal*vPredSmth$se.fit]
datSurv[, EventRate_spline_lower := EventRate_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol3 <- brewer.pal(10, "Paired")[c(4,3)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] estimate'*~italic(f(t))*' with 95% CI'))

# - Create main graph 
(gsurv_ft <- ggplot(datSurv[Time <= maxPeriod,], aes(x=Time, y=EventRate)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=EventRate), colour=vCol3[2]) + 
    geom_line(aes(y=EventRate), colour=vCol3[2], linetype="solid") + 
    # Smoothed quantity
    geom_line(aes(y=EventRate_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=EventRate_spline_lower, ymax=EventRate_spline_upper), fill=vCol3[1], alpha=0.15)+
    # Scales and options
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol3, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol3, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(colour = "none"))

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "Kaplan-Meier/EventProb-", mainEventName,"_SpellLevel_FirstSpell_LatentComp_InclLeftTrunc.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), datSurv)

# - Housekeeping
rm(gsurv_Ft, gsurv_ht, gsurv_ft, km_first, median_survival, datSurv, vlabel, 
   smthEventRate_Act, smthHazard_Act, vPredSmth, dat)






# -------- 3. Comparison of cumulative hazard estimators

# --- Nelson-Aalen estimator vs Kaplan-Meier derived estimator, given that H(t) = -log(S(t))
dat <- datCredit_real[!is.na(PerfSpell_Num), list(PerfSpell_Key, SpellNum=PerfSpell_Num, 
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
  CHaz <- CHaz + d_t / n_t # Nelson-Aalen estimator
  vCHaz[ind] <- CHaz
  ind <- ind + 1
}
plot(vFailTimes, vCHaz, xlim=c(0,600))

# - Compare with cumulative hazard from KM-analysis
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_All <- survfit(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=PerfSpell_Event==1, type="counting") ~ 1, 
                  id=PerfSpell_Key, data=datCredit_real[!is.na(PerfSpell_Num),])
CHaz_KM <- -log(km_All$surv)

# Plot differences
plot(km_All$time, CHaz_KM, type = "s", col = "blue", xlim = c(0, 600),
     ylab = "Cumulative Hazard", xlab = "Time", main = "Nelson-Aalen vs KM-derived Hazard")
lines(vFailTimes, vCHaz, type = "s", col = "red", lty = 2)
legend("topleft", legend = c("KM (from survfit)", "Nelson-Aalen"), col = c("blue", "red"), lty = 1:2)
### RESULTS: Basically equal
### CONCLUSION: KM-derived cumulative hazard is a smoothed version of the Nelson-Aalen estimator

# - Cleanup
rm(km_All, dat, CHaz_KM, vCHaz, vFailTimes, ind)
