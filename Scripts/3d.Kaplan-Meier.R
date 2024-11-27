# =========================== KAPLAN-MEIER ESTIMATION ===================================
# Description of script 
#
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Exploring different approaches to modelling the term-structure 
#                of default risk.
# SCRIPT AUTHOR(S): Dr Arno Botha, Roelinde Bester
# ---------------------------------------------------------------------------------------
# INPUT:
#
# OUTPUT:
#
# =======================================================================================

# ------ 0. Preliminaries

# - Confirm prepared credit data is loaded into memory
if (!exists('datCredit_real')) {
  unpack.ffdf(paste0(survPath,"creditdata_final3"), tempPath)
}

# - Create a path to which the graphs must be saved
graphPath_real <- "C:/Users/R5422965/OneDrive - FRG/Default Survival Modelling/Graphs"

# --- 1. Graphing setup
trainGraphs <- "training_slc" # folder location
portName <- "FNB SLC"

# - Subset only performance spells and execute executions | Memory optimisation to aid subsequent modelling
datCredit_real_perf <- subset(datCredit_real, ExclusionID==0 & !is.na(PerfSpell_Num))
rm(datCredit_real); gc()

# - graphing parameters
col.v <- brewer.pal(10, "Paired")[c(10)]

# --- 2. Default survival modelling with left-truncated spells | Spell-level | Multiple spells
# Each subject has a single spell-level observation that summarize that spell thus far; each subject may have multiple spells
# Competing risks are latent, i.e., treated as right-censored
# Spell correctly specified with left-truncated start time and 0 for unaffected spells

# --- Sample accordingly | First record per spell
datSurv <- subset(datCredit_real_perf, ExclusionID==0 & !is.na(PerfSpell_Num) & PerfSpell_Counter==1 & Partition == "Observed", 
                  select=c("LoanID", "PerfSpell_Key", "PerfSpell_Num", "TimeInPerfSpell","PerfSpell_Age", 
                           "PerfSpellResol_Type_Hist")); gc()

# --- Compute Kaplan-Meier survival estimates (product-limit) for default-event | Spell-level with right-censoring & left-truncation
# Left-truncated spells will have TimeInPerfSpell > 1 while the rest will have TimeInPerfSpell = 1, both of
#   which are adjusted during modelling by subtracting 1 (since stop must be > start for KM-estimation)
#   Note: At TimeInPerfSpell = 1, already 1 month has lapsed, therefore, the true "start" of a normal account is actually 
#     at TimeInPerfSpell = 0. Therefore, subtracting one is defensible.
# All competing events preclude the main event from happening and are therefore considered as censored
# ID is set as performing spell key, with no stratification
kmPerf_default_real_spell1c <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, event=PerfSpellResol_Type_Hist=="Defaulted",
                                            type="counting") ~ 1, id=PerfSpell_Key, data=datSurv)
summary(kmPerf_default_real_spell1c)$table
surv_summary(kmPerf_default_real_spell1c)



# --- Graphing survival and related quantities

# -- Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(10)] # for S(t)
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
sSpan <- 0.1; vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)
mainEventName <- "Default"
chosenFont <- "Cambria"

### AB [2024-11-26]: I made some aesthetic changes given Bernard's MSc-study, though the following was not tested on our data.

# -- Survival probability, S(t)=y
(gsurv1c_a <- ggsurvplot(kmPerf_default_real_spell1c, fun="pct", conf.int=T, legend="none", 
                         break.time.by=round(max(kmPerf_default_real_spell1c$time)/8), palette=vCol,
                         xlab = bquote(Discrete~time~italic(t)*" (months) in performing spell: Multi-spell"),
                         ylab = bquote(Survival~probability~"["*.(mainEventName)*"]"*~italic(S(t))*": spell-level (Kaplan-Meier)"),
                         xlim=c(0, max(kmPerf_default_real_spell1c$time)+1), surv.median.line = "hv", censor=F, 
                         ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                         tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                         cumevents=T, cumevents.title="Cumulative number of events", 
                         cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                         risk.table.title = "Number in (% of) sample at risk of main event", font.family=chosenFont, fontsize=2.5))

# -- Cumulative event/lifetime probability: ;, , so F(t)=1-S(t)
(gsurv1c_b <- ggsurvplot(kmPerf_default_real_spell1c, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                         break.time.by=round(max(kmPerf_default_real_spell1c$time)/8), palette=vCol,
                         xlab = bquote(Discrete~time~italic(t)*" (months) in performing spell: Multi-spell"),
                         ylab = bquote(Cumulative~lifetime~distribution~"["*.(mainEventName)*"]"*~italic(F(t))*": spell-level (Kaplan-Meier)"),
                         xlim=c(0, max(kmPerf_default_real_spell1c$time)+1), censor=F, 
                         ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                         tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                         risk.table.title = "Number in (% of) sample at risk of main event", font.family=chosenFont, fontsize=2.5))

# -- Cumulative hazard: f(t)=-log(y), so H(t) = -log(S(t)). Discrete-variant: f(t)=-sum{ln(1-h(t))}
(gsurv1c_c <- ggsurvplot(kmPerf_default_real_spell1c, fun="cumhaz",conf.int=T, legend="none", 
                         break.time.by=round(max(kmPerf_default_real_spell1c$time)/8), palette=vCol,
                         xlab = bquote(Discrete~time~italic(t)*" (months) in performing spell: Multi-spell"),
                         ylab = bquote(Cumulative~hazard~distribution~"["*.(mainEventName)*"]"*~italic(H(t))*": spell-level (Nelson-Aalen)"),
                         xlim=c(0,  max(kmPerf_default_real_spell1c$time)+1), censor=F, 
                         ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                         tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                         risk.table.title = "Number in (% of) sample at risk of main event", font.family=chosenFont, fontsize=2.5))


# -- Discrete baseline hazard function: h(t)
# - create plotting data object
haz_dat <- data.table(Time=kmPerf_default_real_spell1c$time, AtRisk_n=kmPerf_default_real_spell1c$n.risk, 
                      Event_n = kmPerf_default_real_spell1c$n.event, Censored_n=kmPerf_default_real_spell1c$n.censor,
                      hazard=kmPerf_default_real_spell1c$n.event/kmPerf_default_real_spell1c$n.risk, 
                      CumulHazard = kmPerf_default_real_spell1c$cumhaz, # Nelson-Aalen estimator
                      Group="1",Surv_KM = kmPerf_default_real_spell1c$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1- Surv_KM_Disc]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]
# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
plot(kmExample$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences

# - Graph object for shorter time, informed by previous graphs
(gsurv1c_d <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=hazard)) + theme_minimal() +
    geom_line(linetype="solid", colour=vCol2[1]) + geom_point(colour=vCol2[1]) + 
    geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    labs(y=bquote(plain(Estimated~hazard*" function ["*.(mainEventName)*"]"*~italic(h(t))*": spell-level (Kaplan-Meier)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2[2], labels=label.v) + 
    scale_fill_manual(name="", values=vCol2[2], labels=label.v) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))

### RESULTS: Most continuous-time variants are exactly equal to the discrete-time variants, except cumulative hazard
# Hazard as h(t) = f(t) / n(t) is equal to hazard deduced recursively from S(t) as h(t) = 1 - S(t)/S(t-1), 
# while h(1) = 1 - S(1)

# - save plots
dpi <- 120 # need to decrease size for risk tables' text
ggsave(print(gsurv1c_a,newpage=F), file=paste0(graphPath_real,"/SurvFig1c_a-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(print(gsurv1c_b,newpage=F), file=paste0(graphPath_real,"/SurvFig1c_b-", mainEventName,"_CumulEvent-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
ggsave(print(gsurv1c_c,newpage=F), file=paste0(graphPath_real,"/SurvFig1c_c-", mainEventName,"_CumulHaz-NelsonAalen-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)
dpi <- 140 # reset
ggsave(gsurv1c_d, file=paste0(graphPath_real,"/SurvFig1c_d-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)

# - cleanup
rm(datSurv); gc()