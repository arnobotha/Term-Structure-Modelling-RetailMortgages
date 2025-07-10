# ======================================= GENERAL ANALYTICS ============================================
# Analyses at risk population over unique spell ages
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




# ----------------- 1. Load & prepare data for analysis

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]








# ------ 2. Analytics: At-risk proportion over unique failure times

# --- Preliminaries
# - Create pointer to the appropriate data object 
datCredit <- rbind(datCredit_train, datCredit_valid)
sMaxSpellAge_graph <- 300 # max for [PerfSpell_Age] for graphing purposes


# --- Estimate survival rate of the main event S(t) = P(T >=t) for time-to-event variable T
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_Default <- survfit(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=PerfSpell_Event==1, type="counting") ~ 1, 
                      id=PerfSpell_Key, data=datCredit)

# - Create survival table using surv_summary(), from the subsampled set
(datSurv_sub <- surv_summary(km_Default)) # Survival table
datSurv_sub <- datSurv_sub %>% rename(Time=time, AtRisk_n=`n.risk`, Event_n=`n.event`, Censored_n=`n.censor`, SurvivalProb_KM=`surv`) %>%
  mutate(Hazard_Actual = Event_n/AtRisk_n) %>% 
  mutate(CHaz = cumsum(Hazard_Actual)) %>% # Created as a sanity check
  mutate(EventRate = Hazard_Actual*shift(SurvivalProb_KM, n=1, fill=1)) %>%  # probability mass function f(t)=h(t).S(t-1)
  filter(Event_n > 0 | Censored_n >0) %>% as.data.table()
setorder(datSurv_sub, Time)
datSurv_sub[,AtRisk_perc := AtRisk_n / max(AtRisk_n, na.rm=T)]

# - Create two time series
datGraph <- rbind(datSurv_sub[,list(Time, Value=AtRisk_perc, Type="a_AtRisk")],
                  datSurv_sub[,list(Time, Value=SurvivalProb_KM, Type="b_Survival")])

datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) spell-time"]

# - Aesthetic engineering
chosenFont <- "Cambria"
vCol <- brewer.pal(8, "Set1")[c(3,1)]
vLabel <- c("a_AtRisk"=bquote("At-risk % of max("*italic(n)[t]*")"),
            "b_Survival"=bquote("Survival probability "*italic(S(t))*"(Kaplan-Meier)"))

(gAtRisk <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=Value, group=Type)) + theme_minimal() +
    labs(y="Value (%)", 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    #geom_point(aes(colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(colour=Type, linetype=Type), linewidth=0.5) +
    # Scales and options
    facet_grid(FacetLabel ~ .) + 
    scale_colour_manual(name="", values=vCol, labels=vLabel) + 
    scale_linetype_discrete(name="", labels=vLabel) + 
    scale_shape_discrete(name="", labels=vLabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))

# - Save plot
dpi <- 180 # reset
ggsave(gAtRisk, file=paste0(genFigPath, "Kaplan-Meier/AtRiskPop.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Cleanup
rm(gAtRisk, datGraph, km_Default, datSurv_sub)
