# ================================== Analytics :Competing Risks ==================================
# Analytical script for exploring the exclusion of competing risks vs treating competing risks
# as being right-censored (latent approach). Densities from various datasets are compared, 
# distributional analysis are performed and non-parametric estimators are constructed and compared.
# ------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller
# ------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 2c.Data_Enrich.R
#
# -- INPUT:
#   1) creditdata_final3  | Enriched credit dataset from scirpt 2c.
#
# -- OUTPUT:
#   1) Some graphs used for comparison purposes and statistics.
# ================================================================================================



# ------- 1. Preliminaries

# - Confirm prepared credit data is loaded into memory
if (!exists('datCredit_real')) {
  unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)
}

# - Only selecting one record per performance spell
datSurv <- subset(datCredit_real, PerfSpell_Counter==1 & PerfSpell_Age < 500,
                  select=c("LoanID","PerfSpell_Key","PerfSpell_Age", "PerfSpellResol_Type_Hist", "PerfSpell_Num", "TimeInPerfSpell")); gc()

# - Cleaning
rm(datCredit_real)



# ------ 2. Densities Overlaid
# --- Defaults vs Right Censored vs Competing Risks
# - Group competing risks together in Resol_Type
datSurv[,Resol_Type := case_when(PerfSpellResol_Type_Hist == "Defaulted" ~ "a_Default",
                                 PerfSpellResol_Type_Hist == "Censored" ~ "b_Right-censored",
                                 PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "c_Competing_risk")]
# - Small check
datSurv$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
(Resol_Type.props <- datSurv$Resol_Type %>% table() %>% prop.table())

# - Dataset containing only the competing risks, where the type of competing risk is shown by Resol_Type2
datSurv2 <- datSurv %>% subset(PerfSpellResol_Type_Hist%in%c("Settled", "Paid-up", "Written-off"))
datSurv2[, Resol_Type2 := case_when(PerfSpellResol_Type_Hist=="Settled" ~ "a_Settlement",
                                    PerfSpellResol_Type_Hist=="Paid-up" ~ "b_Paid-up",
                                    TRUE ~ "c_Written-off-Other")]
# - Small check
datSurv$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
(Resol_Type2.props <- datSurv2$Resol_Type2 %>% table() %>% prop.table())


# - Graphing Parameters
chosenFont <- "Cambria"; dpi <- 170
col.v <- brewer.pal(10, "Paired")[c(10,6,4, 2,1,8)]
col.v2 <- brewer.pal(10, "Paired")[c(10,5,3, 1,2,8)]
label.v <- c(paste0("a_Default"="Default (", round(Resol_Type.props[1]*100, digits=1), "%)"), # Need to round to the first decimal place to ensure that the prior add up to one
             paste0("b_Right-censored"="Right-censored (", round(Resol_Type.props[2]*100, digits=1), "%)"),
             paste0("c_Competing_risk"="Competing Risks (", round(Resol_Type.props[3]*100, digits=1), "%)"))
label.v2 <- c(paste0("a_Settlement"="Settlement (", round(Resol_Type2.props[1]*100, digits=1), "%)"),
              paste0("b_Paid-up"="Paid-up (", round(Resol_Type2.props[2]*100, digits=1), "%)"),
              paste0("c_Written-off-Other"="Other/Write-off (", round(Resol_Type2.props[3]*100, digits=1), "%)"))


# - Densities of resolution types overlayed
(g1_Densities_Resol_Type <- ggplot(datSurv[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type)) + theme_minimal() + 
    geom_histogram(aes(y=after_stat(density), colour=Resol_Type, fill=Resol_Type), position="identity",
                   alpha=0.75, size=0.2) + 
    geom_density(aes(colour=Resol_Type, linetype=Resol_Type), size=0.8) + 
    # facets & scale options
    labs(y=bquote(plain(Empirical~failure*' time histogram & density ')~italic(f(t))), 
         x=bquote("Performing spell ages (months)"*~italic(T[ij]))) + 
    theme(text=element_text(family=chosenFont),legend.position=c(0.785,0.2), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
    scale_colour_manual(name=bquote("Resolution Type"*~italic(R[ij])), values=col.v2, labels=label.v) + 
    scale_fill_manual(name=bquote("Resolution Type"*~italic(R[ij])), values=col.v, labels=label.v) + 
    scale_linetype_manual(name=bquote("Resolution Type"*~italic(R[ij])), values=c("solid","dashed", "dotted"), labels=label.v) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Densities of competing risks overlayed
(g2_Densities_Resol_Type2 <- ggplot(datSurv2[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type2)) + theme_bw() +
    geom_histogram(aes(y=after_stat(density), colour=Resol_Type2, fill=Resol_Type2), position="identity",
                   alpha=0.75, size=0.2) + 
    geom_density(aes(colour=Resol_Type2, linetype=Resol_Type2), size=0.6) + 
    # facets & scale options
    labs(y="", x="", title=paste0("Competing risks (", round(Resol_Type.props[3]*100, digits=1), "%)")) + 
    theme(legend.position=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
          #specific for plot-in-plot
          axis.text.y=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.text.x=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.ticks=element_blank(), axis.title.x=element_blank(), #axis.title.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(color="black", fill="white"),
          plot.background = element_rect(color="white"), plot.margin = unit(c(0,0,0,0),"mm"),
          plot.title = element_text(hjust=0.55,vjust=-10,margin=margin(t=-12))) +
    scale_colour_manual(name=bquote("Resolution Type"*~italic(R[ij])), values=col.v[4:6], labels=label.v2) + 
    scale_fill_manual(name=bquote("Resolution Type"*~italic(R[ij])), values=col.v2[4:6], labels=label.v2) + 
    scale_linetype_manual(name=bquote("Resolution Type"*~italic(R[ij])), values=c("solid","dashed", "dotted"), labels=label.v2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Combining the two above plots onto a single graph
ymin <- diff(ggplot_build(g1_Densities_Resol_Type)$layout$panel_params[[1]]$y.range) * 0.3
ymax <- max(ggplot_build(g1_Densities_Resol_Type)$layout$panel_params[[1]]$y.range) * 0.95
(plot.full <- g1_Densities_Resol_Type + annotation_custom(grob = ggplotGrob(g2_Densities_Resol_Type2), xmin=100, xmax=500, ymin=ymin, ymax=ymax))

# - Save plot
ggsave(plot.full, file=paste0(genFigPath,"Default-FailureTime-Densities.png"),width=1350/dpi, height=1000/dpi,dpi=dpi, bg="white")




# ------ 3. Densities of a dataset including competing risks vs densities of a dataset excluding competing risks
# Densities of a dataset including competing risks
# Creating a new dataset with the competing events grouped together
datSurv3 <- datSurv %>% subset()
datSurv3[, Resol_Type3 := case_when(PerfSpellResol_Type_Hist %in% c("Settled","Paid-up","Written-off","Censored") ~ "Right-Censored",        
                                    TRUE ~ "Default")]



# Calculate resolution proportions for labelling purposes
(Resol_Type.props3 <- table(datSurv3$Resol_Type3) %>% proportions())
# 19.10724% of accounts defaulted
# 80.598% of accounts did not default


# Graphing parameters
chosenFont <- "Cambria"
col.v3 <- brewer.pal(12, "Paired")[c(10,4)]
col.v4 <- brewer.pal(12, "Paired")[c(9,3)]
label.v3 <- c(paste0("Default (", round(Resol_Type.props3[1]*100), "%)"),
              paste0("Right-censored (", round(Resol_Type.props3[2]*100), "%)"))

# - Densities of a dataset including competing risks
(g3_Densities_Resol_Type3 <- ggplot(datSurv3[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type3)) + theme_minimal() + 
    geom_histogram(aes(y=..density.., colour=Resol_Type3, fill=Resol_Type3), position="identity",
                   alpha=0.75, size=0.2) + 
    geom_density(aes(colour=Resol_Type3, linetype=Resol_Type3), size=0.8) + 
    # facets & scale options
    labs(y=bquote(plain(Empirical~failure*' time histogram & density ')~italic(f(t))), 
         x=bquote(Discrete~time~italic(t)*" (months) in performing spell")) + 
    theme(text=element_text(family=chosenFont),legend.position=c(0.785,0.2), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
    scale_colour_manual(name="Resolution Type", values=col.v3, labels=label.v3) + 
    scale_fill_manual(name="Resolution Type", values=col.v4, labels=label.v3) + 
    scale_linetype_manual(name="Resolution Type", values=c("solid","dashed"), labels=label.v3) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)


# --- Densities of a dataset including competing risks
# Creating a new dataset with the competing events completely removed
datSurv4 <- datSurv %>% subset(Resol_Type %in% c("a_Default", "b_Right-censored"))


# Calculate resolution proportions for labelling purposes
(Resol_Type.props4 <- table(datSurv4$Resol_Type) %>% proportions())
# 35.560% of accounts defaulted
# 64.444% of accounts did not default

# Graphing parameters
label.v4 <- c(paste0("Default (", round(Resol_Type.props4[1]*100), "%)"),
              paste0("Right-censored (", round(Resol_Type.props4[2]*100), "%)"))

# - Densities of a dataset excluding competing risks
(g4_Densities_Resol_Type4 <- ggplot(datSurv4[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type)) + theme_minimal() + 
    geom_histogram(aes(y=..density.., colour=Resol_Type, fill=Resol_Type), position="identity",
                   alpha=0.75, size=0.2) + 
    geom_density(aes(colour=Resol_Type, linetype=Resol_Type), size=0.8) + 
    # facets & scale options
    labs(y=bquote(plain(Empirical~failure*' time histogram & density ')~italic(f(t))), 
         x=bquote(Discrete~time~italic(t)*" (months) in performing spell")) + 
    theme(text=element_text(family=chosenFont),legend.position=c(0.785,0.2), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
    scale_colour_manual(name="Resolution Type", values=col.v3, labels=label.v4) + 
    scale_fill_manual(name="Resolution Type", values=col.v4, labels=label.v4) + 
    scale_linetype_manual(name="Resolution Type", values=c("solid","dashed"), labels=label.v4) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)




# ------ 4. Distributional Analysis (Quantitative)
# Testing if the default distribution of the dataset including competing risks is the same as the default distribution of the dataset excluding competing risks.
# NOTE: For two distributions to not be equal they need to differ in at least one aspect of their distributions, e.g. E[X], E[X^2], Var(X), or even their expected future lifetime.

# --- Defaulted Accounts
# - Kolmogorov-Smirnov Test
# NOTE: The two-sample KS test uses the maximum discrepency between two distributions as the test statistic - it therefore doesn't take into account the entire distribution (which is one of its drawbacks).
(test1 <- ks.test(datSurv3[PerfSpellResol_Type_Hist=="Defaulted" & PerfSpell_Age<=300, PerfSpell_Age], datSurv4[PerfSpellResol_Type_Hist=="Defaulted" & PerfSpell_Age<=300, PerfSpell_Age]))
# P-value = 1; Null-hypothesis holds and therefore the two distributions are the same.

# - # Computing the MAEP
#   Computing the Kernel-Density estimates of the right-censored accounts in the two datasets
default_density_incl_cmpr <- density(datSurv3[PerfSpellResol_Type_Hist=="Defaulted" & PerfSpell_Age<=300, PerfSpell_Age])
default_density_excl_cmpr <- density(datSurv4[PerfSpellResol_Type_Hist=="Defaulted" & PerfSpell_Age<=300, PerfSpell_Age])                            

# Computing the MAEP
(default_maep_cmpr <- sum(abs(default_density_excl_cmpr$y-default_density_incl_cmpr$y)/(length(default_density_excl_cmpr$y)/2))/max(max(default_density_excl_cmpr$y),max(default_density_incl_cmpr$y)))
# 0% difference in the densities


# --- (Right-) Censored Accounts
# - Kolmogorov-Smirnov Test]
(test2 <- ks.test(datSurv3[Resol_Type3=="Right-Censored" & PerfSpell_Age<=300, PerfSpell_Age], datSurv4[Resol_Type=="Right-Censored" & PerfSpell_Age<=300, PerfSpell_Age]))
# P-value = 0; Null-hypothesis does not hold and therefore the two distributions are not the same.

# - # Computing the MAEP
#   Computing the Kernel-Density estimates of the right-censored accounts in the two datasets
censor_density_incl_cmpr <- density(datSurv3[PerfSpellResol_Type_Hist!="Defaulted" & PerfSpell_Age<=300, PerfSpell_Age])
censor_density_excl_cmpr <- density(datSurv4[PerfSpellResol_Type_Hist!="Defaulted" & PerfSpell_Age<=300, PerfSpell_Age])                            

# Computing the MAEP
(censor_maep_cmpr <- sum(abs(censor_density_excl_cmpr$y-censor_density_incl_cmpr$y)/(length(censor_density_excl_cmpr$y)))/max(max(censor_density_excl_cmpr$y),max(censor_density_incl_cmpr$y)))
# 4.4455% difference in the densities

# --- Two density plots side-by-side (Competing Risks Incl. vs Excls)
(g5_Densities <- ggarrange(g3_Densities_Resol_Type3,g4_Densities_Resol_Type4))
# - save plot
dpi <- 140
ggsave(g5_Densities, file=paste0(genFigPath,"Sie_by_Side_Comparison_of_Densities.png"),width=1200/dpi, height=1000/dpi,dpi=dpi)




# --- Overlaying the distributions of right-censoring (including vs excluding competing events)
# - Creating a variable to use for grouping
datSurv3 <- datSurv3 %>% mutate(Resol_Type5="Incl_CMPR")
datSurv4 <- datSurv4 %>% mutate(Resol_Type5="Excl_CMPR")

# - Combining the dataset including the competing risks and the datast excluding the competing risks
datSurv5 <- rbind(subset(datSurv3[Resol_Type3=="Right-Censored"], select = c(PerfSpell_Age, Resol_Type5)),subset(datSurv4[Resol_Type=="Right-Censored"], select = c(PerfSpell_Age, Resol_Type5)))

# Graphing parameters
label.v5 <- c(paste0("Right-Censored Excl. CMPR (", round(Resol_Type.props4[2]*100), "%)"), paste0("Right-Censored Incl. CMPR (", round(Resol_Type.props3[2]*100), "%)"))
col.v5 <- brewer.pal(12,"Paired")

(g6_Densities_Resol_Type5 <- ggplot(datSurv5[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type5)) + theme_minimal() + 
    geom_histogram(aes(y=..density.., colour=Resol_Type5, fill=Resol_Type5), position="identity",
                   alpha=0.75, size=0.2) + 
    geom_density(aes(colour=Resol_Type5, linetype=Resol_Type5), size=0.8) + 
    # facets & scale options
    labs(y=bquote(plain(Empirical~failure*' time histogram & density ')~italic(f(t))), 
         x=bquote(Discrete~time~italic(t)*" (months) in performing spell")) + 
    theme(text=element_text(family=chosenFont),legend.position=c(0.785,0.2), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) +
    annotate("text", x = 200, y = 0.0075, label = paste0("KS Test P-Value = ", test2$p.value),
             family=chosenFont, size=5, parse=F) +
    annotate("text", x = 200, y = 0.0070, label = paste0("MAEP = ", round(censor_maep_cmpr,4)),
             family=chosenFont, size=5, parse=F) +
    scale_colour_manual(name="Resolution Type", values=col.v5[c(2,8)], labels=label.v5) + 
    scale_fill_manual(name="Resolution Type", values=col.v5[c(1,7)], labels=label.v5) + 
    scale_linetype_manual(name="Resolution Type", values=c("solid","dashed"), labels=label.v5) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)



# ------ 5. Non-Parametric Estimations
# - Creating a factor from PerfSpellResol_Type_Hist to ensure compatibility with the AJ-estimator
datSurv[, PerfSpellResol_Factor := factor(PerfSpellResol_Type_Hist, labels=c("Censored", "Settled", "Paid-up", "Default", "Written-off"))]


# --- KM-Estimator for Default
kmPerf_default <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, event=PerfSpellResol_Factor=="Default",
                               type="counting") ~ 1, id=PerfSpell_Key, data=datSurv)
summary(kmPerf_default)$table

# --- AJ-Estimator
ajPerf <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, event=PerfSpellResol_Factor,
                       type="counting") ~ 1, id=PerfSpell_Key, data=datSurv)
summary(ajPerf)$table

# --- KM- vs AJ Estimator
# - Creating the required datasets for graphing:
# KM Estimator Dataset
haz_dat_def <- data.table(Time=kmPerf_default$time, AtRisk_n=kmPerf_default$n.risk, 
                          Event_n = kmPerf_default$n.event, Censored_n=kmPerf_default$n.censor,
                          hazard=kmPerf_default$n.event/kmPerf_default$n.risk, 
                          CumulHazard = kmPerf_default$cumhaz,
                          Surv = kmPerf_default$surv,
                          CumulIncidence = 1-kmPerf_default$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()

# AJ Estimator Dataset
haz_dat_AFTs <- data.table(Time=ajPerf$time,
                           AtRisk_n=ajPerf$n.risk[,1],
                           Event_Set= ajPerf$n.event[,2], Event_PaidUp = ajPerf$n.event[,3],
                           Event_Def = ajPerf$n.event[,4], Event_WOff = ajPerf$n.event[,5],
                           Censored_n=ajPerf$n.censor,
                           hazard_Set=ajPerf$n.event[,4]/ajPerf$n.risk[,1],
                           hazard_PaidUp=ajPerf$n.event[,3]/ajPerf$n.risk[,1],
                           hazard_Def=ajPerf$n.event[,2]/ajPerf$n.risk[,1],
                           hazard_WOff=ajPerf$n.event[,5]/ajPerf$n.risk[,1],
                           CumulHazard_Set = ajPerf$cumhaz[,3],
                           CumulHazard_PaidUp = ajPerf$cumhaz[,2],
                           CumulHazard_Def = ajPerf$cumhaz[,1],
                           CumulHazard_WOff = ajPerf$cumhaz[,4],
                           Surv = ajPerf$pstate[,1],
                           CumulIncidence_Set = ajPerf$pstate[,4],
                           CumulIncidence_PaidUp = ajPerf$pstate[,3],
                           CumulIncidence_Def = ajPerf$pstate[,2],
                           CumulIncidence_WOff = ajPerf$pstate[,5])

# Getting the data in the correct format for easier plotting
haz_dat_AFTs <-pivot_longer(haz_dat_AFTs, cols=-c(1:2,7,16), names_sep = "_", names_to = c("Col_Name", "Resol_Type")) %>%  mutate(Col_Name=ifelse(Col_Name=="", "value", Col_Name))
haz_dat_AFTs <-pivot_wider(haz_dat_AFTs, id_cols = c(Time, AtRisk_n, Censored_n, Surv, Resol_Type), names_from = Col_Name, values_from = value)
haz_dat_AFTs <- data.table(haz_dat_AFTs)

# Combining KM- and AJ Datsets
npe_plot <- rbind(subset(haz_dat_def[, Estimator:="KM"], select = c(Time, Surv, CumulIncidence, Estimator)),
                  subset(haz_dat_AFTs[, Estimator:="AJ"][Resol_Type=="Def",], select = c(Time, Surv, CumulIncidence, Estimator)))
# Adding a variable to indicate the difference between the two estimators over time
datDiff <- data.frame(Time = seq(1:length(haz_dat_def$Time)),
                      Surv = abs(haz_dat_def$Surv - haz_dat_AFTs[Resol_Type=="Def",Surv]),
                      CumulIncidence = 0,
                      Estimator = rep("Diff",length(haz_dat_def$Time)))
# Combining datasets and ensuring that the columns are of the correct type
npe_plot <- rbind(npe_plot,datDiff); npe_plot$Time <- as.numeric(npe_plot$Time); npe_plot$Surv <- as.numeric(npe_plot$Surv)

# Graphing parameters
label.v6 <- c("Kaplan-Meier", "Aalen-Johansen", "Diffirence")
col.v6 <- brewer.pal(11,"Spectral")

# Annotation - Mean Absolute Error between the KM and AJ estimators
mae_KM_AJ <- sum(datDiff$Surv)/length(datDiff$Surv)

# Plot - Overall survival probabilities
(g7_km_vs_AJ_surv <- ggplot(data=npe_plot[Time<=500,], aes(x=Time)) +
    geom_line(data=subset(npe_plot, Estimator %in% c("KM", "AJ")), aes(y=Surv, colour=Estimator)) +
    geom_ribbon(data=subset(npe_plot, Estimator=="Diff"), aes(ymin=0, ymax=Surv, fill=Estimator), alpha=0.5) +
    theme_minimal() + 
    labs(x="Performance spell age (months)", y="Survival Probability") + 
    theme(legend.position = "bottom", axis.text.x=element_text(angle=90)) +
    annotate("text",x=100, y=0.075, label =paste0("MAE = ", round(mae_KM_AJ,5)), size=5) +
    scale_colour_manual(name="Estimator", values=col.v6[c(1,3,9)], labels=label.v6) + 
    scale_fill_manual(name="Estimator", values=col.v6[c(1,3,9)], labels=label.v6))

# Annotation
mae_KM_AJ2 <- sum(abs(haz_dat_AFTs[Resol_Type=="Def", CumulIncidence]-haz_dat_def[, CumulIncidence]))/length(datDiff$Surv)

# Plot - Cumulative Incidence Function of Default
(g8_km_vs_AJ_cum_inc <- ggplot(data=npe_plot[Time<=500,], aes(x=Time)) +
    geom_line(data=subset(npe_plot, Estimator %in% c("KM", "AJ")), aes(y=CumulIncidence, colour=Estimator)) +
    geom_ribbon(data=subset(npe_plot, Estimator=="Diff"), aes(ymin=0, ymax=CumulIncidence, fill=Estimator)) +
    theme_minimal() + 
    labs(x="Performance spell age (months)", y="Survival Probability") + 
    theme(legend.position = "bottom", axis.text.x=element_text(angle=90)) +
    annotate("text",x=400, y=0.9, label =paste0("MAE = ", round(mae_KM_AJ2,5)), size=5) +
    scale_colour_manual(name="Estimator", values=col.v6[c(1,3,9)], labels=label.v6) + 
    scale_fill_manual(name="Estimator", values=col.v6[c(1,3,9)], labels=label.v6))



