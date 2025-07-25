# ================================= PERFORMANCE SPELL ANALYSES =========================================
# Analysing the performance spells in various ways, including the maximum spell number (histogram),
# the number of risk events over time per spell number, the failure time histogra/densities per 
# resolution type, and the extend of tied events (Proportion tied relative to total event time frequencies)
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
#   - Bar chart of LoanIDs by performance spells
#   - Bar chart of risk events by performance spells
#   - Histogram of failure times f(t)
# ------------------------------------------------------------------------------------------------------


# ----------------- 1. Max spell number ---
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)

# lookup
lookup <- subset(datCredit_real, LoanID == datCredit_real[PerfSpell_Num == 8 & PerfSpell_Counter == 1, LoanID][1])

# - Aggregation to account-level
# NOTE: Assign max conditionally since there are loans that are forever in default and hence will have no 
# information on performing spells
datAggr <- datCredit_real[, list(MaxPerfNum = ifelse(all(is.na(PerfSpell_Num)), 0, 
                                   max(PerfSpell_Num, na.rm=T)) ), by=list(LoanID)]

# - Analysis on Maximum performing spell number
describe(datAggr$MaxPerfNum); hist(datAggr$MaxPerfNum)
### RESULTS: Mean of 1.091 max spells (median: 1), with 5%-95% at [1, 2]. Large outliers of up to 10 spells

# - Rebin
datAggr[, MaxPerfNum_Binned := ifelse(MaxPerfNum >= 5, 5, MaxPerfNum)]
describe(datAggr$MaxPerfNum_Binned); hist(datAggr$MaxPerfNum_Binned)
### REUSLTS: Mean of 1.09 max spells. 1 spell: 92.6%; 2 spells: 5%; 3 spells: 1.2%; 4 Spells: 0.4%; 5+ spells: 0.2%
# Note: rebinning at 4 cap also tried, though this resulted in a mean of 1.89, which is too far from 'true' mean of 0.92,
# at least anecdotally. However, when subsampling, this binning scheme may need to be revisited to allow feasible sample sizes

# - Aesthetic engineering
chosenFont <- "Cambria"; colPalette <- "BrBG"
datAggr[MaxPerfNum_Binned > 0, MaxPerfNum_Binned_Total := .N]
totFreq <- datAggr[MaxPerfNum_Binned > 0, .N]
datAggr2 <- unique(datAggr[MaxPerfNum_Binned > 0, list(MaxPerfNum_Binned_Pc = .N / MaxPerfNum_Binned_Total,
                                  MaxPerfNum_Binned_Freq = .N), by=list(MaxPerfNum_Binned)])
datAggr2[, MaxPerfNum_Binned_Pc_labelY := MaxPerfNum_Binned_Freq + totFreq*0.005]
vLabelX <- c("1"="1", "2"="2", "3"="3", "4"="4", "5"="5+")
vBreaks <- 1:length(vLabelX)
vCol <- brewer.pal(10, "Paired")

# - Graph
(g1 <- ggplot(datAggr[MaxPerfNum_Binned > 0,], aes(x=MaxPerfNum_Binned)) + theme_minimal() + 
  theme(text=element_text(family=chosenFont)) + 
  labs(y="Frequency", x="Maximum Number of Performance Spells (pre-binned)") + 
  geom_bar(fill = vCol[2]) +
  geom_label(data=datAggr2, aes(y=MaxPerfNum_Binned_Pc_labelY, label=paste(percent(MaxPerfNum_Binned_Pc, accuracy=0.1))), 
             family=chosenFont, fill = vCol[1]) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(labels=vLabelX, breaks=vBreaks) )

# - Save graph
dpi <- 200
ggsave(g1, file=paste0(genFigPath, "MaxPerfSpellNum_hist.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")

# House keeping
rm(datAggr,datAggr2,lookup,g1);gc()





# ----------------- 2. Default resolution rate over time per numbered spell | Cohort-end

# --- Merge datasets together for graphing purposes, subset necessary fields, and rename columns for graphing ease
datGraph <- rbind(datCredit_real[PerfSpell_Num==1, list(Date, LoanID, PerfSpell_Key, PerfSpell_Num, 
                                                        PerfSpellResol_Type_Hist, Sample="Spell 1")],
                  datCredit_real[PerfSpell_Num==2, list(Date, LoanID, PerfSpell_Key, PerfSpell_Num, 
                                                        PerfSpellResol_Type_Hist, Sample="Spell 2")],
                  datCredit_real[PerfSpell_Num>=3, list(Date, LoanID, PerfSpell_Key, PerfSpell_Num, 
                                                        PerfSpellResol_Type_Hist, Sample="Spell 3+")])

# - Creating spell-level min/max date variables as stratifiers
datGraph[, PerfSpellDate_End := Date[.N], by=list(Sample, PerfSpell_Key)]

# - Setting some aggregation parameters, purely to facilitate graphing aesthetics
StartDte <- min(datCredit_real$Date, na.rm=T)
EndDte <- max(datCredit_real$Date, na.rm=T)
maxDate <- EndDte %m-% months(1)# A post-hoc filter, used for graphing purposes - left as the end of the sampling window
minDate <- StartDte  #+ month(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window

# - Fixing to spell entry-time, we aggregate to monthly level and observe the time series up to given point
datAggr_cohorts <- merge(datGraph[Date==PerfSpellDate_End, list(Sum_Total = .N), by=list(Sample,Date)],
                         datGraph[Date==PerfSpellDate_End, list(Sum_Resol = .N), by=list(Sample,Date,PerfSpellResol_Type_Hist)],
                         by=c("Sample", "Date"))[Date >= minDate & Date <= maxDate,]
datAggr_cohorts[, Prop := Sum_Resol/Sum_Total]

# - Graphing parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(9, "Set1")

(g1 <- ggplot(datAggr_cohorts[PerfSpellResol_Type_Hist=="Defaulted",], aes(x=Date, y=Prop)) + theme_minimal() + 
    labs(x=bquote("Performing spell cohorts (ccyymm): stop time "*italic(t[s])), y="Default resolution rate (%)") +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Sample, linetype=Sample)) + 
    geom_point(aes(colour=Sample, shape=Sample), size=1) + 
    # scale options
    scale_colour_manual(name="", values=vCol) + 
    scale_shape_discrete(name="") + scale_linetype_discrete(name="") + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))


# - Save graph
dpi <- 220
ggsave(g1, file=paste0(genFigPath, "DefResolRate_SpellNumber.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")

# - Cleanup
rm(datAggr_cohorts, datGraph)




# ----------------- 3. Failure time histogram & densities per resolution type


# ------Preliminaries

# - Select one record per performance spell
datSurv <- subset(datCredit_real, PerfSpell_Counter==1 & PerfSpell_Age < 500,
                  select=c("LoanID","PerfSpell_Key","PerfSpell_Age", "PerfSpellResol_Type_Hist", "PerfSpell_Num", "TimeInPerfSpell")); gc()

# - Group competing risks together in Resol_Type
datSurv[,Resol_Type := case_when(PerfSpellResol_Type_Hist == "Defaulted" ~ "a_Default",
                                 PerfSpellResol_Type_Hist == "Censored" ~ "b_Right-censored",
                                 PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "c_Competing_risk")]
# [SANITY CHECK] Frequency analysis
datSurv$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
(Resol_Type.props <- datSurv$Resol_Type %>% table() %>% prop.table())

# - Create dataset containing only the competing risks, where the type of competing risk is shown by Resol_Type2
datSurv2 <- datSurv %>% subset(PerfSpellResol_Type_Hist%in%c("Settled", "Paid-up", "Written-off"))
datSurv2[, Resol_Type2 := case_when(PerfSpellResol_Type_Hist=="Settled" ~ "a_Settlement",
                                    PerfSpellResol_Type_Hist=="Paid-up" ~ "b_Paid-up",
                                    TRUE ~ "c_Written-off-Other")]
# [SANITY CHECK] Frequency analysis
datSurv$PerfSpellResol_Type_Hist %>% table() %>% prop.table()
(Resol_Type2.props <- datSurv2$Resol_Type2 %>% table() %>% prop.table())


# - Graphing Parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(10, "Paired")[c(10,6,4, 2,1,8)]
vCol2 <- brewer.pal(10, "Paired")[c(10,5,3, 1,2,8)]
vLabels <- c(paste0("a_Default"="Default (", round(Resol_Type.props[1]*100, digits=1), "%)"), # Need to round to the first decimal place to ensure that the prior add up to one
             paste0("b_Right-censored"="Right-censored (", round(Resol_Type.props[2]*100, digits=1), "%)"),
             paste0("c_Competing_risk"="Competing Risks (", round(Resol_Type.props[3]*100, digits=1), "%)"))
vLabels2 <- c(paste0("a_Settlement"="Settlement (", round(Resol_Type2.props[1]*100, digits=1), "%)"),
              paste0("b_Paid-up"="Paid-up (", round(Resol_Type2.props[2]*100, digits=1), "%)"),
              paste0("c_Written-off-Other"="Other/Write-off (", round(Resol_Type2.props[3]*100, digits=1), "%)"))


# - Densities of resolution types overlaid
(g1_Densities_Resol_Type <- ggplot(datSurv[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type)) + theme_minimal() + 
    labs(y=bquote(plain('Empirical histogram & density ')~italic(f(t))), 
         x=bquote("Performing spell age (months)"*~italic(t))) + 
    theme(text=element_text(family=chosenFont),legend.position.inside=c(0.785,0.2), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90),
          legend.position=c(0.75,0.20)) + 
    # Graphs
    geom_histogram(aes(y=after_stat(density), colour=Resol_Type, fill=Resol_Type), position="identity",
                   alpha=0.75, linewidth=0.2) + 
    geom_density(aes(colour=Resol_Type, linetype=Resol_Type), linewidth=0.8) + 
    # facets & scale options
    scale_colour_manual(name=bquote("Resolution Type"*~italic(R)), values=vCol2, labels=vLabels) + 
    scale_fill_manual(name=bquote("Resolution Type"*~italic(R)), values=vCol, labels=vLabels) + 
    scale_linetype_manual(name=bquote("Resolution Type"*~italic(R)), values=c("solid","dashed", "dotted"), labels=vLabels) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Densities of competing risks overlaid
(g2_Densities_Resol_Type2 <- ggplot(datSurv2[PerfSpell_Age<=500,], aes(x=PerfSpell_Age, group=Resol_Type2)) + theme_bw() +
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
    # Graphs
    geom_histogram(aes(y=after_stat(density), colour=Resol_Type2, fill=Resol_Type2), position="identity",
                   alpha=0.75, size=0.2) + 
    geom_density(aes(colour=Resol_Type2, linetype=Resol_Type2), linewidth=0.6) + 
    # facets & scale options
    scale_colour_manual(name="", values=vCol[4:6], labels=vLabels2) + 
    scale_fill_manual(name="", values=vCol2[4:6], labels=vLabels2) + 
    scale_linetype_manual(name="", values=c("solid","dashed", "dotted"), labels=vLabels2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Combining the two above plots onto a single graph
ymin <- diff(ggplot_build(g1_Densities_Resol_Type)$layout$panel_params[[1]]$y.range) * 0.3
ymax <- max(ggplot_build(g1_Densities_Resol_Type)$layout$panel_params[[1]]$y.range) * 0.95
(plot.full <- g1_Densities_Resol_Type + annotation_custom(grob = ggplotGrob(g2_Densities_Resol_Type2), xmin=100, xmax=500, ymin=ymin, ymax=ymax))

# - Save plot
dpi <- 170
ggsave(plot.full, file=paste0(genFigPath,"Default-FailureTime-Densities.png"),width=1350/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Clean-up
rm(datSurv, datSurv2, g1_Densities_Resol_Type, g2_Densities_Resol_Type2, plot.full)





# ----------------- 4. Extent of tied events: Proportion tied relative to total event time frequencies
describe(datCredit_real$PerfSpellResol_Type_Hist)

# - Subset performance spell-level observations towards calculating tied event times
datSpells <- subset(datCredit_real, PerfSpell_Counter==1 & PerfSpell_Age <= 300,
                  select=c("LoanID","PerfSpell_Key","PerfSpell_Age", "PerfSpellResol_Type_Hist", "PerfSpell_Num", "TimeInPerfSpell")); gc()

# - Group competing risks together in Resol_Type
datSpells[,Resol_Type := case_when(PerfSpellResol_Type_Hist == "Defaulted" ~ "a_Default",
                                 PerfSpellResol_Type_Hist == "Censored" ~ "b_Right-censored",
                                 PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "c_Competing_risk")]

# - Apply Latent Risks assumption; competing risks are 'right-censored'
datSpells[,Resol_Type_LR := case_when(PerfSpellResol_Type_Hist == "Defaulted" ~ "a_Default",
                                      TRUE ~ "b_Right-censored")]

# - Aggregate to spell age level and count number of tied events
datAggr <- datSpells[, list(Freq = .N), by=list(PerfSpell_Age, Resol_Type_LR)] %>% setkey(PerfSpell_Age, Resol_Type_LR)

# - Aesthetic engineering
(Resol_Type_LR.props <- datSpells$Resol_Type_LR %>% table() %>% prop.table()) # get proportions of data
(tiedSpells <- datAggr[Resol_Type_LR == "a_Default" & Freq > 1, .N] / datAggr[Resol_Type_LR == "a_Default", .N])

# - Graphing Parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(10, "Paired")[c(10,8, 9, 7)]
vLabels <- c("a_Default"=paste0("Default (", round(Resol_Type_LR.props[1]*100, digits=1), "%)"), # Need to round to the first decimal place to ensure that the prior add up to one
             "b_Right-censored"=paste0("Right-censored & competing risks (", round(Resol_Type_LR.props[2]*100, digits=1), "%)"))

# - Create graph
(g1 <- ggplot(datAggr, aes(x=PerfSpell_Age, y=Freq, group=Resol_Type_LR)) + theme_minimal() +
  labs(y="Tied spell age volumes", 
       x=bquote("Performing spell ages (months)"*~italic(t))) + 
  theme(text=element_text(family=chosenFont),legend.position="inside", 
        legend.position.inside = c(0.65,0.6), 
        strip.background=element_rect(fill="snow2", colour="snow2"),
        strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
  # Main graphs
  geom_line(aes(colour=Resol_Type_LR, linetype=Resol_Type_LR), linewidth=0.5) + 
  geom_point(aes(colour=Resol_Type_LR, shape=Resol_Type_LR), size=1) + 
  annotate(geom="text", x=150, y=5000, label=paste0("% tied spell ages: ", percent(tiedSpells, accuracy=0.01)), 
           family=chosenFont, size=3) + 
  # Facets and scales
  scale_colour_manual(name="Resolution Type", values=vCol, labels=vLabels) + 
  scale_linetype_discrete(name="Resolution Type", labels=vLabels) + 
  scale_shape_discrete(name="Resolution Type", labels=vLabels) + 
  scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
  scale_x_continuous(breaks=breaks_pretty(), label=comma) )

# - Save plot
dpi <- 200
ggsave(g1, file=paste0(genFigPath,"TiedEvents_Extent.png"),width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Cleanup
rm(datAggr, datSpells, g1)






# ----------------- 5. Extent of censored spells over spell age

datAggr <- datCredit_real[, list(CensorRate = mean(ifelse(PerfSpellResol_Type_Hist=="Censored", 1, 0), na.rm=T)), 
                          by=list(PerfSpell_Age)]
datAggr[, Resol_Type := "Right-censoring rate"]

describe(datAggr[PerfSpell_Age<=300, CensorRate]); hist(datAggr[PerfSpell_Age<=300, CensorRate], breaks="FD")

# - Annotations
sMeanCensor <- mean(datAggr[PerfSpell_Age<=300, CensorRate], na.rm=T)

# - Graphing Parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(10, "Paired")[c(2,1)]

(ggOuter <- ggplot(datAggr[PerfSpell_Age<=300, ], aes(x=CensorRate)) + theme_minimal() +
  labs(y="Histogram and empirical density", 
       x="Censoring rate per unique spell age") + 
  theme(text=element_text(family=chosenFont),legend.position="inside", 
        legend.position.inside = c(0.2,0.4), 
        strip.background=element_rect(fill="snow2", colour="snow2"),
        strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
  # Graphs
  geom_histogram(aes(y=after_stat(density), colour=Resol_Type, fill=Resol_Type), position="identity",
                 alpha=0.75, size=0.2, show.legend = F) + 
  geom_density(aes(colour=Resol_Type), linewidth=0.8, show.legend = F) + 
  scale_colour_manual(name="", values=vCol[1]) + 
  scale_fill_manual(name="", values=vCol[2]) + 
  scale_y_continuous(label=comma) + scale_x_continuous(label=percent))


(ggInner <- ggplot(datAggr[PerfSpell_Age<=300,], aes(x=PerfSpell_Age, y=CensorRate)) + theme_bw() +
  labs(y="Mean censoring rate (%)", x="Spell age") + 
  theme(text=element_text(size=12, family="Cambria"),
        #specific for plot-in-plot
        axis.text.y=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
        axis.text.x=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
        axis.ticks=element_blank(), #axis.title.x=element_blank(), #axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill="white"),
        plot.background = element_rect(color="white"), plot.margin = unit(c(0,0,0,0),"mm"),
        plot.title = element_text(hjust=0.55,vjust=-10,margin=margin(t=-12))) + 
  geom_line(aes(colour=Resol_Type), show.legend = F) + 
  geom_point(aes(colour=Resol_Type), show.legend=F, size= 0.5) + 
    geom_hline(aes(colour=Resol_Type), linewidth=0.6, yintercept=sMeanCensor) +
    annotate(geom="text", label=paste0("Mean: ", percent(sMeanCensor,accuracy=0.01)),
             x=75, y=0.6,family=chosenFont, size=4) + 
  scale_colour_manual(name="", values=vCol[1]) + 
  scale_y_continuous(label=percent))


# - Combining the two above plots onto a single graph
ymin <- diff(ggplot_build(ggOuter)$layout$panel_params[[1]]$y.range) * 0.45
ymax <- max(ggplot_build(ggOuter)$layout$panel_params[[1]]$y.range) * 0.975
(plot.full <- ggOuter + annotation_custom(grob = ggplotGrob(ggInner), xmin=0.45, xmax=0.875, ymin=ymin, ymax=ymax))

# - Save plot
dpi <- 170
ggsave(plot.full, file=paste0(genFigPath,"CensoringRate-HistogramDensity.png"),width=1350/dpi, height=1000/dpi,dpi=dpi, bg="white")
