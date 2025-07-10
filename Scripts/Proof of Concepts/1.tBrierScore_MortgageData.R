# ========================================= TIME-DEPENDENT BRIER SCORES ================================
# Calculate and compare the time-dependent Brier scores across survival models
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
#   - 2f.Data_Fusion1.R
#   - 3b.Data_Fusion2_PWP_ST.R
#   - 5a(i).CoxPropHaz_PWP_Advanced.R
#   - 5a(i).CoxPropHaz_PWP_Basic.R
#   - 5b(i).CoxDiscreteTime_Advanced.R
#   - 5b(ii).CoxDiscreteTime_Basic.R
#   - 6d.Analytics_tBrierScore.R

# -- Inputs:
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b
#
# -- Outputs:
#   - <Analytics> | Graphs
# =====================================================================================================
# Note that this script is built within the ambit of a larger project "Term-Structure-Modelling-Retailmortgages"
# ------------------------------------------------------------------------------------------------------




# ----------------- 1. Load & prepare data for analysis

# - Pathing variables
# genPath <- "C:/Data/" # Set this to specific data folder and uncomment
# path_cust <- "C:/Data/" # set this to where all scripts are hosted and uncomment
# genFigPath <- "C:/Data/" # set this to where all figures should live and uncomment

# - Compile custom functions
source(paste0(path_cust,'0e.FunkySurv_tBrierScore.R'))

# - Load necessary packages
require(data.table)
require(survival)
require(ffbase); tempPath <- "C:/TempData"; options("fftempdir"=tempPath)
require(RColorBrewer)
require(ggplot2)
require(extrafont) #remotes::install_version("Rttf2pt1", version = "1.3.8"); Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.55.0/bin/gswin32c.exe"); font_import(); loadfonts(); loadfonts(device="win")
require(scales)

# - Confirm prepared dataset is loaded into memory
unpack.ffdf(paste0(genPath,"datCredit"), tempPath);gc()
datCredit_train <- subset(datCredit, Sample == "Training")
datCredit_validation <- subset(datCredit, Sample == "Validation")





# ----------------- 2. Fit a discrete-time hazard model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic discrete-time hazard model
# - Initialize variables
vars_basic <- c("-1", "Time_Binned", "log(TimeInPerfSpell):PerfSpell_Num_binned",
                "Arrears", "InterestRate_Nom", "M_Inflation_Growth_6")

# - Fit discrete-time hazard model with selected variables
modLR_basic <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars_basic, collapse = " + "))),
                    data=datCredit_train, family="binomial", weights = Weight)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced discrete-time hazard model
# - Initialize variables
vars <- c("-1", "Time_Binned*PerfSpell_Num_binned", #"log(TimeInPerfSpell):PerfSpell_Num_binned",
          "g0_Delinq_SD_4", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom",
          "M_Inflation_Growth_6","M_DTI_Growth")

# - Fit discrete-time hazard model with selected variables
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial", weights = Weight)




# ----------------- 3. Calculate time-dependent Brier scores across survival models

# --- Prentice-Williams-Peterson (PWP) Spell-time definition | Basic discrete-time hazard model

objCoxDisc_bas <- tBrierScore(datCredit, modGiven=modLR_basic, predType="response", spellPeriodMax=300, fldKey="PerfSpell_Key", 
                              fldStart="Start", fldStop="TimeInPerfSpell",fldCensored="PerfSpell_Censored", 
                              fldSpellAge="PerfSpell_Age", fldSpellOutcome="PerfSpellResol_Type_Hist")

# --- Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced discrete-time hazard model

objCoxDisc_adv <- tBrierScore(datCredit, modGiven=modLR, predType="response", spellPeriodMax=300, fldKey="PerfSpell_Key", 
                              fldStart="Start", fldStop="TimeInPerfSpell",fldCensored="PerfSpell_Censored", 
                              fldSpellAge="PerfSpell_Age", fldSpellOutcome="PerfSpellResol_Type_Hist")





# ----------------- 4. Graph tBS-values across models
# --- Discrete-time hazard models

# - Data fusion across models
datGraph <- rbind(data.table(objCoxDisc_bas$tBS, Type="a_Basic"),
                  data.table(objCoxDisc_adv$tBS, Type="b_Advanced"))

# - Aesthetic engineering
datGraph[, FacetLabel := "Discrete-time hazard models"]
zoomedSpellAge <- 120

# - Recalculate Integrated tBS over zoomed spell age
ibs_bas <- mean(objCoxDisc_bas$tBS[TimeInPerfSpell <= zoomedSpellAge, Brier])
ibs_adv <- mean(objCoxDisc_adv$tBS[TimeInPerfSpell <= zoomedSpellAge, Brier])

# - Graphing Parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(8, "Set2")[c(2,1)]
vLabel <- c("a_Basic"="Basic", "b_Advanced"="Advanced")

# - Main graph of tBS
(gOuter <- ggplot(datGraph, aes(x=TimeInPerfSpell, y=Brier, group=Type)) + theme_minimal() + 
    labs(y=bquote("Time-dependent Brier Score "*italic(B)[s](italic(t))), x=bquote("Performing spell age (months)"*~italic(t))) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
    # Main graph
    geom_line(aes(colour=Type, linetype=Type), linewidth=0.5) + 
    # Annotations
    annotate(geom="text", x=50, y=20, label=paste0("IBS (Basic): ", round(objCoxDisc_bas$IBS,3)), 
             family=chosenFont, size=3.5, colour=vCol[1]) + 
    annotate(geom="text", x=43, y=13, label=paste0("IBS (Advanced): ", round(objCoxDisc_adv$IBS,3)), 
             family=chosenFont, size=3.5, colour=vCol[2]) +     
    # Facets & scales
    facet_grid(FacetLabel ~ .) +  
    scale_colour_manual(name="Model", values=vCol, labels=vLabel) + 
    scale_linetype_discrete(name="Model", labels=vLabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Zoomed inset graph of tBS on a smaller time scale
(gInner <- ggplot(datGraph[TimeInPerfSpell<=zoomedSpellAge,], aes(x=TimeInPerfSpell, y=Brier, group=Type)) + 
    theme_bw() + labs(y="", x="") + 
    theme(legend.position.inside=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
          #specific for plot-in-plot
          axis.text.y=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.text.x=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.ticks=element_blank(), axis.title.x=element_blank(), #axis.title.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(color="black", fill="white"),
          plot.background = element_rect(color="white"), plot.margin = unit(c(0,0,0,0),"mm"),
          plot.title = element_text(hjust=0.55,vjust=-10,margin=margin(t=-12))) + 
    # Main graph
    geom_line(aes(colour=Type, linetype=Type), linewidth=0.5, show.legend = F) + 
    # Annotations
    annotate(geom="text", x=50, y=0.85, label=paste0("IBS (Basic): ", round(ibs_bas,3)), 
             family=chosenFont, size=3.5, colour=vCol[1]) + 
    annotate(geom="text", x=45, y=0.7, label=paste0("IBS (Advanced): ", round(ibs_adv,3)), 
             family=chosenFont, size=3.5, colour=vCol[2]) +   
    # Facets & scales
    scale_colour_manual(name="", values=vCol, labels=vLabel) + 
    scale_linetype_discrete(name="", labels=vLabel) +
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Combining the two above plots onto a single graph
(plot.full <- gOuter + annotation_custom(grob = ggplotGrob(gInner), xmin=0, xmax=210, ymin=30, ymax=100))

# - Save plot
dpi <- 280
ggsave(plot.full, file=paste0(genFigPath,"tBrierScores_CoxDisc.png"),width=1350/dpi, height=1000/dpi,dpi=dpi, bg="white")


