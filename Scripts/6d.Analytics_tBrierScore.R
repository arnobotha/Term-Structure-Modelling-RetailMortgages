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

# -- Inputs:
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b
#
# -- Outputs:
#   - <Analytics> | Graphs
# ------------------------------------------------------------------------------------------------------






# ----------------- 1. Load & prepare data for analysis

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# - Use only performance spells
datCredit_train <- datCredit_train_PWPST[!is.na(PerfSpell_Num),]
datCredit_valid <- datCredit_valid_PWPST[!is.na(PerfSpell_Num),]

# Create start and stop columns
datCredit_train[, Start := TimeInPerfSpell - 1]
datCredit_valid[, Start := TimeInPerfSpell - 1]

# - Weigh default cases heavier. as determined interactively based on calibration success (script 6e)
datCredit_train[, Weight := ifelse(DefaultStatus1==1,10,1)]
datCredit_valid[, Weight := ifelse(DefaultStatus1==1,10,1)] # for merging purposes




# ----------------- 2a. Fit a Cox regression model on the resampled prepared data

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic Cox-regression model
# - Initialize variables
vecVars_PWPST_bas <- c("Arrears", "InterestRate_Nom")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_basic <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                           paste(vecVars_PWPST_bas,collapse=" + "), 
                                           " + strata(PerfSpell_Num_binned)")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T, x=T)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced Cox-regression model
# - Initialize variables
vecVars_PWPST_adv <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_adv <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                         paste(vecVars_PWPST_adv,collapse=" + "), 
                                         " + strata(PerfSpell_Num_binned)")),
                       id=PerfSpell_Key, datCredit_train, ties="efron", model=T, x=T)




# ----------------- 2b. Fit a discrete-time hazard model on the resampled prepared data


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

# - Create pointer to the appropriate data object 
datCredit <- rbind(datCredit_train, datCredit_valid)


# --- Prentice-Williams-Peterson (PWP) Spell-time definition | Basic Cox-regression model

objCoxPH_bas <- tBrierScore(datCredit, modGiven=cox_PWPST_basic, predType="exp", spellPeriodMax=300, fldKey="PerfSpell_Key", 
                              fldStart="Start", fldStop="TimeInPerfSpell",fldCensored="PerfSpell_Censored", 
                              fldSpellAge="PerfSpell_Age", fldSpellOutcome="PerfSpellResol_Type_Hist")

# --- Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced Cox-regression model

objCoxPH_adv <- tBrierScore(datCredit, modGiven=cox_PWPST_adv, predType="exp", spellPeriodMax=300, fldKey="PerfSpell_Key", 
                            fldStart="Start", fldStop="TimeInPerfSpell",fldCensored="PerfSpell_Censored", 
                            fldSpellAge="PerfSpell_Age", fldSpellOutcome="PerfSpellResol_Type_Hist")

# --- Prentice-Williams-Peterson (PWP) Spell-time definition | Basic discrete-time hazard model

objCoxDisc_bas <- tBrierScore(datCredit, modGiven=modLR_basic, predType="response", spellPeriodMax=300, fldKey="PerfSpell_Key", 
                              fldStart="Start", fldStop="TimeInPerfSpell",fldCensored="PerfSpell_Censored", 
                              fldSpellAge="PerfSpell_Age", fldSpellOutcome="PerfSpellResol_Type_Hist")

# --- Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced discrete-time hazard model

objCoxDisc_adv <- tBrierScore(datCredit, modGiven=modLR, predType="response", spellPeriodMax=300, fldKey="PerfSpell_Key", 
                              fldStart="Start", fldStop="TimeInPerfSpell",fldCensored="PerfSpell_Censored", 
                              fldSpellAge="PerfSpell_Age", fldSpellOutcome="PerfSpellResol_Type_Hist")




# ----------------- 4. Graph tBS-values across models

# --- Cox Proportional Hazard models

# - Data fusion across models
datGraph <- rbind(data.table(objCoxPH_bas$tBS, Type="a_Basic"),
                  data.table(objCoxPH_adv$tBS, Type="b_Advanced"))

# - Aesthetic engineering
datGraph[, FacetLabel := "Cox Proportional Hazard (CPH) models"]
zoomedSpellAge <- 120

# - Recalculate Integrated tBS over zoomed spell age
ibs_bas <- mean(objCoxPH_bas$tBS[TimeInPerfSpell <= zoomedSpellAge, Brier])
ibs_adv <- mean(objCoxPH_adv$tBS[TimeInPerfSpell <= zoomedSpellAge, Brier])

# - Graphing Parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(8, "Set2")[c(2,1)]
vLabel <- c("a_Basic"="Basic", "b_Advanced"="Advanced")

# - Main graph of tBS
(gOuter <- ggplot(datGraph, aes(x=TimeInPerfSpell, y=Brier, group=Type)) + theme_minimal() + 
    labs(y=bquote("Time-dependent Brier Score "*italic(B)[s](italic(t))), x=bquote("Performing spell age (months) "*~italic(t))) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) + 
    # Main graph
    geom_line(aes(colour=Type, linetype=Type), linewidth=0.5) + 
    # Annotations
    annotate(geom="text", x=50, y=0.15, label=paste0("IBS (Basic): ", round(objCoxPH_bas$IBS,3)), 
             family=chosenFont, size=3.5, colour=vCol[1]) + 
    annotate(geom="text", x=43, y=0.1, label=paste0("IBS (Advanced): ", round(objCoxPH_adv$IBS,3)), 
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
    theme(legend.position=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
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
    annotate(geom="text", x=60, y=0.01, label=paste0("IBS (Basic): ", round(ibs_bas,3)), 
             family=chosenFont, size=3.5, colour=vCol[1]) + 
    annotate(geom="text", x=55, y=0.009, label=paste0("IBS (Advanced): ", round(ibs_adv,3)), 
             family=chosenFont, size=3.5, colour=vCol[2]) +   
    # Facets & scales
    scale_colour_manual(name="", values=vCol, labels=vLabel) + 
    scale_linetype_discrete(name="", labels=vLabel) +
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Combining the two above plots onto a single graph
(plot.full <- gOuter + annotation_custom(grob = ggplotGrob(gInner), xmin=0, xmax=200, ymin=0.25, ymax=0.9))

# - Save plot
dpi <- 280
ggsave(plot.full, file=paste0(genFigPath,"tBrierScores_CoxPH.png"),width=1350/dpi, height=1000/dpi,dpi=dpi, bg="white")



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
  annotate(geom="text", x=50, y=18, label=paste0("IBS (Basic): ", round(objCoxDisc_bas$IBS,3)), 
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
  theme(legend.position=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
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
  annotate(geom="text", x=50, y=0.8, label=paste0("IBS (Basic): ", round(ibs_bas,3)), 
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


# -- Cleanup
rm(gInner, datGraph, gOuter, objCoxDisc_adv, objCoxDisc_bas, objCoxPH_adv, objCoxPH_bas, plot.full)



##### TO-DO
# - Compare tBS for Cox PH model with that obtained from pec(). And debug!

library(riskRegression)

y <- Score(list("Cox.bas"=cox_PWPST_basic),
        formula=as.formula(paste0("Hist(TimeInPerfSpell,DefaultStatus1) ~ 1")),
        data=datCredit_train, conf.int=FALSE, summary=c("risks","IPA","ibs"),
        times=1:300,ncpus=4, parallel="multicore") 


library(pec)

test <- pec(
  object = list("Cox.bas" = cox_PWPST_basic),
  formula = Surv(TimeInPerfSpell, PerfSpell_Event) ~ 1,
  data = datCredit_train,
  times = 1:300,
  exact = T,
  cens.model = "marginal"
)
### REESULTS: This does not work. Need to apply mind.