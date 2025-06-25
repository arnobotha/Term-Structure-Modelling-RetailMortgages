# ============================== TERM-STUCTURE: Cox Discrete-time Hazard ===============================
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




# ----------------- 2a. Fit a Cox regression model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic Cox-regression model
# - Initialize variables
vecVars_PWPST_bas <- c("Arrears", "InterestRate_Nom")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_basic <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                           paste(vecVars_PWPST_bas,collapse=" + "), 
                                           " + strata(PerfSpell_Num_binned)")),
                         id=PerfSpell_Key, data=datCredit_train, ties="efron")
summary(cox_PWPST_basic); AIC(cox_PWPST_basic); concordance(cox_PWPST_basic)



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
                       id=PerfSpell_Key, datCredit_train, ties="efron")
summary(cox_PWPST_adv); AIC(cox_PWPST_adv); concordance(cox_PWPST_adv)




# ----------------- 2b. Fit a discrete-time hazard model on the resampled prepared data


# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Basic discrete-time hazard model
# - Initialize variables
vars_basic <- c("-1", "Time_Binned", "log(TimeInPerfSpell):PerfSpell_Num_binned",
                "Arrears", "InterestRate_Nom", "M_Inflation_Growth_6")

# - Fit discrete-time hazard model with selected variables
modLR_basic <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars_basic, collapse = " + "))),
                    data=datCredit_train, family="binomial")
summary(modLR_basic)



# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced discrete-time hazard model
# - Initialize variables
vars <- c("-1", "Time_Binned*PerfSpell_Num_binned", #"log(TimeInPerfSpell):PerfSpell_Num_binned",
          "g0_Delinq_SD_4", "g0_Delinq_Lag_1", "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_mean",
          "AgeToTerm_Aggr_Mean", "InstalmentToBalance_Aggr_Prop", "NewLoans_Aggr_Prop",
          "pmnt_method_grp", "InterestRate_Nom", "BalanceToPrincipal",
          "M_Inflation_Growth_6","M_DTI_Growth")

# - Fit discrete-time hazard model with selected variables
modLR <- glm( as.formula(paste("PerfSpell_Event ~", paste(vars, collapse = " + "))),
              data=datCredit_train, family="binomial")
summary(modLR)




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








tBrierScore <- function(datGiven, modGiven, predType="response", spellPeriodMax=300,
                        fldKey="PerfSpell_Key", fldStart = "Start", fldStop="TimeInPerfSpell",
                        fldCensored="PerfSpell_Censored", fldSpellAge="PerfSpell_Age",
                        fldSpellOutcome="PerfSpellResol_Type_Hist") {
  # Testing conditions
  # datGiven <- datCredit; modGiven <- modLR; predType <- "response"; spellPeriodMax <- 300
  # fldKey <- "PerfSpell_Key"; fldStart <- "Start"; fldStop<-"TimeInPerfSpell";
  # fldCensored<-"PerfSpell_Censored"; fldSpellAge<-"PerfSpell_Age"; fldSpellOutcome<-"PerfSpellResol_Type_Hist"
  
  
  # --- Estimate survival rate of the censoring event G(t) = P(C >= t) for time-to-censoring variable C
  # This prepares for implementing the Inverse Probability of Censoring Weighting (IPCW) scheme,
  # as an adjustment to event rates (discrete density), as well as for the time-dependent Brier score
  
  # - Compute Kaplan-Meier survival estimates (product-limit) for censoring-event | Spell-level with right-censoring & left-truncation
  km_Censoring <- survfit(
    as.formula(paste0("Surv(time=",fldStart,", time2=", fldStop, ", event=", fldCensored, "==1, type='counting') ~ 1")), 
                          id=get(fldKey), data=datGiven)
  # summary(km_Censoring)$table # overall summary statistics
  # plot(km_Censoring, conf.int = T) # survival curve
  datG <- data.table(summary(km_Censoring)$time, summary(km_Censoring)$surv)
  setnames(datG, c(fldStop, "G_t"))
  datG_ti <- data.table(EventTime = summary(km_Censoring)$time, G_Ti = summary(km_Censoring)$surv)
  
  
  # --- Calculate survival quantities of interest
  # Predict hazard h(t) = P(T=t | T>= t) in discrete-time
  datGiven[, Hazard := predict(modGiven, newdata=datGiven, type = predType)]
  if (predType=="response") {
    # Derive survival probability S(t) = \prod ( 1- hazard), based on output of predict()
    datGiven[, Survival := cumprod(1-Hazard), by=list(get(fldKey))] 
  } else if (predType=="exp") {
    # Calculate survival probability S(t,x)=exp(-H(t,x)), based on output of predict()
    # NOTE: Hazard is actually the cumulative hazard in this context
    datGiven[, Survival := exp(-Hazard), by=list(get(fldKey))] 
  }
  # - Merge censoring survival probability unto main set
  datGiven <- merge(datGiven, datG, by=fldStop, all.x=T)
  datGiven[is.na(G_t), G_t := 1] # Fill any missing G_t with 1 (no censoring information)
  
  
  # --- Compute Censoring survival probability at age T_i for each subject (i,j), i..e, G_hat(T_i)
  # This requires mapping each subject's observed event time (Spell age) to g_hat(T_i)
  datGiven <- merge(datGiven, datG_ti, by.x=fldSpellAge, by.y="EventTime", all.x=T)
  datGiven[is.na(G_Ti), G_Ti := 1]
  
  
  # --- Compute constituent quantities of time-dependent Brier score (tBS)
  # - Compute binary event indicator at time t: y(t) = I(T>t)
  datGiven[, y_t := as.integer(get(fldSpellAge) > get(fldStop))]
  # - Compute weighted Brier score at each time point
  datGiven[, Weight_w := fifelse(get(fldSpellAge) > get(fldStop), 1/G_t,
                                 fifelse(get(fldSpellAge) <= get(fldStop) & get(fldSpellOutcome)=="Defaulted", 1/G_Ti,
                                         0)), by=list(get(fldKey))]
  # -- Compute squared error loss per row ---
  datGiven[, SquaredError := (y_t - Survival)^2]
  
  
  # --- Aggregate and compute weighted Brier score per time point
  datTBS <- datGiven[Weight_w > 0 & get(fldStop)<=spellPeriodMax, .(
    mean(Weight_w * SquaredError, na.rm = TRUE)
  ), by=list(get(fldStop))]
  setnames(datTBS, c(fldStop, "Brier"))
  # Subset only to specified max
  IBS <- mean(datTBS$Brier, na.rm=T)
  
  return(list(tBS=datTBS, IBS=IBS))

  # - Cleanup (useful during interactive debugging of function)
  rm(datGiven, modGiven, fldKey, fldStart, fldStop, fldCensored, fldSpellAge, predType,
     km_Censoring, datTBS,datG, datG_ti)
}



# ----------------- 4. Graph tBS-values across models

# --- Discrete-time hazard models

# - Data fusion across models
datGraph <- rbind(data.table(objCoxDisc_bas$tBS, Type="a_Basic"),
                  data.table(objCoxDisc_adv$tBS, Type="b_Advanced"))

# - Aesthetic engineering
datGraph[, FacetLabel := "Discrete-time hazard models"]
zoomedSpellAge <- 120

# - Recalcualte Integrated tBS over zoomed spell age
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
  annotate(geom="text", x=50, y=1.5, label=paste0("IBS (Basic): ", round(objCoxDisc_bas$IBS,3)), 
           family=chosenFont, size=3.5, colour=vCol[1]) + 
  annotate(geom="text", x=43, y=1, label=paste0("IBS (Advanced): ", round(objCoxDisc_adv$IBS,3)), 
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
  annotate(geom="text", x=50, y=0.07, label=paste0("IBS (Basic): ", round(ibs_bas,3)), 
           family=chosenFont, size=3.5, colour=vCol[1]) + 
  annotate(geom="text", x=45, y=0.063, label=paste0("IBS (Advanced): ", round(ibs_adv,3)), 
           family=chosenFont, size=3.5, colour=vCol[2]) +   
  # Facets & scales
  scale_colour_manual(name="", values=vCol, labels=vLabel) + 
  scale_linetype_discrete(name="", labels=vLabel) +
  scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
  scale_x_continuous(breaks=breaks_pretty(), label=comma)
)

# - Combining the two above plots onto a single graph
(plot.full <- gOuter + annotation_custom(grob = ggplotGrob(gInner), xmin=0, xmax=200, ymin=3, ymax=13))


# - Save plot
dpi <- 170
ggsave(plot.full, file=paste0(genFigPath,"tBrierScores_CoxDisc.png"),width=1350/dpi, height=1000/dpi,dpi=dpi, bg="white")





plot(objCoxDisc_adv$tBS$TimeInPerfSpell, objCoxDisc_adv$tBS$Brier, type = "s", lwd = 2, col = "blue",
     ylab = "Brier Score", xlab = "Time", main = "Time-dependent Brier Score")
lines(objCoxDisc_bas$tBS$TimeInPerfSpell, objCoxDisc_bas$tBS$Brier, type="s", lwd=2, col="red")
abline(h = objCoxDisc_adv$IBS, col = "blue", lty = 2)
abline(h = objCoxDisc_bas$IBS, col = "red", lty = 2)


##### TO-DO
# - Compare tBS for Cox PH model with that obtained from pec(). And debug!

plot(objCoxPH_bas$tBS$TimeInPerfSpell, objCoxPH_bas$tBS$Brier, type = "s", lwd = 2, col = "red",
     ylab = "Brier Score", xlab = "Time", main = "Time-dependent Brier Score")
abline(h = objCoxPH_bas$IBS, col = "red", lty = 2)
