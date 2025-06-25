# ========================= ANALYTICS: TERM-STRUCTURE ==========================
# Ancillary script to compare different approaches to deriving the term-structure 
# of default risk across time from the various Cox regression models
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# --------------------------------------------------------------------------------
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

# -- Inputs:
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b

#
# -- Outputs:
#   - <Analytics> | tROC-graphs
# ================================================================================





# ----------------- 1. Load data

# - General parameters
sMaxSpellAge <- 300 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))

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
setDT(datCredit_train, key=c("PerfSpell_Key", "TimeInPerfSpell"))
setDT(datCredit_valid, key=c("PerfSpell_Key", "TimeInPerfSpell"))





# ----------------- 2. Fit a Cox regression model on the resampled prepared data

# ------ Prentice-Williams-Peterson (PWP) Spell-time definition | Advanced
# - Initialize variables
vecVars_PWPST_adv <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST_adv <- coxph(as.formula(paste0("Surv(TimeInPerfSpell-1,TimeInPerfSpell,DefaultStatus1) ~ ", 
                                         paste(vecVars_PWPST_adv,collapse=" + "), 
                                         " + strata(PerfSpell_Num_binned)")),
                       id=PerfSpell_Key, datCredit_train, ties="efron", model=T)
summary(cox_PWPST_adv); AIC(cox_PWPST_adv); concordance(cox_PWPST_adv)




# ----------------- 3. Term-structure of default risk
# Compare various approaches towards deriving the term-structure of default risk,
# using a single Cox regression model (PWP-definition)


# ------ Approach 1: Looped individual survfit()-calls by ID (multi-threaded)

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid$PerfSpell_Key)


# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
cl.port <- makeCluster(round(8)); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_PWPST <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                       .packages=c('data.table', 'survival'), .export=c('survQuants')) %dopar%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_PWPST_adv,
                          it=j, numKeys=numSpellKeys, genPath=genPath,
                          timeVar="TimeInPerfSpell", startVar="Start")
  } # ----------------- End of Inner Loop -----------------
stopCluster(cl.port); 
proc.time() - ptm; # 54h

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_PWPST"), datSurv_PWPST)



# ------- Approach 1: Graphing the event density / probability mass function f(t)

# -- 1. Kaplan-Meier estimation towards constructing empirical term-structure of default risk

# --- Preliminaries
# - Create pointer to the appropriate data object 
datCredit <- rbind(datCredit_train, datCredit_valid)


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

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_PWPST')) unpack.ffdf(paste0(genPath,"datSurv_PWPST"), tempPath);gc()
setDT(datSurv_PWPST, key="End")
#if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), tempPath);gc()


# - Determine population average survival and event rate across loans per time period
datAggr <- datSurv_PWPST[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datAggr[End <= sMaxSpellAge, End], datAggr[End <= sMaxSpellAge, EventRate], type="b")
plot(datSurv_sub[Time <= 300, Time], datSurv_sub[Time <= 300, EventRate], type="b")

# - General parameters
sMaxSpellAge <- 300 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))
sMaxSpellAge_graph <- 300 # max for [PerfSpell_Age] for graphing purposes

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 12
smthEventRate_Act <- lm(EventRate ~ ns(Time, df=sDf_Act), data=datSurv_sub[Time <= sMaxSpellAge,])
smthEventRate_Exp <- lm(EventRate ~ ns(End, df=sDf_Exp), data=datAggr[End <= sMaxSpellAge])
summary(smthEventRate_Act); summary(smthEventRate_Exp)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth_Act <- predict(smthEventRate_Act, newdata=datSurv_sub, se=T)
vPredSmth_Exp <- predict(smthEventRate_Exp, newdata=datAggr, se=T)

# - Add smoothed estimate to graphing object
datSurv_sub[, EventRate_spline := vPredSmth_Act$fit]
datAggr[, EventRate_spline := vPredSmth_Exp$fit]

# - Create graphing data object
datGraph <- rbind(datSurv_sub[,list(Time, EventRate, Type="a_Actual")],
                  datSurv_sub[,list(Time, EventRate=EventRate_spline, Type="b_Actual_spline")],
                  datAggr[, list(Time=End, EventRate, Type="c_Expected")],
                  datAggr[, list(Time=End, EventRate=EventRate_spline, Type="d_Expected_spline")]
)

# - Create different groupings for graphing purposes
datGraph[Type %in% c("a_Actual","c_Expected"), EventRatePoint := EventRate ]
datGraph[Type %in% c("b_Actual_spline","d_Expected_spline"), EventRateLine := EventRate ]
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) spell-time"]

# - Aesthetic engineering
chosenFont <- "Cambria"
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"

# - Calculate MAE between event rates
datFusion <- merge(datSurv_sub[Time <= sMaxSpellAge], 
                   datAggr[End <= sMaxSpellAge,list(Time=End, EventRate_Exp=EventRate)], by="Time")
MAE_eventProb <- mean(abs(datFusion$EventRate - datFusion$EventRate_Exp), na.rm=T)

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(3,4,1,2)]
vLabel2 <- c("b_Actual_spline"=paste0("Actual spline (df=",sDf_Act,")"), 
             "d_Expected_spline"=paste0("Expected spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Expected")
vSize <- c(0.5,0.6,0.5,0.6)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.0025,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
             size = 3) + 
    # Scales and options
    facet_grid(FacetLabel ~ .) + 
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(color=guide_legend(nrow=2))
)

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "Approach1_EventProb_", mainEventName,"_ActVsExp_CoxPH.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")





# ------ Approach 2: Generic survfit()-call by ID

ptm <- proc.time() #IGNORE: for computation time calculation
survFit_pred_PWPST <- survfit(cox_PWPST_adv, centered=F, newdata=datCredit_valid, id=PerfSpell_Key)
proc.time() - ptm # 86m

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"survFit_PWPST"), survFit_pred_PWPST)

# - Confirm prepared datasets are loaded into memory
if (!exists('survFit_pred_PWPST')) unpack.ffdf(paste0(genPath,"survFit_PWPST"), tempPath);gc()

# - Create survival data object for graphing and comparison purposes
(datSurv_test <- data.table(End=survFit_pred_PWPST$time,
                            CHaz=survFit_pred_PWPST$cumhaz, RiskSetSize=survFit_pred_PWPST$n.risk,
                            NumEvents=survFit_pred_PWPST$n.event, NumCensored=survFit_pred_PWPST$n.censor,
                            Survival=round(survFit_pred_PWPST$surv,digits=15)))
# NOTE: The times-object within survFit_pred_PWPST is only an index-object, not true stop-times.


# --- Melt vectors into ID-specific matrix-form: rows (unique failure times) by columns (SpellKey)
vOnes <- which(survFit_pred_PWPST$time == 1) # designates the starting index of new IDs within the broader vector
matSurv <- matrix(NA, nrow=max(survFit_pred_PWPST$time), ncol=length(vOnes))
matCHaz <- matrix(NA, nrow=max(survFit_pred_PWPST$time), ncol=length(vOnes))
matHaz <- matrix(NA, nrow=max(survFit_pred_PWPST$time), ncol=length(vOnes))
matEventProb <- matrix(NA, nrow=max(survFit_pred_PWPST$time), ncol=length(vOnes))
k <-2
t <- 1

for (i in 1:length(survFit_pred_PWPST$time)) {
  if (k <= length(vOnes)) {
    
    if (i < vOnes[k]) {
      matSurv[t,k-1] <- survFit_pred_PWPST$surv[i]
      matCHaz[t,k-1] <- survFit_pred_PWPST$cumhaz[i]
      if (t != 1) matHaz[t,k-1] <- (matCHaz[t,k-1] - matCHaz[t-1,k-1]) else matHaz[t,k-1] <- matCHaz[t,k-1]
      if (t != 1) matEventProb[t,k-1] <- matHaz[t,k-1] * matSurv[t-1,k-1] else matEventProb[t,k-1] <- matHaz[t,k-1] * 1
      t <- t+1 # increment period
    } else {
      # resetting indices
      k <- k+1; t <- 1
      matSurv[t,k-1] <- survFit_pred_PWPST$surv[i]
      matCHaz[t,k-1] <- survFit_pred_PWPST$cumhaz[i]
      matHaz[t,k-1] <- matCHaz[t,k-1]
      matEventProb[t,k-1] <- matHaz[t,k-1] * 1
      t <- t+1 # increment period
    }
  }
}

# [TEST] Case 1
matSurv[!is.na(matSurv[,1]), 1]
survFit_pred_PWPST$surv[1: (vOnes[2]-1)]
# [TEST] Case 2
matSurv[!is.na(matSurv[,2]), 2]
survFit_pred_PWPST$surv[vOnes[2]: (vOnes[3]-1)]



# ------ Approach 2:Graphing the event density / probability mass function f(t)

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), tempPath);gc()

# - Determine population average survival and event rate across loans per time period
datAggr <- data.table(End=unique(survFit_pred_PWPST$time), EventRate=rowMeans(matEventProb, na.rm=T))
plot(datAggr[End <= sMaxSpellAge, End], datAggr[End <= sMaxSpellAge, EventRate], type="b")

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 18
smthEventRate_Act <- lm(EventRate ~ ns(Time, df=sDf_Act), data=datSurv[Time <= sMaxSpellAge,])
smthEventRate_Exp <- lm(EventRate ~ ns(End, df=sDf_Exp), data=datAggr[End <= sMaxSpellAge])
summary(smthEventRate_Act); summary(smthEventRate_Exp)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth_Act <- predict(smthEventRate_Act, newdata=datSurv, se=T)
vPredSmth_Exp <- predict(smthEventRate_Exp, newdata=datAggr, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth_Act$fit]
datAggr[, EventRate_spline := vPredSmth_Exp$fit]

# - Create graphing data object
datGraph <- rbind(datSurv[,list(Time, EventRate, Type="a_Actual")],
                  datSurv[,list(Time, EventRate=EventRate_spline, Type="b_Actual_spline")],
                  datAggr[, list(Time=End, EventRate, Type="c_Expected")],
                  datAggr[, list(Time=End, EventRate=EventRate_spline, Type="d_Expected_spline")]
)

# - Create different groupings for graphing purposes
datGraph[Type %in% c("a_Actual","c_Expected"), EventRatePoint := EventRate ]
datGraph[Type %in% c("b_Actual_spline","d_Expected_spline"), EventRateLine := EventRate ]
datGraph[, FacetLabel := "Approach 2"]

# - Aesthetic engineering
chosenFont <- "Cambria"
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"

# - Calculate MAE between event rates
datFusion <- merge(datSurv[Time <= sMaxSpellAge], 
                   datAggr[End <= sMaxSpellAge,list(Time=End, EventRate_Exp=EventRate)], by="Time")
MAE_eventProb <- mean(abs(datFusion$EventRate - datFusion$EventRate_Exp), na.rm=T)

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(3,4,1,2)]
vLabel2 <- c("b_Actual_spline"=paste0("Actual natural spline (df=",sDf_Act,")"), 
             "d_Expected_spline"=paste0("Scored natural spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Scored")
vSize <- c(0.5,0.75,0.5,0.75)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.005,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
             size = 3) + 
    # Scales and options
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma)
)

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "Approach2_EventProb_", mainEventName,"_ActVsExp_CoxPH.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")





# ------ Approach 3: Using predict() with type = "expected"

# -- Assumption: The eventual lagging of [Survival] assumes TimeinPerfSpell=1 at the start, since only then will S(t0)=1.
# But left-truncated accounts violate this assumption.
(sum(datCredit_valid[PerfSpell_Counter == 1, .(Faulty=ifelse(TimeInPerfSpell > 1, 1, 0)), by=list(PerfSpell_Key)]$Faulty) / 
    datCredit_valid[PerfSpell_Counter == 1, .N])
### RESULTS: 38% of accounts violate this assumption

# - Extract those left-truncated cases, add a record, decrement the time variable, and merge back
datAdd <- datCredit_valid[PerfSpell_Counter == 1 & TimeInPerfSpell>1, ]
datAdd[, TimeInPerfSpell := TimeInPerfSpell - 1]
datAdd[, Added := T]; datCredit_valid[, Added := F]
datCredit_valid <- rbind(datAdd, datCredit_valid); rm(datAdd); gc()
setDT(datCredit_valid, key=c("PerfSpell_Key", "TimeInPerfSpell"))
#test <- subset(datCredit_valid, PerfSpell_Key == vSpellKeys[j])

# - Score cases by estimating the cumulative hazard function H(t,x)
datCredit_valid[, CHaz:=predict(cox_PWPST_adv, new=datCredit_valid, type="expected", id=PerfSpell_Key)]
# Calculate survival probability S(t,x)=exp(-H(t,x))
datCredit_valid[, Survival := exp(-CHaz)]
datCredit_valid[, Survival_1 := shift(Survival, n=1, type="lag", fill=1), by=list(PerfSpell_Key)]
# Calculate key survival quantities from S(t,x)
datCredit_valid[, Hazard := (Survival_1 - Survival)/Survival_1]
#datCredit_valid[, EventRate := Survival_1 * Hazard]
datCredit_valid[, EventRate := Survival_1 - Survival]

# Remove added cases
datCredit_valid <- subset(datCredit_valid, Added == F)
datCredit_valid[, Added := NULL]



# ------ Approach 3: Graphing the event density / probability mass function f(t)

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), tempPath);gc()
setDT(datSurv_PWPST, key="End")

# - Determine population average survival and event rate across loans per time period
datAggr <- datCredit_valid[, list(EventRate = mean(EventRate,na.rm=T), Freq=.N),by=list(TimeInPerfSpell)]
plot(datAggr[TimeInPerfSpell <= 300, TimeInPerfSpell], datAggr[TimeInPerfSpell <= 300, EventRate], type="b")
plot(datSurv[Time <= 300, Time], datSurv[Time <= 300, EventRate], type="b")


# - General parameters
sMaxSpellAge <- 240 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))
sMaxSpellAge_graph <- 240 # max for [PerfSpell_Age] for graphing purposes

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 12
smthEventRate_Act <- lm(EventRate ~ ns(Time, df=sDf_Act), data=datSurv[Time <= sMaxSpellAge,])
smthEventRate_Exp <- lm(EventRate ~ ns(TimeInPerfSpell, df=sDf_Exp), data=datAggr[TimeInPerfSpell <= sMaxSpellAge])
summary(smthEventRate_Act); summary(smthEventRate_Exp)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth_Act <- predict(smthEventRate_Act, newdata=datSurv, se=T)
vPredSmth_Exp <- predict(smthEventRate_Exp, newdata=datAggr, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth_Act$fit]
datAggr[, EventRate_spline := vPredSmth_Exp$fit]

# - Create graphing data object
datGraph <- rbind(datSurv[,list(Time, EventRate, Type="a_Actual")],
                  datSurv[,list(Time, EventRate=EventRate_spline, Type="b_Actual_spline")],
                  datAggr[, list(Time=TimeInPerfSpell, EventRate, Type="c_Expected")],
                  datAggr[, list(Time=TimeInPerfSpell, EventRate=EventRate_spline, Type="d_Expected_spline")]
)

# - Create different groupings for graphing purposes
datGraph[Type %in% c("a_Actual","c_Expected"), EventRatePoint := EventRate ]
datGraph[Type %in% c("b_Actual_spline","d_Expected_spline"), EventRateLine := EventRate ]
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP): Advanced"]

# - Aesthetic engineering
chosenFont <- "Cambria"
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"

# - Calculate MAE between event rates
datFusion <- merge(datSurv[Time <= sMaxSpellAge], 
                   datAggr[TimeInPerfSpell <= sMaxSpellAge,list(Time=TimeInPerfSpell, EventRate_Exp=EventRate)], by="Time")
MAE_eventProb <- mean(abs(datFusion$EventRate - datFusion$EventRate_Exp), na.rm=T)

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(3,4,1,2)]
vLabel2 <- c("b_Actual_spline"=paste0("Actual spline (df=",sDf_Act,")"), 
             "d_Expected_spline"=paste0("Expected spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Expected")
vSize <- c(0.5,0.6,0.5,0.6)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.0025,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
             size = 3) + 
    # Scales and options
    facet_grid(FacetLabel ~ .) + 
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(color=guide_legend(nrow=2))
)

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "Approach3_EventProb_", mainEventName,"_ActVsExp_CoxPH.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")


# -- Test case comparison
# Extract test case
test1 <- subset(datSurv_PWPST, PerfSpell_Key %in% vSpellKeys[1])
# Extract test case
test2 <- subset(datCredit_valid, PerfSpell_Key %in% vSpellKeys[1],
               select=c("LoanID", "Date", vecVars_PWPST_adv, "PerfSpell_Key", "PerfSpell_Num","PerfSpell_Counter","Start", 
                        "TimeInPerfSpell", "DefaultStatus1", "CHaz", "Survival", "EventRate"))
plot(test1$Survival); plot(test2$Survival)
### RESULTS: Despite the simplicity of Approach 3 (and runtime reduction), it produces unintuitive results in 
# the survival probability for a single test case, as opposed to those from survQuants() from Approach 1.
# The resulting term-structure looks similar to the KM-analysis (event probabilities), though has greater discrepancy.
# Lastly, the event probabilities are negative in some cases, which violates the axioms of probabilities.







# ------ Approach 4: Looped individual survfit()-calls by ID | Custom calculation

# --- Compute baseline hazard
datHaz <- survQuants.data(datGiven_train=datCredit_train, datGiven_score=test, 
                          vars=vecVars_PWPST_adv, beta=cox_PWPST_adv$coefficients, fldID="PerfSpell_Key", 
                          fldStart="Start", fldStop="End", fldEvent="Default_Ind",
                          centered = T, ties="Breslow", coxGiven=cox_PWPST_adv)
### RESULTS: Method is flawed since it provides negative event probabilities

# - Compare
datFusion <- merge(datSurv_PWPST, datHaz[,list(PerfSpell_Key=Key, End=Time, Survival_custom=Survival,
                                               CHaz_custom=CHaz, EventProb_Custom=EventProb)], by=c("PerfSpell_Key", "End"))
plot(datFusion$Survival, type="b"); plot(datFusion$Survival_custom, col="red", type="b")
plot(datFusion$Survival-datFusion$Survival_custom, type="b")