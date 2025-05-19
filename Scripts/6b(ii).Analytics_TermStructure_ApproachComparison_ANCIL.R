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
setDT(datCredit_train, key="PerfSpell_Key")
setDT(datCredit_valid, key="PerfSpell_Key")





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
# using a single Cox regression model (TFD-definition)


# ------ Approach 1: Looped individual survfit()-calls by ID

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid$PerfSpell_Key)


# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
#cl.port <- makeCluster(round(6)); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_TFD <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                       .packages=c('data.table', 'survival'), .export=c('survQuants')) %do%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_TFD,
                          it=j, numKeys=numSpellKeys, genPath=genPath)
  } # ----------------- End of Inner Loop -----------------
#stopCluster(cl.port); 
proc.time() - ptm; # 54h

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_TFD"), datSurv_TFD)



# --- Graphing the event density / probability mass function f(t)

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_TFD')) unpack.ffdf(paste0(genPath,"datSurv_TFD"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), tempPath);gc()
setDT(datSurv_TFD, key="End")

# - Determine population average survival and event rate across loans per time period
datAggr <- datSurv_TFD[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
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
datGraph[, FacetLabel := "Approach 1"]

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
    annotate("text", y=0.006,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
              size = 3) + 
    # Scales and options
    facet_grid(FacetLabel ~ .) + 
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma)
)

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "TFD/Approach1_EventProb_", mainEventName,"_SpellLevel_FirstSpell_ActVsExp.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")





# ------ Approach 2: Generic survfit()-call by ID

ptm <- proc.time() #IGNORE: for computation time calculation
survFit_pred_TFD <- survfit(cox_TFD, centered=F, newdata=datCredit_valid, id=PerfSpell_Key)
proc.time() - ptm # 86m

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"survFit_TFD"), survFit_pred_TFD)

# - Confirm prepared datasets are loaded into memory
if (!exists('survFit_pred_TFD')) unpack.ffdf(paste0(genPath,"survFit_TFD"), tempPath);gc()

# - Create survival data object for graphing and comparison purposes
(datSurv_test <- data.table(End=survFit_pred_TFD$time,
                            CHaz=survFit_pred_TFD$cumhaz, RiskSetSize=survFit_pred_TFD$n.risk,
                            NumEvents=survFit_pred_TFD$n.event, NumCensored=survFit_pred_TFD$n.censor,
                            Survival=round(survFit_pred_TFD$surv,digits=15)))
# NOTE: The times-object within survFit_pred_TFD is only an index-object, not true stop-times.


# --- Melt vectors into ID-specific matrix-form: rows (unique failure times) by columns (SpellKey)
vOnes <- which(survFit_pred_TFD$time == 1) # designates the starting index of new IDs within the broader vector
matSurv <- matrix(NA, nrow=max(survFit_pred_TFD$time), ncol=length(vOnes))
matCHaz <- matrix(NA, nrow=max(survFit_pred_TFD$time), ncol=length(vOnes))
matHaz <- matrix(NA, nrow=max(survFit_pred_TFD$time), ncol=length(vOnes))
matEventProb <- matrix(NA, nrow=max(survFit_pred_TFD$time), ncol=length(vOnes))
k <-2
t <- 1

for (i in 1:length(survFit_pred_TFD$time)) {
  if (k <= length(vOnes)) {
    
    if (i < vOnes[k]) {
      matSurv[t,k-1] <- survFit_pred_TFD$surv[i]
      matCHaz[t,k-1] <- survFit_pred_TFD$cumhaz[i]
      if (t != 1) matHaz[t,k-1] <- (matCHaz[t,k-1] - matCHaz[t-1,k-1]) else matHaz[t,k-1] <- matCHaz[t,k-1]
      if (t != 1) matEventProb[t,k-1] <- matHaz[t,k-1] * matSurv[t-1,k-1] else matEventProb[t,k-1] <- matHaz[t,k-1] * 1
      t <- t+1 # increment period
    } else {
      # resetting indices
      k <- k+1; t <- 1
      matSurv[t,k-1] <- survFit_pred_TFD$surv[i]
      matCHaz[t,k-1] <- survFit_pred_TFD$cumhaz[i]
      matHaz[t,k-1] <- matCHaz[t,k-1]
      matEventProb[t,k-1] <- matHaz[t,k-1] * 1
      t <- t+1 # increment period
    }
  }
}

# [TEST] Case 1
matSurv[!is.na(matSurv[,1]), 1]
survFit_pred_TFD$surv[1: (vOnes[2]-1)]
# [TEST] Case 2
matSurv[!is.na(matSurv[,2]), 2]
survFit_pred_TFD$surv[vOnes[2]: (vOnes[3]-1)]



# --- Graphing the event density / probability mass function f(t)

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_TFD')) unpack.ffdf(paste0(genPath,"datSurv_TFD"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), tempPath);gc()
setDT(datSurv_TFD, key="End")

# - Determine population average survival and event rate across loans per time period
datAggr <- data.table(End=unique(survFit_pred_TFD$time), EventRate=rowMeans(matEventProb, na.rm=T))
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
ggsave(gsurv_ft, file=paste0(genFigPath, "TFD/Approach2_EventProb_", mainEventName,"_SpellLevel_FirstSpell_ActVsExp.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")






# ------ Approach 3: Looped individual survfit()-calls by ID | Custom calculation

test <- subset(datCredit_valid, LoanID %in% unique(datCredit_valid[PerfSpell_Age > 5 & PerfSpell_Num > 2,LoanID])[1],
               select=c("LoanID", "Date", vecVars_TFD, "PerfSpell_Key", "PerfSpell_Num","PerfSpell_Counter","Start", "End", "Default_Ind"))


# --- Compute baseline hazard
datHaz <- survQuants.data(datGiven_train=datCredit_train, datGiven_score=test, 
                          vars=vecVars_TFD, beta=cox_TFD$coefficients, fldID="PerfSpell_Key", 
                          fldStart="Start", fldStop="End", fldEvent="Default_Ind",
                          centered = T, ties="Breslow", coxGiven=cox_TFD)
### RESULTS: Method is flawed since it provides negative event probabilities

# - Compare
datFusion <- merge(datSurv_TFD, datHaz[,list(PerfSpell_Key=Key, End=Time, Survival_custom=Survival,
                                             CHaz_custom=CHaz, EventProb_Custom=EventProb)], by=c("PerfSpell_Key", "End"))
plot(datFusion$Survival, type="b"); plot(datFusion$Survival_custom, col="red", type="b")
plot(datFusion$Survival-datFusion$Survival_custom, type="b")





# ------ Approach 4: Looped individual survfit()-calls by ID (multi-threaded)

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid$PerfSpell_Key)

rm(datCredit_train)

# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
cl.port <- makeCluster(round(6)); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_TFD2 <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                       .packages=c('data.table', 'survival'), .export=c('survQuants')) %dopar%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_TFD,
                          it=j, numKeys=numSpellKeys, genPath=genPath)
  } # ----------------- End of Inner Loop -----------------
stopCluster(cl.port); 
proc.time() - ptm # 54h

