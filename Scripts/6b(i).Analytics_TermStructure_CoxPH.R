# ====================================== TERM-STRUCTURES: ACTUAL VS EXPECTED ===========================
# Calculate and compare the term-structures of default risk
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
                         id=PerfSpell_Key, data=datCredit_train, ties="efron", model=T)



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
                       id=PerfSpell_Key, datCredit_train, ties="efron", model=T)





# ----------------- 3. Term-structure of default risk
# Compare various approaches towards deriving the term-structure of default risk,
# using a single Cox regression model (PWP-definition)

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid$PerfSpell_Key)



# ------ Approach 1: Looped individual survfit()-calls by ID (multi-threaded) | Basic CPH-model

# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
cl.port <- makeCluster(round(8)); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_PWPST_basic <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                         .packages=c('data.table', 'survival'), .export=c('survQuants')) %dopar%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_PWPST_basic,
                          it=j, numKeys=numSpellKeys, genPath=genPath,
                          timeVar="TimeInPerfSpell", startVar="Start")
  } # ----------------- End of Inner Loop -----------------
stopCluster(cl.port); 
proc.time() - ptm; # 

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_PWPST_basic"), datSurv_PWPST_basic)




# ------ Approach 1: Looped individual survfit()-calls by ID (multi-threaded) | Advanced CPH-model

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
proc.time() - ptm;

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


# -- Aggregate data by Time in Performing Spell

# - General parameters
sMaxSpellAge <- 300 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))
sMaxSpellAge_graph <- 300 # max for [PerfSpell_Age] for graphing purposes

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_PWPST')) unpack.ffdf(paste0(genPath,"datSurv_PWPST"), tempPath);gc()
if (!exists('datSurv_PWPST_basic')) unpack.ffdf(paste0(genPath,"datSurv_PWPST_basic"), tempPath);gc()
setDT(datSurv_PWPST, key="End"); setDT(datSurv_PWPST_basic, key="End")
#if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), tempPath);gc()

# - Determine population average survival and event rate across loans per time period
datAggr <- datSurv_PWPST[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
datAggr_basic <- datSurv_PWPST_basic[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datSurv_sub[Time <= 300, Time], datSurv_sub[Time <= 300, EventRate], type="b") # actual
plot(datAggr_basic[End <= sMaxSpellAge, End], datAggr_basic[End <= sMaxSpellAge, EventRate], type="b") # expected
plot(datAggr[End <= sMaxSpellAge, End], datAggr[End <= sMaxSpellAge, EventRate], type="b") # expected

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 12
smthEventRate_Act <- lm(EventRate ~ ns(Time, df=sDf_Act), data=datSurv_sub[Time <= sMaxSpellAge,])
smthEventRate_Exp_bas <- lm(EventRate ~ ns(End, df=sDf_Exp), data=datAggr_basic[End <= sMaxSpellAge])
smthEventRate_Exp_adv <- lm(EventRate ~ ns(End, df=sDf_Exp), data=datAggr[End <= sMaxSpellAge])
summary(smthEventRate_Act); summary(smthEventRate_Exp_bas); summary(smthEventRate_Exp_adv)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth_Act <- predict(smthEventRate_Act, newdata=datSurv_sub, se=T)
vPredSmth_Exp_bas <- predict(smthEventRate_Exp_bas, newdata=datAggr_basic, se=T)
vPredSmth_Exp_adv <- predict(smthEventRate_Exp_adv, newdata=datAggr, se=T)

# - Add smoothed estimate to graphing object
datSurv_sub[, EventRate_spline := vPredSmth_Act$fit]
datAggr_basic[, EventRate_spline := vPredSmth_Exp_bas$fit]
datAggr[, EventRate_spline := vPredSmth_Exp_adv$fit]

# - Create graphing data object
datGraph <- rbind(datSurv_sub[,list(Time, EventRate, Type="a_Actual")],
                  datSurv_sub[,list(Time, EventRate=EventRate_spline, Type="b_Actual_spline")],
                  datAggr_basic[, list(Time=End, EventRate, Type="c_Expected_bas")],
                  datAggr_basic[, list(Time=End, EventRate=EventRate_spline, Type="d_Expected_spline_bas")],
                  datAggr[, list(Time=End, EventRate, Type="e_Expected_adv")],
                  datAggr[, list(Time=End, EventRate=EventRate_spline, Type="f_Expected_spline_adv")]
)

# - Create different groupings for graphing purposes
datGraph[Type %in% c("a_Actual","c_Expected_bas", "e_Expected_adv"), EventRatePoint := EventRate ]
datGraph[Type %in% c("b_Actual_spline","d_Expected_spline_bas", "f_Expected_spline_adv"), EventRateLine := EventRate ]
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) spell-time"]

# - Aesthetic engineering
chosenFont <- "Cambria"
vCol_Point <- brewer.pal(8, "Pastel1")[c(3,1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(3,1,2)]
mainEventName <- "Default"

# - Calculate MAE between event rates
datFusion_bas <- merge(datSurv_sub[Time <= sMaxSpellAge], 
                   datAggr_basic[End <= sMaxSpellAge,list(Time=End, EventRate_Exp=EventRate)], by="Time")
datFusion_adv <- merge(datSurv_sub[Time <= sMaxSpellAge], 
                       datAggr[End <= sMaxSpellAge,list(Time=End, EventRate_Exp=EventRate)], by="Time")
MAE_eventProb_bas <- mean(abs(datFusion_bas$EventRate - datFusion_bas$EventRate_Exp), na.rm=T)
MAE_eventProb_adv <- mean(abs(datFusion_adv$EventRate - datFusion_adv$EventRate_Exp), na.rm=T)

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(3,4,5,6, 1,2)]
vLabel2 <- c("b_Actual_spline"=paste0("Actual spline"), 
             "d_Expected_spline_bas"=paste0("Exp spline: Basic"),
             "f_Expected_spline_adv"=paste0("Exp spline: Advanced"),
             "a_Actual"="Actual", "c_Expected_bas"="Exp: Basic", "e_Expected_adv"="Exp: Advanced" )
vSize <- c(0.2,0.3,0.2,0.3,0.2,0.3)
vLineType <- c("dashed", "solid", "dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote("Performing spell age (months)"*~italic(t)),
         #subtitle="Term-structures of default risk: Cox PH"
         ) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=0.6) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.0065,x=100, label=paste0("MAE (Basic): ", percent(MAE_eventProb_bas, accuracy=0.0001)), family=chosenFont,
             size = 3) + 
    annotate("text", y=0.0055,x=100, label=paste0("MAE (Advanced): ", percent(MAE_eventProb_adv, accuracy=0.0001)), family=chosenFont,
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
dpi <- 280 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "EventProb_", mainEventName,"_ActVsExp_CoxPH.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")


# - Cleanup
rm(datAggr, datAggr_basic, datFusion_adv, datFusion_bas, datGraph)
