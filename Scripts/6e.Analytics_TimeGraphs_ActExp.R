# ======================= TIME GRAPHS OF EVENT RATES: ACTUAL VS EXPECTED ===============================
# Calculate and compare time graphs of the actual vs expected resolution & 12-month default rates
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
#   - 6b(i).Analytics_TermStructure_CoxPH.R

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

# - Weigh default cases heavier. as determined interactively based on calibration success (script 6e)
datCredit_train[, Weight := ifelse(DefaultStatus1==1,10,1)]
datCredit_valid[, Weight := ifelse(DefaultStatus1==1,10,1)] # for merging purposes





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




# ----------------- 3. Actual vs expected 12-month default rates | Discrete-time hazard models

# --- Preliminaries
# - Create pointer to the appropriate data object 
#datCredit <- rbind(datCredit_train_PWPST, datCredit_valid_PWPST)
datCredit <- datCredit_valid_PWPST


# --- Handle left-truncated spells by adding a starting record 
# This is necessary for calculating certain survival quantities later
datAdd <- subset(datCredit, Counter == 1 & TimeInPerfSpell > 1)
datAdd[, TimeInPerfSpell := TimeInPerfSpell-1]
datAdd[, Counter := 0]
datCredit <- rbind(datCredit, datAdd)
setorder(datCredit, PerfSpell_Key, TimeInPerfSpell) # Preparing for survival model scoring


# --- Calculate survival quantities of interest
# Predict hazard h(t) = P(T=t | T>= t) in discrete-time
datCredit[!is.na(PerfSpell_Num), Hazard_adv := predict(modLR, newdata=.SD[], type = "response")]
datCredit[!is.na(PerfSpell_Num), Hazard_bas := predict(modLR_basic, newdata=.SD[], type = "response")]
# Derive survival probability S(t) = \prod ( 1- hazard)
datCredit[!is.na(PerfSpell_Num), Survival_adv := cumprod(1-Hazard_adv), by=list(PerfSpell_Key)]
datCredit[!is.na(PerfSpell_Num), Survival_bas := cumprod(1-Hazard_bas), by=list(PerfSpell_Key)]
# Derive discrete density, or event probability f(t) = S(t-1) . h(t)
datCredit[!is.na(PerfSpell_Num), EventRate_adv := shift(Survival_adv, type="lag", n=1, fill=1) - Survival_adv, by=list(PerfSpell_Key)]
datCredit[!is.na(PerfSpell_Num), EventRate_bas := shift(Survival_bas, type="lag", n=1, fill=1) - Survival_bas, by=list(PerfSpell_Key)]

# - Remove added rows
datCredit <- subset(datCredit, Counter > 0) %>% setorder(LoanID, Date)


# --- Create 12-month PD by summing event probabilities across a rolling window
# NOTE: This step deliberately spans both performing and default spells
# NOTE: Need to specify a (k+1)-window for the "frollapply()" function, e.g., a 12-month outcome implies 13 elements
# Uses the custom function "imputLastKnown" defined in script 0 to fill in NAs caused towards the end of the window
maxDate <- max(datCredit[,Date], na.rm=T) - years(1) # Dates larger than maxDate do not have 12-month default because of the end of the sampling window
datCredit[, PD_12_bas := ifelse(Date<=maxDate,imputeLastKnown(frollapply(
  x=EventRate_bas, n=13, align="left", FUN=function(x){
    sum(x,na.rm=T)
  })),NA), by=list(LoanID)]
datCredit[, PD_12_adv := ifelse(Date<=maxDate,imputeLastKnown(frollapply(
  x=EventRate_adv, n=13, align="left", FUN=function(x){
    sum(x,na.rm=T)
  })),NA), by=list(LoanID)]

# Assign event probability 1 to defaulted cases
# NOTE: This does NOT influence the eventual default rate time graphs
datCredit[is.na(EventRate_bas), PD_12_bas:=1]
datCredit[is.na(EventRate_adv), PD_12_adv:=1]


# --- Graph 12-month expected vs actual default rate

# - Aggregate to reporting time level
datAggr_act <- datCredit[DefaultStatus1==0, list(EventRate = mean(DefaultStatus1_lead_12_max,na.rm=T), Type="a_Act"),, by=list(Date)]
datAggr_bas <- datCredit[DefaultStatus1==0, list(EventRate = mean(PD_12_bas,na.rm=T), Type="b_ExpBasic"),, by=list(Date)]
datAggr_adv <- datCredit[DefaultStatus1==0, list(EventRate = mean(PD_12_adv,na.rm=T), Type="c_ExpAdvanced"),, by=list(Date)]
datAggr <- rbind(datAggr_act, datAggr_bas, datAggr_adv)[Date<=maxDate,]

# - Annotations
(MAE_bas <- mean(abs(datAggr_act$EventRate - datAggr_bas$EventRate), na.rm=T))
(MAE_adv <- mean(abs(datAggr_act$EventRate - datAggr_adv$EventRate), na.rm=T))
datAggr[, FacetLabel := "Discrete-time hazard models"]

# - Graphing parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(9, "Set1")[c(3,1,2)]
vLabel <- c("a_Act" = "Actual rate", 
            "b_ExpBasic" = "Expected rate: Basic",
            "c_ExpAdvanced" = "Expected rate: Advanced")

(gPlot <- ggplot(datAggr, aes(x=Date, y=EventRate, group=Type)) + theme_minimal() + 
    labs(x="Reporting time (ccyymm)", y="12-month default rate (%)") +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Type, linetype=Type), linewidth=0.5) + 
    #geom_point(aes(colour=Type, shape=Type), size=1) + 
    #annotations
    annotate(geom="text", y=0.05, x=as.Date("2015-04-30"), family=chosenFont, size=3,
             label=paste0("MAE (basic): ", round(MAE_bas*100,digits=3), "%")) + 
    annotate(geom="text", y=0.04, x=as.Date("2014-12-31"), family=chosenFont, size=3,
             label=paste0("MAE (advanced): ", round(MAE_adv*100,digits=3), "%")) + 
    # scale options
    facet_grid(FacetLabel ~ .) +
    scale_colour_manual(name="", values=vCol, labels=vLabel) + 
    scale_shape_discrete(name="", labels=vLabel) + 
    scale_linetype_discrete(name="", labels=vLabel) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# Save graph
dpi <- 240
ggsave(gPlot, file=paste0(genFigPath, "DefaultRate_ActVsExp_CoxDisc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")





# ----------------- 4. Actual vs expected 12-month default rates | Cox Proportional Hazards models

# --- Preliminaries
# - Confirm prepared datasets are loaded into memory
### NOTE: Due to computation constraints, we rely on previously scored datasets.
if (!exists('datSurv_PWPST')) unpack.ffdf(paste0(genPath,"datSurv_PWPST"), tempPath);gc()
if (!exists('datSurv_PWPST_basic')) unpack.ffdf(paste0(genPath,"datSurv_PWPST_basic"), tempPath);gc()

# - Merge back with validation set
datCredit <- merge(datCredit_valid_PWPST, 
                   datSurv_PWPST_basic[,list(PerfSpell_Key, TimeInPerfSpell=End, EventProb_bas=EventProb)],
                   by=c("PerfSpell_Key", "TimeInPerfSpell"), all.x=T)
datCredit <- merge(datCredit, 
                    datSurv_PWPST[,list(PerfSpell_Key, TimeInPerfSpell=End, EventProb_adv=EventProb)],
                   by=c("PerfSpell_Key", "TimeInPerfSpell"), all.x=T)


# --- Create 12-month PD by summing event probabilities across a rolling window
# NOTE: This step deliberately spans both performing and default spells
# NOTE: Need to specify a (k+1)-window for the "frollapply()" function, e.g., a 12-month outcome implies 13 elements
# Uses the custom function "imputLastKnown" defined in script 0 to fill in NAs caused towards the end of the window
maxDate <- max(datCredit[,Date], na.rm=T) - years(1) # Dates larger than maxDate do not have 12-month default because of the end of the sampling window
datCredit[, PD_12_bas := ifelse(Date<=maxDate,imputeLastKnown(frollapply(
  x=EventProb_bas, n=13, align="left", FUN=function(x){
    sum(x,na.rm=T)
  })),NA), by=list(LoanID)]
datCredit[, PD_12_adv := ifelse(Date<=maxDate,imputeLastKnown(frollapply(
  x=EventProb_adv, n=13, align="left", FUN=function(x){
    sum(x,na.rm=T)
  })),NA), by=list(LoanID)]


# --- Graph 12-month expected vs actual default rate

# - Aggregate to reporting time level
datAggr_act <- datCredit[DefaultStatus1==0, list(EventRate = mean(DefaultStatus1_lead_12_max,na.rm=T), Type="a_Act"),, by=list(Date)]
datAggr_bas <- datCredit[DefaultStatus1==0, list(EventRate = mean(PD_12_bas,na.rm=T), Type="b_ExpBasic"),, by=list(Date)]
datAggr_adv <- datCredit[DefaultStatus1==0, list(EventRate = mean(PD_12_adv,na.rm=T), Type="c_ExpAdvanced"),, by=list(Date)]
datAggr <- rbind(datAggr_act, datAggr_bas, datAggr_adv)[Date<=maxDate,]

# - Annotations
(MAE_bas <- mean(abs(datAggr_act$EventRate - datAggr_bas$EventRate), na.rm=T))
(MAE_adv <- mean(abs(datAggr_act$EventRate - datAggr_adv$EventRate), na.rm=T))
datAggr[, FacetLabel := "Cox Proportional Hazard (CPH) models"]

# - Graphing parameters
chosenFont <- "Cambria"
vCol <- brewer.pal(9, "Set1")[c(3,1,2)]
vLabel <- c("a_Act" = "Actual rate", 
            "b_ExpBasic" = "Expected rate: Basic",
            "c_ExpAdvanced" = "Expected rate: Advanced")

(gPlot <- ggplot(datAggr, aes(x=Date, y=EventRate, group=Type)) + theme_minimal() + 
    labs(x="Reporting time (ccyymm)", y="12-month default rate (%)") +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Type, linetype=Type), linewidth=0.5) + 
    #geom_point(aes(colour=Type, shape=Type), size=1) + 
    #annotations
    annotate(geom="text", y=0.05, x=as.Date("2015-04-30"), family=chosenFont, size=3,
             label=paste0("MAE (basic): ", round(MAE_bas*100,digits=3), "%")) + 
    annotate(geom="text", y=0.04, x=as.Date("2014-12-31"), family=chosenFont, size=3,
             label=paste0("MAE (advanced): ", round(MAE_adv*100,digits=3), "%")) + 
    # scale options
    facet_grid(FacetLabel ~ .) +
    scale_colour_manual(name="", values=vCol, labels=vLabel) + 
    scale_shape_discrete(name="", labels=vLabel) + 
    scale_linetype_discrete(name="", labels=vLabel) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# Save graph
dpi <- 240
ggsave(gPlot, file=paste0(genFigPath, "DefaultRate_ActVsExp_CoxPH.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


# --- Cleanup
rm(datAggr_act, datAggr_bas, datAggr_adv, datAdd, datAggr, datCredit, datSurv_PWPST, datSurv_PWPST_basic,
  gPlot, datCredit_train, datCredit_valid)

