# ===================================== DATA FUSION Experiments =====================================
# Comparing the distribution of a (winsorised) variable, taken from the full credit dataset and 
# a subsampled version of said dataset. This script serves to inform on whether feature engineering
# can proceed after subsampling is appleid.
# ---------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller
# ---------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 2d(i).Data_Fusion.R

# -- Inputs:
#   - datCredit_real | prepared credit data from script 2d(i)
#   - various parameters set in the setup script 0
#
# -- Outputs:
#   - A graph for comparing the distribution of (a winsorised version of) [Principal].
# ---------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# --- Loading in the full dataset
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

# --- Feature engineering an input space field
(desc_Full_Raw <- describe(datCredit_real$Principal)) # [Principal] has scale [0.01;208477000], with 507600 at 50% quantile and 850000 at 75% quantile and mean of 648267.
hist(datCredit_real$Principal, breaks='FD'); skewness(datCredit_real$Principal, na.rm = T); datCredit_real[Principal==0,.N]/datCredit_real[,.N] # Distribution is skewed to the right; Skewness = 6.606161; 0% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_real[!is.na(PerfSpell_Key) & !(is.na(Principal) | Principal == ""), Principal], probs = 0.999) # 99.9% quantile = 4 440 000
datCredit_real[, Principal_wins := ifelse(Principal < wins_quant, Principal, wins_quant)] # Note that this is applied to the top 0.1%, this is since the scale is still very large when applying Winsorisation to the top 0.01%/
(desc_Full_Wins <- describe(datCredit_real$Principal_wins)) # [Principal_wins] has scale [0.01;4440000.00] and 0 observations have missing values; with 507600 at 50% quantile and 850000 at 75% quantile and mean of 646190.
hist(datCredit_real$Principal_wins, breaks='FD'); skewness(datCredit_real$Principal_wins, na.rm = T); datCredit_real[Principal_wins==0,.N]/datCredit_real[,.N] # [Principal_wins] is skewed to the right; Skewness = 2.016036; 0% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with shifting)
datCredit_real[, Principal_wins_mm_scaled := scaler(Principal_wins,shift = T)]
(desc_Full_Wins_Scaled <- describe(datCredit_real$Principal_wins_mm_scaled))
hist(datCredit_real$Principal_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile.
###            Winsorised variable scaled using min-max scaling with no shifting.
###            [Principal_wins_mms_scaled] has scale [0;1] and 0 observations have missing values; with 0.11432 at 50% quantile and 0.19144 at 75% quantile and mean of 0.1455.
###            No further treatment necessary

# --- Some feature engineering for applying spell-level clustering
# - Max date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Max_Date := max(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date:= NA]
# - Max date of each default spell (used as a stratifier and for identifying FALSE default spells)
datCredit_real[!is.na(DefSpell_Key), DefSpell_Max_Date := max(Date, na.rm=T), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), DefSpell_Max_Date:= NA]
# - Min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Min_Date := min(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Min_Date:= NA]
# - Min date of each default spell (used as a stratifier)
datCredit_real[!is.na(DefSpell_Key), DefSpell_Min_Date := min(Date, na.rm=T), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), DefSpell_Min_Date:= NA]
# - Creating new performance spell resolution types
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  TRUE ~ "Settled & Other"))


# --- Identifying single observation default spells as they are not taken into account directly through the subsampling/ resampling scheme
datCredit_real[, DefSpell_Exc := F] # Creating a variable for identifying these observations
datCredit_real[DefSpell_Counter==1 & Date==DefSpell_Max_Date, DefSpell_Exc := T] # Identifying FALSE default spells
# - Checking how many single-observation default spells exist and assessing their impact
(check.1 <- datCredit_real[DefSpell_Exc==T, .N] / datCredit_real[DefSpell_Counter==1,.N])
cat(sprintf("%.4f", check.1), "% (", datCredit_real[DefSpell_Exc==T, .N], " of ", datCredit_real[DefSpell_Counter==1,.N], ")", "of default spells are to be excluded and are indirectly taken in to account via the subsampling/ resampling scheme.")

# --- Creating a dataset for graphing the distribution of the feature engineering of [Principal]

# --- Clean up
rm(check.1); gc()




# ------ 2. Subsampling
# --- Preliminaries
# - Field names
stratifiers_Perf <- c("PerfSpell_Max_Date", "PerfSpellResol_Type_Hist2") #First variable should be of type "date"
stratifiers_Def <- c("DefSpell_Max_Date", "DefSpellResol_Type_Hist")
targetVar <- c("DefaultStatus1_lead_12_max") # Field name of the main target (i.e., the 12-month default rate)
CurStatus <- "DefaultStatus1"
resolPerf <- "PerfSpellResol_Type_Hist2" # Field name of performance spell resolution types - first level should be the target event (default)
resolDef <- "DefSpellResol_Type_Hist" # Field name of default spell resolution types
excDef <- "DefSpell_Exc" # Field for identifying single-observation default spells (they are indirectly excluded from the resampling scheme)
clusVar <- "LoanID" # Field name of unique account identifier
clusVar_Perf <- "PerfSpell_Key" # Field name of unique performance spell identifier
clusVar_Def <- "DefSpell_Key" # Field name of unique default spell identifier
timeVar <- "Date" # Field name of time variable
counter <- "Counter" # Field name of counter for loan observations
perfCounter <- "PerfSpell_Counter" # Field name of counter for performance spell observations
defCounter <- "DefSpell_Counter" # Field name of counter for default spell observations
perfMax_timeVar <- "PerfSpell_Max_Date" # Field name of last observed time of the associated performance spell
defMax_timeVar <- "DefSpell_Max_Date" # Field name of last observed time of the associated default spell
# - Other parameters
smp_size <- 65000 # fixed size of downsampled set in terms of the number of unique performance spells
cat(smp_size, " is ", sprintf("%.4f", smp_size/length(unique(datCredit_real[,get(clusVar_Perf)]))*100), "% of all performance spells.")
smp_frac <- 0.7 # sampling fraction for resampling scheme
minStrata_size_Perf <- 20 # Minimum size for performance spells
minStrata_size_Def <- 5 # Minimum size for default spells
# - Implied sampling fraction for downsampling step
smp_perc <- smp_size/length(unique(datCredit_real[,get(clusVar_Perf)]))

# --- Subsampling
# - Obtain the first observation of each performance- and default spells (in doing so we obtain all the unique PerfSpell_Keys and DefSpell_Keys)
dat_keys_perf <- unique(datCredit_real[!is.na(clusVar_Perf), mget(c(clusVar_Perf, stratifiers_Perf))])
dat_keys_def <- unique(datCredit_real[get(excDef)==F & !is.na(get(clusVar_Def)), mget(c(clusVar_Def, stratifiers_Def))])

# - Setting the seed
set.seed(1)
# - Use simple random sampling with the stratifiers to select the performance- and default spell IDs that ought to be in the subsampled dataset
dat_smp_keys_perf <- dat_keys_perf %>% group_by(across(all_of(stratifiers_Perf))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
dat_smp_keys_def <- dat_keys_def %>% group_by(across(all_of(stratifiers_Def))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
# - Create two subsampled datasets from the sampled performance- and default spell IDs
datCredit_smp_perf <- datCredit_real %>% subset(get(clusVar_Perf) %in% dat_smp_keys_perf[, get(clusVar_Perf)])
datCredit_smp_def <- datCredit_real %>% subset(get(clusVar_Def) %in% dat_smp_keys_def[, get(clusVar_Def)])
# - Fuse the subsampled datasets
datCredit_smp <- funion(datCredit_smp_perf, datCredit_smp_def, all = F) # Using a union to concatenate the two training datasets (duplicate rows are removed)
# - Arranging the subsample and removing defExc (this will be recreated for the resampling scheme)
datCredit_smp <- datCredit_smp %>% arrange(get(clusVar), get(timeVar)) %>% setDT()
# - Creating a variable for identifying FALSE performance- and default spells
datCredit_smp <- False_Perf_Def(datCredit_smp, LoanID=clusVar, Date=timeVar, PerfSpellID=clusVar_Perf, DefSpellID=clusVar_Def,
                                Counter=counter, PerfSpell_Counter=perfCounter, DefSpell_Counter=defCounter,
                                PerfSpell_Max_Date = perfMax_timeVar, DefSpell_Max_Date = defMax_timeVar)

# - [CHECK] Proportion of default spells in the full dataset and the subsequent training- and validation datasets
check.2_a <- datCredit_real[DefSpell_Counter==1,.N]/(datCredit_real[PerfSpell_Counter==1, .N] + datCredit_real[DefSpell_Counter==1, .N])
check.2_b <- datCredit_smp[DefSpell_Counter==1, .N]/(datCredit_smp[PerfSpell_Counter==1, .N] + datCredit_smp[DefSpell_Counter==1 & !is.na(DefSpell_Key), .N])
cat(sprintf("%.4f", check.2_a*100), "% of all spells in the full dataset are default spells. \n",
    sprintf("%.4f", check.2_b*100), "% of all spells in the subsampled dataset are performance spells. \n")

# - Identifying single observation default spells as they are not taken into account directly through the resampling scheme
datCredit_smp[, (excDef) := F] # Creating a variable for identifying these observations
datCredit_smp[get(defCounter)==1 & get(timeVar)==get(defMax_timeVar), (excDef) := T] # Identifying FALSE default spells
# Checking how many single observations default spells exist and assessing their impact
(check.3 <- datCredit_smp[get(excDef)==T, .N] / datCredit_smp[get(defCounter)==1,.N])
cat(sprintf("%.4f", check.3), "% (", datCredit_smp[get(excDef)==T, .N], " of ", datCredit_smp[get(defCounter)==1,.N], ")", "of default spells are to be excluded and are indirectly taken into account via the subsampling/ resampling scheme.")

# - Clean up
rm(dat_smp_keys_perf, dat_smp_keys_def,  datCredit_smp_perf, datCredit_smp_def, check.2_a, check.2_b, check.3)


# --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold for PERFORMANCE SPELLS
# - Obtaining the stratum that are below the minimum
selectionVar_smp_Perf <- c(clusVar_Perf, timeVar, stratifiers_Perf)
datStrata_smp_min_Perf <- datCredit_smp[get(perfCounter)==1, mget(selectionVar_smp_Perf)][, list(Freq = .N), by=stratifiers_Perf][Freq<minStrata_size_Perf,]
cat(sum(datStrata_smp_min_Perf[,Freq]), "accounts of ", datCredit_smp[get(perfCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Perf[,Freq])/datCredit_smp[get(perfCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Perf)

# - Conditionally applying the exclusions
if (sum(datStrata_smp_min_Perf[,Freq]) > 0){
  # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
  datCredit_smp_old_n_Perf <- datCredit_smp[!is.na(get(clusVar_Perf)),.N]; datCredit_smp_prior_Perf <- (datCredit_smp[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table())[[2]] # NOTE: Change this indexing number manually for the target variable
  # Initiating a vector which will contain the exclusion IDs
  dat_keys_exc_Perf <- NA
  # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
  for (i in 1:datStrata_smp_min_Perf[,.N]){
    class_type_Perf <- sapply(datStrata_smp_min_Perf[,1:length(stratifiers_Perf)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
    
    excCond_Perf <- datStrata_smp_min_Perf[i,1:length(stratifiers_Perf)] # Getting the values of the ith minimum strata
    excCond_Perf <- data.table(Stratifier = colnames(excCond_Perf), # Building a dataset
                               Value = unname(t(excCond_Perf)),
                               Class = class_type_Perf)
    excCond_Perf[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
    excCond_Perf[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
    excCond2_Perf <- parse(text = paste0(paste0(excCond_Perf$Condition, collapse = ""), perfCounter, "==1")) # Compiling the ith condition
    
    dat_keys_exc_Perf <- c(dat_keys_exc_Perf, as.vector(datCredit_smp[eval(excCond2_Perf), get(clusVar_Perf)]))
  }
  dat_keys_exc_Perf <- dat_keys_exc_Perf[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
  
  # Applying the exclusions to the subsampled dataset
  datCredit_smp <- copy(datCredit_smp[!(get(clusVar_Perf) %in% dat_keys_exc_Perf),])
  
  cat(datCredit_smp_old_n_Perf-datCredit_smp[!is.na(get(clusVar_Perf)),.N], " observations removed (", sprintf("%.4f", (datCredit_smp_old_n_Perf-datCredit_smp[!is.na(get(clusVar_Perf)),.N])/datCredit_smp_old_n_Perf*100), "% ) \n",
      "Prior probability = ", sprintf("%.4f", datCredit_smp_prior_Perf*100), "% comapred to ", sprintf("%.4f", (datCredit_smp[get(perfCounter)==1, get(resolPerf)] %>% table() %>% prop.table())[[2]]*100), "%")
}
# - Obtaining the stratum that are below the minimum
datStrata_smp_min_Perf <- datCredit_smp[get(perfCounter)==1, mget(selectionVar_smp_Perf)][, list(Freq = .N), by=stratifiers_Perf][Freq<minStrata_size_Perf,]
cat(sum(datStrata_smp_min_Perf[,Freq]), "accounts of ", datCredit_smp[get(perfCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Perf[,Freq])/datCredit_smp[get(perfCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Perf)

# - Cleanup
suppressWarnings(rm(selectionVar_smp_Perf, datStrata_smp_min_Perf, datCredit_smp_old_n_Perf, datCredit_smp_prior_Perf, class_type_Perf, excCond_Perf, excCond2_Perf))


# --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold for DEFAULT SPELLS
# - Obtaining the stratum that are below the minimum
selectionVar_smp_Def <- c(clusVar_Def, timeVar, stratifiers_Def)
datStrata_smp_min_Def <- datCredit_smp[get(defCounter)==1, mget(selectionVar_smp_Def)][, list(Freq = .N), by=stratifiers_Def][Freq<minStrata_size_Def,]
cat(sum(datStrata_smp_min_Def[,Freq]), "accounts of ", datCredit_smp[get(defCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Def[,Freq])/datCredit_smp[get(defCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Def)

# - Conditionally applying the exclusions
if (sum(datStrata_smp_min_Def[,Freq]) > 0){
  # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
  datCredit_smp_old_n_Def <- datCredit_smp[!is.na(get(clusVar_Def)),.N]; datCredit_smp_prior_Def <- (datCredit_smp[get(defCounter)==1, get(resolDef)] %>% table() %>% prop.table())[[2]] # NOTE: Change this indexing number manually for the target variable
  # Initiating a vector which will contain the exclusion IDs
  dat_keys_exc_Def <- NA
  # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
  for (i in 1:datStrata_smp_min_Def[,.N]){
    class_type_Def <- sapply(datStrata_smp_min_Def[,1:length(stratifiers_Def)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
    
    excCond_Def <- datStrata_smp_min_Def[i,1:length(stratifiers_Def)] # Getting the values of the ith minimum strata
    excCond_Def <- data.table(Stratifier = colnames(excCond_Def), # Building a dataset
                              Value = unname(t(excCond_Def)),
                              Class = class_type_Def)
    excCond_Def[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
    excCond_Def[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
    excCond2_Def <- parse(text = paste0(paste0(excCond_Def$Condition, collapse = ""), defCounter, "==1")) # Compiling the ith condition
    
    dat_keys_exc_Def <- c(dat_keys_exc_Def, as.vector(datCredit_smp[eval(excCond2_Def), get(clusVar_Def)]))
  }
  dat_keys_exc_Def <- dat_keys_exc_Def[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
  
  # Applying the exclusions to the subsampled dataset
  datCredit_smp <- copy(datCredit_smp[!(get(clusVar_Def) %in% dat_keys_exc_Def),])
  
  cat(datCredit_smp_old_n_Def-datCredit_smp[!is.na(get(clusVar_Def)),.N], " observations removed (", sprintf("%.4f", (datCredit_smp_old_n_Def-datCredit_smp[!is.na(get(clusVar_Def)),.N])/datCredit_smp_old_n_Def*100), "% ) \n",
      "Prior probability = ", sprintf("%.4f", datCredit_smp_prior_Def*100), "% comapred to ", sprintf("%.4f", (datCredit_smp[get(defCounter)==1, get(resolDef)] %>% table() %>% prop.table())[[2]]*100), "%")
}
# - Obtaining the stratum that are below the minimum
datStrata_smp_min_Def <- datCredit_smp[get(defCounter)==1, mget(selectionVar_smp_Def)][, list(Freq = .N), by=stratifiers_Def][Freq<minStrata_size_Def,]
cat(sum(datStrata_smp_min_Def[,Freq]), "accounts of ", datCredit_smp[get(defCounter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min_Def[,Freq])/datCredit_smp[!is.na(get(clusVar_Def)) & get(defCounter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size_Def)

# - Clean up
suppressWarnings(rm(selectionVar_smp_Def, datStrata_smp_min_Def, datCredit_smp_old_n_Def, datCredit_smp_prior_Def, class_type_Def, excCond_Def, excCond2_Def,
                    dat_keys_smp_perf, dat_keys_smp_def, dat_keys_exc_Def, dat_keys_exc_Perf, defCounter, defMax_timeVar, CurStatus, counter, excDef,
                    dat_keys_def, dat_keys_perf, datCredit_smp, datStrata_smp_min_Def, datStrata_smp_min_Perf, excCond_Def,
                    smp_frac, smp_perc, resolDef, resolPerf, perfCounter, perfMax_timeVar, timeVar, targetVar, minStrata_size_Def, minStrata_size_Perf,
                    clusVar, clusVar_Def, clusVar_Perf, class_type_Def))




# ------ 3. Feature-engineering on the subsampled dataset
(desc_Sub_Raw <- describe(datCredit_smp$Principal)) # [Principal] has scale [0.01;208477000], with 509000 at 50% quantile and 850000 at 75% quantile and mean of 646260.
hist(datCredit_smp$Principal, breaks="FD"); skewness(datCredit_smp$Principal, na.rm = T); datCredit_smp[Principal==0,.N]/datCredit_smp[,.N] # Distribution is skewed to the right; Skewness = 15.13721; 0% of variables have zero values.
# [TREATMENT] Apply Winsorisation to the upper 0.1% quantile to remove the influence of outliers (which for this variable does not seem to be influential values)
wins_quant <- quantile(datCredit_smp[!is.na(PerfSpell_Key) & !(is.na(Principal) | Principal == ""), Principal], probs = 0.999) # 99.9% quantile = 4 680 480
datCredit_smp[, Principal_wins := ifelse(Principal < wins_quant, Principal, wins_quant)] # Note that this is applied to the top 0.1%, this is since the scale is still very large when applying Winsorisation to the top 0.01%/
(desc_Sub_Wins <- describe(datCredit_smp$Principal_wins)) # [Principal_wins] has scale [0.62;4750000.00] and 0 observations have missing values; with 510000 at 50% quantile and 850000 at 75% quantile and mean of 644103.
hist(datCredit_smp$Principal_wins, breaks='FD'); skewness(datCredit_smp$Principal_wins, na.rm = T); datCredit_smp[Principal_wins==0,.N]/datCredit_smp[,.N] # [Principal_wins] is skewed to the right; Skewness = 2.028359; 0% of variables have zero values.
# [TREATMENT] Scale the data using min-max scaling (with shifting)
datCredit_smp[, Principal_wins_mm_scaled := scaler(Principal_wins,shift = T)]
(desc_Sub_Wins_Scaled <- describe(datCredit_smp$Principal_wins_mm_scaled))
hist(datCredit_smp$Principal_wins_mm_scaled, breaks='FD')
### RESULTS: ~ Extreme outliers detected and Winsorisation applied to deal with the upper 0.1% quantile.
###            Winsorised variable scaled using min-max scaling with no shifting.
###            [Principal_wins_mms_scaled] has scale [0;1] and 0 observations have missing values; with 0.11486 at 50% quantile and 0.19144 at 75% quantile and mean of 0.1451.
###            No further treatment necessary




# ------ 4. Distributional comparison of [Principal]
# --- Creating a long dataset containing both the full- and subsampled datasets
dat_graph <- rbind(datCredit_real %>% subset(select=c("LoanID", "Date", "Principal", "Principal_wins", "Principal_wins_mm_scaled")) %>% mutate(Sample = "Full"),
                   datCredit_smp %>% subset(select=c("LoanID", "Date", "Principal", "Principal_wins", "Principal_wins_mm_scaled")) %>% mutate(Sample = "Subsample"))

# --- Graphing
# - Graphing parameters
chosenFont <- "Cambri"; dpi <- 340
col.v <- brewer.pal(8, "Dark2")[c(1,7)]
fill.v <- brewer.pal(8, "Set2")[c(1,7)]
linetype.v <- c("solid", "dashed")
label.v <- c("Full", "Subsample")

# - Graph for winsorised variable
(g1 <- ggplot(dat_graph, aes(x=Principal_wins, group=Sample)) +
    theme_minimal() +
    theme(text=element_text(family=chosenFont),legend.position = "bottom", 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text = element_text(size=8, colour="gray50"), strip.text.y.right = element_text(angle=90)) +
    geom_histogram(aes(y=..density.., colour=Sample, fill=Sample), position="identity", alpha=0.5, size=0.2, bins = 50) +
    geom_density(aes(colour=Sample, linetype=Sample), size=0.8) +
    scale_colour_manual(name="Sample", values=col.v, labels=label.v) + 
    scale_fill_manual(name="Sample", values=col.v, labels=label.v) + 
    scale_linetype_manual(name="Sample", values=linetype.v, labels=label.v))
### RESULTS:~ Visually, the distributions look the same.
# - Save graph
ggsave(g1, file=paste0(genFigPath, "/Principal_wins_Density_Comparison_", round(datCredit_real[,.N]/1000000),"m.png"), width=2500/(dpi/2), height=2000/(dpi/1.5), dpi=dpi, bg="white")
# - Basic statistics
desc_Full_Wins; desc_Sub_Wins
### RESULTS:~  25% Quantile = 259 000 || 257 000; 50% Quantile = 507 600 || 510 000;75% Quantile = || Mean = 646 190 || 644 103
# - KS-test
ks.test(dat_graph[Sample=="Full", Principal_wins], dat_graph[Sample=="Subsample", Principal_wins]) # P-value = 0
### RESULTS:~ According to a two-sample KS-test, the distributions are not the same; although non-significant p-values in these big sample sizes are prevalent in almost all project/ scripts,
###            we do not put too must trust in them.

### CONSLUCION:~ SAFE to subsample and then use outlier handling, transformations, and scaling.

# --- Clean up
rm(wins_quant, desc_Full_Raw, desc_Full_Wins, desc_Full_Wins_Scaled, desc_Sub_Raw, desc_Sub_Wins, desc_Sub_Wins_Scaled, dat_graph, chosenFont, col.v, fill.v, linetype.v, label.v, g1)

