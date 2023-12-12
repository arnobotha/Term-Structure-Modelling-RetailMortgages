# ============================= RESAMPLING SCHEMES FOR SURVIVAL MODELS - SCS =================================
# A tool for implementing subsampling & resampling on a survival (panel) dataset using simple clustered
# sampling, with either random sampling or n-way stratified sampling. In the case of n-way stratified sampling,
# a frequency analysis is conducted where observations/ records from stratum who have sizes smaller than a
# specified threshold are excluded.
# ------------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Marcel Muller

# DESCRIPTION:
# This ancillary & exploratory script implements a given subsampling- and resampling scheme using simple
# clustered sampling. Either random sampling or n-way stratified sampling is applied as the inner sampling
# technique. For the latter a frequency analysis is conducted and stratum with sizes below a specified
# threshold are excluded. The resulting cross-validation scheme is screened for time-dependent sampling biased
# using the MAEs from a 12-month conditional default rate. In addition to screening for time-dependent
# sampling bias, resolution rates defined for spell entry time is used to investigate the structure of the
# datasets. These graphs are inadequate when investigating the TTC-means though, and subsequent resolution
# rates for spell exit (stop) are with associated MAEs are also computed. These three measures,
# (12-month conditional default rate, resolution rate for spell entry, and resolution rate for spell exit)
# give conclusive insight on the representativeness of the chosen cross-validation scheme and chosen
# inner sampling technique.
# ------------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup
#   - 2d.Data_Fusion

# -- Inputs:
#   - datCredit_real | Prepared from script 2d.Data_Fusion.
#
# -- Outputs:
#   - Series of graphs for testing the time-dependent sampling bias of the chosen (simple clustered sampling)
#     resampling scheme and inner sampling technique:
#       - Frequency analysis of strata size graph (stratification specified as the inner sampling technique);
#       - 12-month conditional default rate graph;
#       - Resolution rates for spell entry (split for performance- and default spells);
#       - Resolution rates for spell exit/ stop (split for performance- and default spells).
# ------------------------------------------------------------------------------------------------------------




# ------ 1. Preliminaries
# --- Load in Dataset
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

# --- Some feature engineering
# - Creating a variable for the first observation of a loan (used as stratification variable)
datCredit_real[, Date_First := Date[1], by=LoanID]
# [SANITY CHECK] Checking if the variable was created correctly
(check.1 <- datCredit_real[is.na(Date_First),.N] == 0) # Should be TRUE
cat(check.1 %?% 'SAFE: variable [Date_First] was successfully created.\n' %:% 'WARNING: variable [Date_First] was not successfully created!\n')

# - Loan ever defaulted
# Creating a default spell "maximum" date as to facilitate subsequent feature engineering
datCredit_real[!is.na(DefSpell_Key), DefSpell_Max_Date := max(Date, na.rm=T), by=list(DefSpell_Key)]
datCredit_real[is.na(DefSpell_Key), DefSpell_Max_Date:= NA]
# Feature engineering the variable
datCredit_real[, HasDefault_Ever := ifelse(max(DefSpell_Num)>0, "Has default(s)", "Has no default(s)"), by=list(LoanID)]

# - Creating a curing indicator for default spells (enables graphing of curing events in default spells)
datCredit_real[, Cured_Ind := ifelse(!is.na(DefSpell_Key) & DefSpellResol_Type_Hist=="Cured" & Date==DefSpell_Max_Date,1,0)] # | Reduce conditions to only the second and third (check feasibility)

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  PerfSpellResol_Type_Hist %in% c("Settled", "Paid-up", "Written-off") ~ "Settled & Other",
                                                                                  TRUE ~ NA))
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist3 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist %in% c("Censored", "Settled", "Paid-up", "Written-off") ~ "Other",
                                                                                  TRUE ~ NA))
# Sanity check - Should be TRUE
datCredit_real[is.na(PerfSpell_Key),.N] == datCredit_real[is.na(PerfSpellResol_Type_Hist2),.N] # TRUE, field created successfully
datCredit_real[is.na(PerfSpell_Key),.N] == datCredit_real[is.na(PerfSpellResol_Type_Hist3),.N] # TRUE, field created successfully


# - Creating new default spell resolution types
datCredit_real <- datCredit_real %>% mutate(DefSpellResol_Type_Hist2 = case_when(DefSpellResol_Type_Hist=="WOFF" ~ "WOFF",
                                                                                 DefSpellResol_Type_Hist %in% c("Cured", "Censored") ~ "Other",
                                                                                 TRUE ~ NA))
# Sanity check - Should be TRUE
datCredit_real[is.na(DefSpell_Key),.N] == datCredit_real[is.na(DefSpellResol_Type_Hist2),.N] # TRUE, field created successfully


# --- Field specification and subsetting
# - Confidence interval parameter
confLevel <- 0.95

# - Required field names
targetVar <- c("DefaultStatus1_lead_12_max") # Field name of the main target (i.e., the 12-month default rate)
CurStatus <- "DefaultStatus1" # Field name of the current status of an account (default vs non-default)
resolPerf_targetVar <- "Defaulted" # Reference level in the performance spell resolution type (as specified by [resolPerf_start]) for the target variable
clusVar <- "LoanID"
clusVar_Perf <- "PerfSpell_Key"
clusVar_Def <- "DefSpell_Key"
timeVar <- "Date"
counter <- "Counter"
PerfSpell_counter <- "PerfSpell_Counter"
DefSpell_counter <- "DefSpell_Counter"
# - For stratification (optional)
stratifiers <- NA # c("Date_First", "Event_Type") # First variable should be of type "date" | Assign "NA" for no stratifiers | Other good stratifier candidates are "Event_Type", "LN_TPE", and HasDefaulted_Ever
# Facet specification (for graphing purposes of the resolution rates)
resolPerf <- "PerfSpellResol_Type_Hist2" # Field name of performance spell resolution types - first level should be the target event (default)
resolDef <- "DefSpellResol_Type_Hist" # Field name of default spell resolution types
resolPerf_stop <- "PerfSpellResol_Type_Hist3" # Field name for performance spell resolution rate; specific for the stopping time cohort | Assign [resolPerf] if not interested in controlling the resolution rate facet for stop dates
resolPerf_stop2 <- "Defaulted" # Name of the main resolution type (typically default) used to obtain a single facet (this level needs to be within the resolPerf_stop variable) | Set to NA if not interested in creating additional facets for the performance spells using stopping time
resolDef_stop <-  "DefSpellResol_Type_Hist2" # Field name for default spell resolution rate; specific for the stopping time cohort | Assign [resolDef] if not interested in controlling the resolution rate facet for stop dates
# - Final selection
selectionVar <- unique(c(clusVar, clusVar_Perf, clusVar_Def, timeVar, counter,
                         resolPerf, resolDef, resolPerf_stop, resolDef_stop, CurStatus, stratifiers, targetVar)) # Variables to subset
selectionVar <- selectionVar[!is.na(selectionVar)] # Facilitating cases where the variables are left unspecified (specifically for use of no stratifiers)


# - Subset given dataset accordingly; an efficiency enhancement
datCredit <- subset(datCredit_real, select=selectionVar)

# - Subsampling & resampling parameters
smp_size <- 65000 # fixed size of downsampled set in terms of the number of unique loan accounts
cat(smp_size, " is ", percent(smp_size/length(unique(datCredit[,get(clusVar)]))), "of all loans.")
smp_frac <- 0.7 # sampling fraction for resampling scheme

minStrata_size <- 0 # Minimum strata size specified for subsample


# --- Clean up
rm(datCredit_real, check.1); gc()




# ------ 2. Subsampled resampling scheme: basic cross-validation with simple random sampling
# --- Feature engineering spell level date variables
# - Max- and min date of each performance spell (used as a stratifier)
datCredit[!is.na(get(clusVar_Perf)), c("timeVar_Perf_Min","timeVar_Perf_Max") := as.list(range(get(timeVar), na.rm=TRUE)), by=list(get(clusVar_Perf))]
datCredit[is.na(get(clusVar_Perf)), c("timeVar_Perf_Min","timeVar_Perf_Max") := as.list(c(NA,NA)), by=list(get(clusVar_Perf))]
# - Max- and min date of each default spell (used as a stratifier and for identifying FALSE default spells)
datCredit[!is.na(get(clusVar_Def)), c("timeVar_Def_Min","timeVar_Def_Max") := as.list(range(get(timeVar), na.rm=TRUE)), by=list(get(clusVar_Def))]
datCredit[is.na(get(clusVar_Def)), c("timeVar_Def_Min","timeVar_Def_Max") := as.list(c(NA,NA)), by=list(get(clusVar_Def))]

# --- Preliminaries
# - Implied sampling fraction for downsampling step
smp_perc <- smp_size/length(unique(datCredit[,get(clusVar)]))

# --- Downsample data into a set with a fixed size (using stratified sampling) before implementing resampling scheme
# - Set seed
set.seed(1)
# - Conditional loop for stratifiers
if (all(is.na(stratifiers))){ # - No stratifiers
  # Get unique loan account IDs from the full dataset
  dat_keys <- unique(datCredit[, mget(c(clusVar))])
  # Use simple random sampling to select the loan IDs that ought to be in the subsampled dataset
  dat_smp_keys <- dat_keys %>% slice_sample(prop=smp_perc) %>% as.data.table()
} else { # - Stratifiers
  # Get unique loan account IDs from the full dataset
  dat_keys <- unique(datCredit[, mget(c(clusVar, stratifiers))])
  # Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the subsampled dataset
  dat_smp_keys <- dat_keys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_perc) %>% as.data.table()
}
# - Obtain the associated loan records as to create the subsampled dataset
datCredit_smp <- copy(datCredit[get(clusVar) %in% dat_smp_keys[, get(clusVar),]])


# --- Minimum stratum analysis and subsequent exclusions to ensure adherence to specified threshold
# - Obtaining the stratum that are below the minimum
if (all(!is.na(stratifiers))){ # - Conditional loop for strata
  selectionVar_smp <- c(clusVar, timeVar, stratifiers)
  datStrata_smp_min <- datCredit_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
  cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datCredit_smp[get(counter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datCredit_smp[get(counter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
  
  # - Conditionally applying the exclusions
  if (sum(datStrata_smp_min[,Freq]) > 0){
    # Saving the number of records and the prior probability, in the subsampled dataset, for reporting
    datCredit_smp_old_n <- datCredit_smp[,.N]
    datCredit_smp_prior <- datCredit_smp[get(timeVar)==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table(); datCredit_smp_prior <- datCredit_smp_prior[names(datCredit_smp_prior)[names(datCredit_smp_prior) == resolPerf_targetVar]][[1]] # Computing the prior probabilities of the performance spell resolution outcomes
    # Initiating a vector which will contain the exclusion IDs
    dat_keys_exc <- NA
    # Looping through the minimum strata dataset and building an exclusion condition (filter) for each row therein
    for (i in 1:datStrata_smp_min[,.N]){
      class_type <- sapply(datStrata_smp_min[,1:length(stratifiers)], function(x) {class(x[[1]])}) # Getting the type of class of each stratifier (used for building the ith condition)
      
      excCond <- datStrata_smp_min[i,1:length(stratifiers)] # Getting the values of the ith minimum strata
      excCond <- data.table(Stratifier = colnames(excCond), # Building a dataset
                            Value = unname(t(excCond)), # Ensure that the column name is Value instead of Value.V1
                            Class = class_type)
      excCond[, Value.V1 := ifelse(Class %in% c("numeric", "Date"), paste0("as.",Class,"(",'"',Value.V1,'"',")"), paste0('"', Value.V1, '"'))]
      excCond[, Condition := paste0(Stratifier, " == ", Value.V1, " & ")] # Adding an "and" operator to enable multiple conditions
      excCond2 <- parse(text = paste0(paste0(excCond$Condition, collapse = ""), counter,"==1")) # Compiling the ith condition
    
      dat_keys_exc <- c(dat_keys_exc, as.vector(datCredit_smp[eval(excCond2), get(clusVar)]))
    }
    dat_keys_exc <- dat_keys_exc[-1] # Removing the first value (as it is a missing value stemming from the vector's creation)
    
    # Applying the exclusions to the subsampled dataset
    datCredit_smp <- copy(datCredit_smp[!(get(clusVar) %in% dat_keys_exc),])
    
    cat(datCredit_smp_old_n-datCredit_smp[,.N], " observations removed (", sprintf("%.4f", (datCredit_smp_old_n-datCredit_smp[,.N])/datCredit_smp_old_n*100), "% ) \n",
        "Prior probability = ", sprintf("%.4f", datCredit_smp_prior*100), "% comapred to ", sprintf("%.4f", (datCredit_smp[get(timeVar)==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table())[[2]]*100), "%")
  }
  # - Obtaining the stratum that are below the minimum
  datStrata_smp_min <- datCredit_smp[get(counter)==1, mget(selectionVar_smp)][, list(Freq = .N), by=stratifiers][Freq<minStrata_size,]
  cat(sum(datStrata_smp_min[,Freq]), "accounts of ", datCredit_smp[get(counter)==1,.N], "(", sprintf("%.4f", sum(datStrata_smp_min[,Freq])/datCredit_smp[get(counter)==1,.N]*100), "%) need to be excluded to ensure a minimum strata size of ", minStrata_size)
}



# --- Implement resampling scheme using given main sampling fraction
# - Set seed
set.seed(1)
# - Use simple random sampling with the stratifiers to select the loan IDs that ought to be in the training dataset
if (all(!is.na(stratifiers))){ # Stratifiers
  dat_train_keys <- dat_smp_keys %>% group_by(across(all_of(stratifiers))) %>% slice_sample(prop=smp_frac) %>% as.data.table() 
} else { # No stratifiers
  dat_train_keys <- dat_smp_keys %>% slice_sample(prop=smp_frac) %>% as.data.table()
}
# - Obtain the associated loan records as to create the training dataset
datCredit_train <- copy(datCredit_smp[get(clusVar) %in% dat_train_keys[, get(clusVar)],])
# - Obtain the associated loan records of the validation dataset
datCredit_valid <- copy(datCredit_smp[!(get(clusVar) %in% dat_train_keys[, get(clusVar)]),])

# --- [SANITY CHECK]
(check.2 <- datCredit_smp[,.N] == datCredit_train[,.N] + datCredit_valid[,.N]) # Should be TRUE

# --- Clean up
suppressWarnings(rm(smp_perc, dat_keys, dat_smp_keys, dat_train_keys, check.2, datCredit_smp_old_n, datCredit_smp_prior, dat_keys_exc, class_type, excCond, excCond2))




# ------ 3. Preliminaries to graphing
# --- Restructuring the datasets (renaming fields as to use standardised naming convention)
# - Full dataset
colnames(datCredit)[1:5] <- c("ClusVar", "ClusVar_Perf", "ClusVar_Def", "timeVar", "Counter")
# - Subsampled dataset
colnames(datCredit_smp)[1:5] <- c("ClusVar", "ClusVar_Perf", "ClusVar_Def", "timeVar", "Counter")
# - Training dataset
colnames(datCredit_train)[1:5] <- c("ClusVar", "ClusVar_Perf", "ClusVar_Def", "timeVar", "Counter")
# - Validation dataset
colnames(datCredit_valid)[1:5] <- c("ClusVar", "ClusVar_Perf", "ClusVar_Def", "timeVar", "Counter")

# --- Merge datasets together for graphing purposes
datGraph <- rbind(datCredit[, Sample:="a_Full"],
                  datCredit_train[, Sample:="b_Train"],
                  datCredit_valid[, Sample:="c_Valid"])




# ------ 4. Frequency analysis on sampling technique
if (all(!is.na(stratifiers))){
  # - Determine subsampling window given cross-sectional design
  StartDte <- min(datCredit_smp$timeVar, na.rm=T)
  EndDte <- max(datCredit_smp$timeVar, na.rm=T)
  minDate <- StartDte %m+% months(1)
  maxDate <- EndDte - years(1)# A post-hoc filter, used for graphing purposes, given a 12-month outcome window
  
  # - Aggregate data according to the same n-way stratified sampling technique used within subsampling/resampling scheme
  selectionVar_train <- c("ClusVar", "timeVar", stratifiers)
  datStrata <- copy(datCredit_train[Counter==1, ..selectionVar_train][, list(Freq = .N), by=stratifiers])
  
  # - Aesthetics engineering
  datStrata[, Facet_label := "Strata Frequency Analysis"]
  
  # - Create summaries for annotations within graph
  datStrata_aggr <- datStrata[, list(StratumSize_N = .N, StratumSize_Min = min(Freq,na.rm=T), StratumSize_Mean = mean(Freq,na.rm=T),
                                     StratumSize_SD = sd(Freq,na.rm=T))]
  
  datStrata_aggr[, StrataSize_Margin := qnorm(1-(1-confLevel)/2) * datStrata_aggr$StratumSize_SD / sqrt(datStrata_aggr$StratumSize_N)]
  
  # - Graphing parameters
  chosenFont <- "Cambria"; dpi <- 340
  col.v <- brewer.pal(8, "Dark2") # [1:length(targetVars_Perf)]
  fill.v <- brewer.pal(8, "Set2") # [1:length(targetVars_Perf)]
  x_pos <- min(as.vector(data.table(datStrata[,..stratifiers])[,1])[[1]]) + round((max(as.vector(data.table(datStrata[,..stratifiers])[,1])[[1]]) - min(as.vector(data.table(datStrata[,..stratifiers])[,1])[[1]]))/2) # x-position fot the annotation
  
  # - Create graph to evidence minimum strata sizes
  (g1 <- ggplot(datStrata[get(stratifiers[1])>=minDate & get(stratifiers[1])<=maxDate,], aes(x=get(stratifiers[1]), y=Freq)) + theme_minimal() + 
      labs(x=bquote("Date "*italic(t)), y=bquote("Proporionate volume stratifiers (%) within "*italic(D[T])~"("*.(round(datCredit_train[,.N]/1000))*"k)")) + 
      theme(text=element_text(family=chosenFont),legend.position = "bottom",
            axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
            strip.background=element_rect(fill="snow2", colour="snow2"),
            strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) +
      # Creating conditional geoms based on the number of stratification variables used
      {if (length(stratifiers)>1){
        # main bar graph
        list(geom_bar(position="stack", stat="identity", aes(colour=get(stratifiers[-1]), fill=get(stratifiers[-1]))),
        # scale options
             scale_colour_manual(name=stratifiers[-1], values=col.v),
             scale_fill_manual(name=stratifiers[-1], values=fill.v))
      } else {
       # main bar graph
        list(geom_bar(position="stack", stat="identity", colour=col.v[1], fill=fill.v[1]))
          }} +
      # annotations
      annotate("text", x=x_pos, y=Inf, size=3, hjust=0.5, vjust=4, family=chosenFont,
                 label=paste0(datStrata_aggr$StratumSize_N, " total strata with a mean cell size of ", 
                              comma(datStrata_aggr$StratumSize_Mean, accuracy=0.1),
                              " ± ", sprintf("%.1f", datStrata_aggr$StrataSize_Margin), " and a minimum size of ", 
                              sprintf("%.0f", datStrata_aggr$StratumSize_Min))) + 
        # Rest of the facet & scale options
        facet_grid(Facet_label ~ .) + 
        scale_y_continuous(breaks=pretty_breaks(), label=comma) + 
        scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))
  
  # - Save graph
  ggsave(g1, file=paste0(genFigPath_Res_anc, "StrataDesign_Train_", round(datCredit_smp[,.N]/1000),"k.png"), width=2550/dpi, height=2000/dpi, dpi=dpi, bg="white")
  
  # --- Clean up
  rm(datStrata, datStrata_aggr, selectionVar_train, col.v, fill.v, chosenFont, dpi)
}




# ------ 5. Graphing the main event rate (12-month conditional default rate) over time given resampled subsample
# - Check representatives | dataset-level proportions should be similar
datCredit[get(CurStatus)==0, ..targetVar] %>% table() %>% prop.table()
datCredit_train[get(CurStatus)==0, ..targetVar] %>% table() %>% prop.table()
datCredit_valid[get(CurStatus)==0, ..targetVar] %>% table() %>% prop.table()

# - Setting some aggregation parameters, purely to facilitate graphing aesthetics
StartDte <- min(datCredit_smp$timeVar, na.rm=T)
EndDte <- max(datCredit_smp$timeVar, na.rm=T)
maxDate <- EndDte - years(1) # A post-hoc filter, used for graphing purposes, given a 12-month outcome window (scenario specific)
minDate <- StartDte

# - Aggregate to monthly level and observe up to given point
# Aggregation for main target event (maximum time imposed to cater for the last 12-months)
port.aggr <- datGraph[get(CurStatus)==0, sum(get(targetVar), na.rm=T)/.N, by=list(Sample, timeVar)][timeVar >= StartDte & timeVar <= maxDate,] %>% setkey(timeVar)
# Renaming the column containing the conditional event rate
colnames(port.aggr)[3] <- targetVar
# Creating a faceting label
port.aggr[, Facet:="Worst-ever aggregation approach"]
# - Calculate MAE over time by sample, by target event
port.aggr2 <- port.aggr %>% pivot_wider(id_cols = c(timeVar), names_from = c(Sample), values_from = targetVar) %>% as.data.table()

# - Initiating the annotation dataset
x_pos <- min(datCredit$timeVar) + round((max(datCredit$timeVar) - min(datCredit$timeVar))/2) # x-position fot the annotation
dat_anno <- data.table(MAE = rep(0,4),
                            Mean_EventRate = rep(0, 4),
                            stdError_EventRate = rep(0, 4),
                            margin_EventRate = rep(0, 4),
                            name = rep(targetVar,4),
                            Dataset = c("A-B","A-C","B-C", "B"),
                            Label = c(paste0("'MAE between '*italic(A[t])*' and '*italic(B[t])*'"),
                                      paste0("'MAE between '*italic(A[t])*' and '*italic(C[t])*'"),
                                      paste0("'MAE between '*italic(B[t])*' and '*italic(C[t])*'"), 
                                      paste0("'TTC-mean '*E(italic(B[t]))*'")),
                            x = rep(x_pos,4),
                            y = rep(Inf, 4), # c(0.06,0.0585, 0.0570, 0.0540),
                            vjust = c(2,3,4,6),
                            hjust=rep(0.5,4))
# - Getting the column names to help compute the MAEs
colnames <- colnames(port.aggr2)
# - Populating the annotation dataset
dat_anno[1, MAE := mean(abs(port.aggr2[,2][[1]] - port.aggr2[,3][[1]]), na.rm = T)] # MAE between the full- and training dataset
dat_anno[2, MAE := mean(abs(port.aggr2[,2][[1]] - port.aggr2[,4][[1]]), na.rm = T)] # MAE between the training- and validation dataset
dat_anno[3, MAE := mean(abs(port.aggr2[,3][[1]] - port.aggr2[,4][[1]]), na.rm = T)] # MAE between the training- and validation dataset
dat_anno[4, mean_EventRate := mean(port.aggr2[,3][[1]], na.rm=T)]
dat_anno[4, stdError_EventRate := sd(port.aggr2[,3][[1]], na.rm=T)/ sqrt(length(port.aggr2[,3][[1]]))]
# - Last adjustments to the annotation dataset before plotting
dat_anno[4, margin_EventRate := qnorm(1-(1-confLevel)/2) * stdError_EventRate]
dat_anno[4, Label := paste0(Label, " = ", sprintf("%.2f", mean_EventRate*100) , "% +-", sprintf("%.3f", margin_EventRate*100), "%'")]
dat_anno[1:3, Label := paste0(Label, " = ", sprintf("%.4f",MAE*100), "%'")]

# - Graphing parameters
chosenFont <- "Cambria"; dpi <- 170
col.v <- brewer.pal(9, "Set1")
label.v <- c("a_Full"=expression(italic(A)[t]*": Full set "*italic(D)),
             "b_Train"=bquote(italic(B)[t]*": Training set "*italic(D)[italic(T)]~"("*.(round(datCredit_train[,.N]/1000))*"k)"),
             "c_Valid"=bquote(italic(C)[t]*": Validation set "*italic(D)[italic(V)]~"("*.(round(datCredit_valid[,.N]/1000))*"k)"))

# - Create graph 1 (all sets)
(g2 <- ggplot(port.aggr[timeVar>=minDate & timeVar<=maxDate,], aes(x=timeVar, y=get(targetVar))) + theme_minimal() + 
    labs(x="Reporting date (months)", y=bquote("Conditional 12-month default rate (%) across sample "*italic(bar(D)))) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Sample, linetype=Sample)) + 
    geom_point(aes(colour=Sample, shape=Sample), size=1) + 
    # facets
    facet_wrap(.~Facet, scales = "free", strip.position="right") + 
    #annotations
    geom_text(data=dat_anno, aes(x=x, y=y, hjust=hjust, vjust=vjust, label = Label), family=chosenFont, size=3, parse=T) + 
    # scale options
    scale_colour_manual(name=bquote("Sample "*italic(bar(D))), values=col.v, labels=label.v) + 
    scale_shape_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + scale_linetype_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# - Save graph
ggsave(g2, file=paste0(genFigPath_Res_anc, "DefaultRate_Subsample-", round(datCredit_smp[,.N]/1000),"k.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")



# - Cleanup
rm(dat_anno, chosenFont, col.v, label.v, colnames, port.aggr, port.aggr2, maxDate, minDate)




# ------ 6. Graphing performance spell resolution rates over time given the cross-validation scheme | Spell entry time (t_e)
# - Check representatives | dataset-level proportions should be similar
datCredit[timeVar==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table()
datCredit_train[timeVar==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table()
datCredit_valid[timeVar==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table()
# Checking the proportions of the subsampled dataset and creating the corresponding faceting labels
(Facet_Label_Perf <- datCredit_smp[timeVar==timeVar_Perf_Min, get(resolPerf)] %>% table() %>% prop.table() %>% data.table()) # Saving these proportions as they are used in the facets
colnames(Facet_Label_Perf) <- c("PerfSpell_Resol", "Prior") # Renaming the columns
# - Subsetting the main long dataset to only include the necessary variables for performance spells
datGraph_Perf <- datGraph %>% subset(!is.na(ClusVar_Perf), select = c("ClusVar_Perf", "timeVar", "timeVar_Perf_Min", resolPerf, "Sample"))
colnames(datGraph_Perf) <- c("ClusVar_Perf", "timeVar", "timeVar_Perf_Min", "PerfSpell_Resol", "Sample")

# - Setting some aggregation parameters, purely to facilitate graphing aesthetics
StartDte <- min(datCredit_smp$timeVar, na.rm=T)
EndDte <- max(datCredit_smp$timeVar, na.rm=T)
maxDate <- EndDte # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
minDate <- StartDte %m+% months(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window

# - Aggregate to monthly level and observe up to given point
port.aggr_perf <- merge(datGraph_Perf[timeVar==timeVar_Perf_Min, list(Sum_Total = .N), by=list(Sample,timeVar)],
                        datGraph_Perf[timeVar==timeVar_Perf_Min, list(Sum_Resol = .N), by=list(Sample,timeVar,PerfSpell_Resol)],
                        by=c("Sample", "timeVar"))[timeVar >= minDate & timeVar <= maxDate,]
port.aggr_perf[, Prop := Sum_Resol/Sum_Total]

# - Calculate MAE over time by sample, by performance event(s)
port.aggr_perf2 <- port.aggr_perf %>% pivot_wider(id_cols = c(timeVar), names_from = c(Sample, PerfSpell_Resol), values_from = Prop) %>% data.table()

# - Get the unique factors of the performance spell resolution type (used as facets)
resolPerf_levels <- unique(datCredit[!is.na(get(resolPerf)), get(resolPerf)])
# - Number of annotations "sets" to create
anno_n <- length(resolPerf_levels)
# - Initiating the annotation dataset
x_pos <- min(datCredit$timeVar) + round((max(datCredit$timeVar) - min(datCredit$timeVar))/2) # x-position fot the annotation
dat_anno_perf <- data.table(MAE = rep(0,anno_n*3),
                            Mean_EventRate = rep(0, anno_n*3),
                            stdError_EventRate = rep(0, anno_n*3),
                            margin_EventRate = rep(0, anno_n*3),
                            PerfSpell_Resol = unlist(lapply(resolPerf_levels, function(x){rep(x,3)})),
                            Dataset = rep(c("A-B","A-C","B-C"), anno_n),
                            Label = rep(c(paste0("'MAE between '*italic(A[t])*' and '*italic(B[t])*'"),
                                          paste0("'MAE between '*italic(A[t])*' and '*italic(C[t])*'"),
                                          paste0("'MAE between '*italic(B[t])*' and '*italic(C[t])*'")),
                                        anno_n),
                            x = rep(x_pos,anno_n*3),
                            y = rep(Inf, anno_n*3), # c(c(0.9,0.83,0.76),c(0.6,0.55,0.5),c(0.9,0.83,0.76)),
                            vjust = rep(c(1,2,3), anno_n),
                            hjust=rep(0.5, anno_n*3))
# - Getting the column names to help compute the MAEs
colnames <- colnames(port.aggr_perf2)
# - Populating the annotation dataset
for (i in 1:anno_n){
  dat_anno_perf[i*3-2, MAE := mean(abs(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("a_Full_", resolPerf_levels[i])])[[1]] - subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and training dataset
  dat_anno_perf[i*3-1, MAE := mean(abs(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("a_Full_", resolPerf_levels[i])])[[1]] - subset(port.aggr_perf2, select=colnames[colnames %in% paste0("c_Valid_", resolPerf_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and validation dataset
  dat_anno_perf[i*3, MAE := mean(abs(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]] - subset(port.aggr_perf2, select=colnames[colnames %in% paste0("c_Valid_", resolPerf_levels[i])])[[1]]), na.rm = T)] # MAE between the training- and validation dataset
}
# - Finessing the annotation dataset for plotting
dat_anno_perf[, Label := paste0(Label, " = ", sprintf("%.4f",MAE*100), "%'")]
# - Adding an column to accommodate the facets
port.aggr_perf <- merge(port.aggr_perf, Facet_Label_Perf, by="PerfSpell_Resol")
port.aggr_perf[, Facet:=paste0('"', PerfSpell_Resol, ' (', sprintf("%.2f", Prior*100), '%)"')]
dat_anno_perf <- merge(dat_anno_perf, unique(subset(port.aggr_perf, select=c("PerfSpell_Resol", "Facet"))), by="PerfSpell_Resol")

# - Graphing parameters
chosenFont <- "Cambria"; dpi <- 340
col.v <- brewer.pal(9, "Set1")
label.v <- c("a_Full"=expression(italic(A)[t]*": Full set "*italic(D)),
             "b_Train"=bquote(italic(B)[t]*": Training set "*italic(D)[italic(T)]~"("*.(round(datCredit_train[,.N]/1000))*"k)"),
             "c_Valid"=bquote(italic(C)[t]*": Validation set "*italic(D)[italic(V)]~"("*.(round(datCredit_valid[,.N]/1000))*"k)"))

# - Create graph
(g3 <- ggplot(port.aggr_perf, aes(x=timeVar, y=Prop)) + theme_minimal() + 
    labs(x=bquote("Performing spell cohorts (ccyymm): entry time "*italic(t[e])), y=bquote("Resolution rate (%) of type "*~italic(kappa))) +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Sample, linetype=Sample)) + 
    geom_point(aes(colour=Sample, shape=Sample), size=1) + 
    # facets
    facet_wrap(Facet~., labeller = label_parsed, scales = "free", nrow=length(resolPerf_levels), strip.position="right") + 
    #annotations
    geom_text(data=dat_anno_perf, aes(x=x, y=y, hjust=hjust, vjust=vjust, label = Label), family=chosenFont, size=3, parse=T) + 
    # scale options
    scale_colour_manual(name=bquote("Sample "*italic(bar(D))), values=col.v, labels=label.v) + 
    scale_shape_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + scale_linetype_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# - Save graph
ggsave(g3, file=paste0(genFigPath_Res_anc, "ResolutionRates_Perf_te_Subsample-", round(datCredit_smp[,.N]/1000),"k.png"), width=5000/(dpi*2.25), height=4000/(dpi*1.4), dpi=dpi, bg="white")

# - Cleanup
rm(dat_anno_perf, resolPerf_levels, chosenFont, col.v, label.v, colnames, datGraph_Perf, port.aggr_perf, port.aggr_perf2, maxDate, minDate, Facet_Label_Perf)




# ------ 7. Graphing default spell resolution rates over time given resampled | Spell entry time (t_e)
# - Check representatives | dataset-level proportions should be similar
datCredit[timeVar==timeVar_Def_Min, get(resolDef)] %>% table() %>% prop.table()
datCredit_train[timeVar==timeVar_Def_Min, get(resolDef)] %>% table() %>% prop.table()
datCredit_valid[timeVar==timeVar_Def_Min, get(resolDef)] %>% table() %>% prop.table()
# Checking the proportions of the subsampled dataset and creating the corresponding faceting labels
(Facet_Label_Def <- datCredit_smp[timeVar==timeVar_Def_Min, get(resolDef)] %>% table() %>% prop.table() %>% data.table()) # Saving these proportions as they are used in the facets
colnames(Facet_Label_Def) <- c("DefSpell_Resol", "Prior") # Renaming the columns

# - Subsetting the main long dataset to only include the necessary variables for performance spells
datGraph_Def <- datGraph %>% subset(!is.na(ClusVar_Def), select = c("ClusVar_Def", "timeVar", "timeVar_Def_Min", "timeVar_Def_Min", resolDef, "Sample"))
colnames(datGraph_Def) <- c("ClusVar_Def", "timeVar", "timeVar_Def_Min", "timeVar_Def_Min", "DefSpell_Resol", "Sample")

# - Setting some aggregation parameters, purely to facilitate graphing aesthetics
StartDte <- min(datCredit_smp$timeVar, na.rm=T)
EndDte <- max(datCredit_smp$timeVar, na.rm=T)
maxDate <- EndDte # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
minDate <- StartDte # %m+% months(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window

# - Aggregate to monthly level and observe up to given point
port.aggr_def <- merge(datGraph_Def[timeVar==timeVar_Def_Min, list(Sum_Total = .N), by=list(Sample,timeVar)],
                       datGraph_Def[timeVar==timeVar_Def_Min, list(Sum_Resol = .N), by=list(Sample,timeVar,DefSpell_Resol)],
                       by=c("Sample", "timeVar"))[timeVar >= minDate & timeVar <= maxDate,]
port.aggr_def[, Prop := Sum_Resol/Sum_Total]

# - Calculate MAE over time by sample, by competing event(s)
port.aggr_def2 <- port.aggr_def %>% pivot_wider(id_cols = c(timeVar), names_from = c(Sample, DefSpell_Resol), values_from = Prop) %>% as.data.table()

# - Get the unique factors of the default spell resolution type (used as facets)
resolDef_levels <- unique(datCredit[!is.na(get(resolDef)), get(resolDef)])
# - Number of annotations "sets" to create
anno_n <- length(resolDef_levels)

# - Initiating the annotation dataset
x_pos <- min(datCredit$timeVar) + round((max(datCredit$timeVar) - min(datCredit$timeVar))/2) # x-position fot the annotation
dat_anno_def <- data.table(MAE = rep(0,anno_n*3),
                           Mean_EventRate = rep(0, anno_n*3),
                           stdError_EventRate = rep(0, anno_n*3),
                           margin_EventRate = rep(0, anno_n*3),
                           DefSpell_Resol = unlist(lapply(resolDef_levels, function(x){rep(x,3)})),
                           Dataset = rep(c("A-B","A-C","B-C"), anno_n),
                           Label = rep(c(paste0("'MAE between '*italic(A[t])*' and '*italic(B[t])*'"),
                                         paste0("'MAE between '*italic(A[t])*' and '*italic(C[t])*'"),
                                         paste0("'MAE between '*italic(B[t])*' and '*italic(C[t])*'")),
                                       anno_n),
                           x = rep(x_pos,anno_n*3),
                           y = rep(Inf, anno_n*3),
                           vjust = rep(c(1,2,3), anno_n),
                           hjust=rep(0.5, anno_n*3))
# - Getting the column names to help compute the MAEs
colnames <- colnames(port.aggr_def2)
# - Populating the annotation dataset
for (i in 1:anno_n){
  dat_anno_def[i*3-2, MAE := mean(abs(subset(port.aggr_def2, select=colnames[colnames %in% paste0("a_Full_", resolDef_levels[i])])[[1]] - subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and training dataset
  dat_anno_def[i*3-1, MAE := mean(abs(subset(port.aggr_def2, select=colnames[colnames %in% paste0("a_Full_", resolDef_levels[i])])[[1]] - subset(port.aggr_def2, select=colnames[colnames %in% paste0("c_Valid_", resolDef_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and validation dataset
  dat_anno_def[i*3, MAE := mean(abs(subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]] - subset(port.aggr_def2, select=colnames[colnames %in% paste0("c_Valid_", resolDef_levels[i])])[[1]]), na.rm = T)] # MAE between the training- and validation dataset
}
# - Finessing the annotation dataset for plotting
dat_anno_def[, Label := paste0(Label, " = ", sprintf("%.4f",MAE*100), "%'")]
# - Adding an column to accommodate the facets
port.aggr_def <- merge(port.aggr_def, Facet_Label_Def, by="DefSpell_Resol")
port.aggr_def[, Facet:=paste0('"', DefSpell_Resol, ' (', sprintf("%.2f", Prior*100), '%)"')]
dat_anno_def <- merge(dat_anno_def, unique(subset(port.aggr_def, select=c("DefSpell_Resol", "Facet"))), by="DefSpell_Resol")


# - Graphing parameters
chosenFont <- "Cambria"; dpi <- 340
col.v <- brewer.pal(9, "Set1")
label.v <- c("a_Full"=expression(italic(A)[t]*": Full set "*italic(D)),
             "b_Train"=bquote(italic(B)[t]*": Training set "*italic(D)[italic(T)]~"("*.(round(datCredit_train[,.N]/1000))*"k)"),
             "c_Valid"=bquote(italic(C)[t]*": Validation set "*italic(D)[italic(V)]~"("*.(round(datCredit_valid[,.N]/1000))*"k)"))

# - Create graph
(g4 <- ggplot(port.aggr_def, aes(x=timeVar, y=Prop)) + theme_minimal() + 
    labs(x=bquote("Default spell cohorts (ccyymm): entry time "*italic(t[e])), y=bquote("Resolution rate (%) of type "*~italic(kappa))) +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Sample, linetype=Sample)) + 
    geom_point(aes(colour=Sample, shape=Sample), size=1) + 
    # facets
    facet_wrap(Facet~., labeller = label_parsed, scales = "free", nrow=length(resolDef_levels), strip.position="right") + 
    #annotations
    geom_text(data=dat_anno_def, aes(x=x, y=y, hjust=hjust, vjust=vjust, label = Label), family=chosenFont, size=3, parse=T) + 
    # scale options
    scale_colour_manual(name=bquote("Sample "*italic(bar(D))), values=col.v, labels=label.v) + 
    scale_shape_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + scale_linetype_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# - Save graph
ggsave(g4, file=paste0(genFigPath_Res_anc, "ResolutionRates_Def_te_Subsample-", round(datCredit_smp[,.N]/1000),"k.png"), width=5000/(dpi*2.25), height=4000/(dpi*1.4), dpi=dpi, bg="white")

# - Cleanup
rm(dat_anno_def, chosenFont, col.v, label.v, colnames, datGraph_Def, port.aggr_def, port.aggr_def2, maxDate, minDate)




# ------ 8. Graphing performance spell resolution rates over time given resampled sets | Spell stop time (t_s)
# - Check representatives | dataset-level proportions should be similar
datCredit[timeVar==timeVar_Perf_Max, get(resolPerf_stop)] %>% table() %>% prop.table()
datCredit_train[timeVar==timeVar_Perf_Max, get(resolPerf_stop)] %>% table() %>% prop.table()
datCredit_valid[timeVar==timeVar_Perf_Max, get(resolPerf_stop)] %>% table() %>% prop.table()
# Checking the proportions of the subsampled dataset and creating the corresponding faceting labels
(Facet_Label_Perf <- datCredit_smp[timeVar==timeVar_Perf_Max, get(resolPerf_stop)] %>% table() %>% prop.table() %>% data.table()) # Saving these proportions as they are used in the facets
colnames(Facet_Label_Perf) <- c("PerfSpell_Resol_Stop", "Prior") # Renaming the columns
# - Subsetting the main long dataset to only include the necessary variables for performance spells
datGraph_Perf <- datGraph %>% subset (!is.na(ClusVar_Perf), select=c("ClusVar_Perf", "timeVar", "timeVar_Perf_Max", resolPerf_stop, "Sample"))
colnames(datGraph_Perf) <- c("ClusVar_Perf", "timeVar", "timeVar_Perf_Max", "PerfSpell_Resol_Stop", "Sample")


# - Setting some aggregation parameters, purely to facilitate graphing aesthetics
StartDte <- min(datCredit_smp$timeVar, na.rm=T)
EndDte <- max(datCredit_smp$timeVar, na.rm=T)
maxDate <- EndDte %m-% months(1)# A post-hoc filter, used for graphing purposes - left as the end of the sampling window
minDate <- StartDte # %m+% months(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window

# - Aggregate to monthly level and observe up to given point
port.aggr_perf <- merge(datGraph_Perf[timeVar==timeVar_Perf_Max, list(Sum_Total = .N), by=list(Sample,timeVar)],
                        datGraph_Perf[timeVar==timeVar_Perf_Max, list(Sum_Resol = .N), by=list(Sample,timeVar,PerfSpell_Resol_Stop)],
                        by=c("Sample", "timeVar"))[timeVar >= minDate & timeVar <= maxDate,]
port.aggr_perf[, Prop := Sum_Resol/Sum_Total]

# - Calculate MAE over time by sample, by performance event(s)
port.aggr_perf2 <- port.aggr_perf %>% pivot_wider(id_cols = c(timeVar), names_from = c(Sample, PerfSpell_Resol_Stop), values_from = Prop) %>% data.table()

# - Get the unique factors of the performance spell resolution type (used as facets)
resolPerf_levels <- unique(datCredit[!is.na(get(resolPerf_stop)), get(resolPerf_stop)])
# - Number of annotations "sets" to create
anno_n <- length(resolPerf_levels)
# - Initiating the annotation dataset
x_pos <- min(datCredit$timeVar) + round((max(datCredit$timeVar) - min(datCredit$timeVar))/2) # x-position fot the annotation
dat_anno_perf <- data.table(MAE = rep(0,anno_n*4),
                            Mean_EventRate = rep(0, anno_n*4),
                            stdError_EventRate = rep(0, anno_n*4),
                            margin_EventRate = rep(0, anno_n*4),
                            PerfSpell_Resol_Stop = unlist(lapply(resolPerf_levels, function(x){rep(x,4)})),
                            Dataset = rep(c("A-B","A-C","B-C", "B"), anno_n),
                            Label = rep(c(paste0("'MAE between '*italic(A[t])*' and '*italic(B[t])*'"),
                                          paste0("'MAE between '*italic(A[t])*' and '*italic(C[t])*'"),
                                          paste0("'MAE between '*italic(B[t])*' and '*italic(C[t])*'"), 
                                          paste0("'TTC-mean '*E(italic(B[t]))*'")), anno_n),
                            x = rep(x_pos,anno_n*4),
                            y = rep(Inf, anno_n*4), # c(c(0.48,0.45,0.42,0.38),c(0.65,0.62,0.59,0.55)), 
                            vjust = rep(c(1,2,3, 4), anno_n),
                            hjust=rep(0.5, anno_n*4))
# - Getting the column names to help compute the MAEs
colnames <- colnames(port.aggr_perf2)
# - Populating the annotation dataset
for (i in 1:anno_n){
  dat_anno_perf[i*4-3, MAE := mean(abs(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("a_Full_", resolPerf_levels[i])])[[1]] - subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and training dataset
  dat_anno_perf[i*4-2, MAE := mean(abs(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("a_Full_", resolPerf_levels[i])])[[1]] - subset(port.aggr_perf2, select=colnames[colnames %in% paste0("c_Valid_", resolPerf_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and validation dataset
  dat_anno_perf[i*4-1, MAE := mean(abs(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]] - subset(port.aggr_perf2, select=colnames[colnames %in% paste0("c_Valid_", resolPerf_levels[i])])[[1]]), na.rm = T)] # MAE between the training- and validation dataset
  dat_anno_perf[i*4,   mean_EventRate := mean(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]], na.rm=T)]
  dat_anno_perf[i*4,   stdError_EventRate := sd(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]], na.rm=T)/ sqrt(length(subset(port.aggr_perf2, select=colnames[colnames %in% paste0("b_Train_", resolPerf_levels[i])])[[1]]))]
}
# - Finessing the annotation dataset for plotting
ind <- seq(from=4, to=4*anno_n, by=4)
dat_anno_perf[ind, margin_EventRate := qnorm(1-(1-confLevel)/2) * stdError_EventRate]
dat_anno_perf[ind, Label := paste0(Label, " = ", sprintf("%.2f", mean_EventRate*100) , "% +-", sprintf("%.3f", margin_EventRate*100), "%'")]
dat_anno_perf[seq(from=1, to=4*anno_n)[!(seq(from=1, to=4*anno_n) %in% ind)], Label := paste0(Label, " = ", sprintf("%.4f",MAE*100), "%'")]
# - Adding an column to accommodate the facets
port.aggr_perf <- merge(port.aggr_perf, Facet_Label_Perf, by="PerfSpell_Resol_Stop")
port.aggr_perf[, Facet:=paste0('"', PerfSpell_Resol_Stop, ' (', sprintf("%.2f", Prior*100), '%)"')]
dat_anno_perf <- merge(dat_anno_perf, unique(subset(port.aggr_perf, select=c("PerfSpell_Resol_Stop", "Facet"))), by="PerfSpell_Resol_Stop")

# - Graphing parameters
chosenFont <- "Cambria"; dpi <- 340
col.v <- brewer.pal(9, "Set1")
label.v <- c("a_Full"=expression(italic(A)[t]*": Full set "*italic(D)),
             "b_Train"=bquote(italic(B)[t]*": Training set "*italic(D)[italic(T)]~"("*.(round(datCredit_train[,.N]/1000))*"k)"),
             "c_Valid"=bquote(italic(C)[t]*": Validation set "*italic(D)[italic(V)]~"("*.(round(datCredit_valid[,.N]/1000))*"k)"))

# - Create graph
(g5_1 <- ggplot(port.aggr_perf, aes(x=timeVar, y=Prop)) + theme_minimal() + 
      labs(x=bquote("Performing spell cohorts (ccyymm): stop time "*italic(t[s])), y=bquote("Resolution rate (%) of type "*italic(kappa))) +
      theme(text=element_text(family=chosenFont),legend.position = "bottom",
            axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
            strip.background=element_rect(fill="snow2", colour="snow2"),
            strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
      # main line graph with overlaid points
      geom_line(aes(colour=Sample, linetype=Sample)) + 
      geom_point(aes(colour=Sample, shape=Sample), size=1) + 
      # facets
      facet_wrap(Facet~., labeller = label_parsed, scales = "free", nrow=length(resolPerf_levels), strip.position="right") + 
      #annotations
      geom_text(data=dat_anno_perf, aes(x=x, y=y, hjust=hjust, vjust=vjust, label = Label), family=chosenFont, size=3, parse=T) + 
      # scale options
      scale_colour_manual(name=bquote("Sample "*italic(bar(D))), values=col.v, labels=label.v) + 
      scale_shape_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + scale_linetype_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + 
      scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
      scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# - Save graph
ggsave(g5_1, file=paste0(genFigPath_Res_anc, "ResolutionRates_Perf_ts_Subsample-", round(datCredit_smp[,.N]/1000),"k.png"), width=5000/(dpi*2.25), height=4000/(dpi*1.4), dpi=dpi, bg="white")

# - Create graph using only the first facet (conditional on the faceting variable having more than one level)
if (!is.na(resolPerf_stop2)){
  # Create graph
  (g5_2 <- ggplot(port.aggr_perf[PerfSpell_Resol_Stop==resolPerf_stop2,], aes(x=timeVar, y=Prop)) + theme_minimal() + 
     labs(x=bquote("Performing spell cohorts (ccyymm): stop time "*italic(t[s])), y=bquote("Resolution rate (%) of type "*italic(kappa))) +
     theme(text=element_text(family=chosenFont),legend.position = "bottom",
           axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
           strip.background=element_rect(fill="snow2", colour="snow2"),
           strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
     # main line graph with overlaid points
     geom_line(aes(colour=Sample, linetype=Sample)) + 
     geom_point(aes(colour=Sample, shape=Sample), size=1) + 
     # facets
     facet_wrap(Facet~., labeller = label_parsed, scales = "free", strip.position="right") + 
     #annotations
     geom_text(data=dat_anno_perf[PerfSpell_Resol_Stop==resolPerf_stop2, ], aes(x=x, y=y, hjust=hjust, vjust=vjust, label = Label), family=chosenFont, size=3, parse=T) + 
     # scale options
     scale_colour_manual(name=bquote("Sample "*italic(bar(D))), values=col.v, labels=label.v) + 
     scale_shape_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + scale_linetype_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + 
     scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
     scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))
  
  # Save graph
  dpi <- 170
  ggsave(g5_2, file=paste0(genFigPath_Res_anc, "ResolutionRates_Perf_ts_Subsample_Single_Facet-", round(datCredit_smp[,.N]/1000),"k.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
}

# - Cleanup
rm(dat_anno_perf, resolPerf_levels, ind, chosenFont, col.v, label.v, colnames, datGraph_Perf, port.aggr_perf, port.aggr_perf2, maxDate, minDate, Facet_Label_Perf)




# ------ 9. Graphing default spell resolution rates over time given resampled sets | Spell exit time (t_s)
# - Check representatives | dataset-level proportions should be similar
datCredit[timeVar==timeVar_Def_Max, get(resolDef_stop)] %>% table() %>% prop.table()
datCredit_train[timeVar==timeVar_Def_Max, get(resolDef_stop)] %>% table() %>% prop.table()
datCredit_valid[timeVar==timeVar_Def_Max, get(resolDef_stop)] %>% table() %>% prop.table()
# Checking the proportions of the subsampled dataset and creating the corresponding faceting labels
(Facet_Label_Def <- datCredit_smp[timeVar==timeVar_Def_Max, get(resolDef_stop)] %>% table() %>% prop.table() %>% data.table()) # Saving these proportions as they are used in the facets
colnames(Facet_Label_Def) <- c("DefSpell_Resol_Stop", "Prior") # Renaming the columns

# - Subsetting the main long dataset to only include the necessary variables for performance spells
datGraph_Def <- datGraph %>% subset(!is.na(ClusVar_Def), select=c("ClusVar_Def", "timeVar", "timeVar_Def_Max", resolDef_stop, "Sample"))
colnames(datGraph_Def) <- c("ClusVar_Def", "timeVar", "timeVar_Def_Max", "DefSpell_Resol_Stop", "Sample")

# - Setting some aggregation parameters, purely to facilitate graphing aesthetics
StartDte <- min(datCredit_smp$timeVar, na.rm=T)
EndDte <- max(datCredit_smp$timeVar, na.rm=T)
maxDate <- EndDte # A post-hoc filter, used for graphing purposes - left as the end of the sampling window
minDate <- StartDte # %m+% months(1) # A post-hoc filter, used for graphing purposes - set as one month after the sampling window

# - Aggregate to monthly level and observe up to given point
port.aggr_def <- merge(datGraph_Def[timeVar==timeVar_Def_Max, list(Sum_Total = .N), by=list(Sample,timeVar)],
                       datGraph_Def[timeVar==timeVar_Def_Max, list(Sum_Resol = .N), by=list(Sample,timeVar,DefSpell_Resol_Stop)],
                       by=c("Sample", "timeVar"))[timeVar >= minDate & timeVar <= maxDate,]
port.aggr_def[, Prop := Sum_Resol/Sum_Total]

# - Calculate MAE over time by sample, by competing event(s)
port.aggr_def2 <- port.aggr_def %>% pivot_wider(id_cols = c(timeVar), names_from = c(Sample, DefSpell_Resol_Stop), values_from = Prop) %>% as.data.table()

# - Get the unique factors of the default spell resolution type (used as facets)
resolDef_levels <- unique(datCredit[!is.na(get(resolDef_stop)), get(resolDef_stop)])
# - Number of annotations "sets" to create
anno_n <- length(resolDef_levels)

# - Initiating the annotation dataset
x_pos <- min(datCredit$timeVar) + round((max(datCredit$timeVar) - min(datCredit$timeVar))/2) + months(16)# x-position fot the annotation
dat_anno_def <- data.table(MAE = rep(0,anno_n*4),
                           Mean_EventRate = rep(0, anno_n*4),
                           stdError_EventRate = rep(0, anno_n*4),
                           margin_EventRate = rep(0, anno_n*4),
                           DefSpell_Resol_Stop = unlist(lapply(resolDef_levels, function(x){rep(x,4)})),
                           Dataset = rep(c("A-B","A-C","B-C", "B"), anno_n),
                           Label = rep(c(paste0("'MAE between '*italic(A[t])*' and '*italic(B[t])*'"),
                                         paste0("'MAE between '*italic(A[t])*' and '*italic(C[t])*'"),
                                         paste0("'MAE between '*italic(B[t])*' and '*italic(C[t])*'"), 
                                         paste0("'TTC-mean '*E(italic(B[t]))*'")), anno_n),
                           x = rep(x_pos,anno_n*4),
                           y =   rep(Inf, anno_n*4), # c(c(0.5,0.48,0.46,0.42),c(0.64,0.62,0.60,0.56)),
                           vjust = rep(c(1,2,3, 4), anno_n),
                           hjust=rep(0.5, anno_n*4))
# - Getting the column names to help compute the MAEs
colnames <- colnames(port.aggr_def2)
# - Populating the annotation dataset
for (i in 1:anno_n){
  dat_anno_def[i*4-3, MAE := mean(abs(subset(port.aggr_def2, select=colnames[colnames %in% paste0("a_Full_", resolDef_levels[i])])[[1]] - subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and training dataset
  dat_anno_def[i*4-2, MAE := mean(abs(subset(port.aggr_def2, select=colnames[colnames %in% paste0("a_Full_", resolDef_levels[i])])[[1]] - subset(port.aggr_def2, select=colnames[colnames %in% paste0("c_Valid_", resolDef_levels[i])])[[1]]), na.rm = T)] # MAE between the full- and validation dataset
  dat_anno_def[i*4-1, MAE := mean(abs(subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]] - subset(port.aggr_def2, select=colnames[colnames %in% paste0("c_Valid_", resolDef_levels[i])])[[1]]), na.rm = T)] # MAE between the training- and validation dataset
  dat_anno_def[i*4,   mean_EventRate := mean(subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]], na.rm=T)]
  dat_anno_def[i*4,   stdError_EventRate := sd(subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]], na.rm=T)/ sqrt(length(subset(port.aggr_def2, select=colnames[colnames %in% paste0("b_Train_", resolDef_levels[i])])[[1]]))]
}
# - Finessing the annotation dataset for plotting
ind <- seq(from=4, to=4*anno_n, by=4)
dat_anno_def[ind, margin_EventRate := qnorm(1-(1-confLevel)/2) * stdError_EventRate]
dat_anno_def[ind, Label := paste0(Label, " = ", sprintf("%.2f", mean_EventRate*100) , "% +-", sprintf("%.3f", margin_EventRate*100), "%'")]
dat_anno_def[seq(from=1, to=4*anno_n)[!(seq(from=1, to=4*anno_n) %in% ind)], Label := paste0(Label, " = ", sprintf("%.4f",MAE*100), "%'")]
# - Adding an column to accommodate the facets
port.aggr_def <- merge(port.aggr_def, Facet_Label_Def, by="DefSpell_Resol_Stop")
port.aggr_def[, Facet:=paste0('"', DefSpell_Resol_Stop, ' (', sprintf("%.2f", Prior*100), '%)"')]
dat_anno_def <- merge(dat_anno_def, unique(subset(port.aggr_def, select=c("DefSpell_Resol_Stop", "Facet"))), by="DefSpell_Resol_Stop")


# - Graphing parameters
chosenFont <- "Cambria"; dpi <- 340
col.v <- brewer.pal(9, "Set1")
label.v <- c("a_Full"=expression(italic(A)[t]*": Full set "*italic(D)),
             "b_Train"=bquote(italic(B)[t]*": Training set "*italic(D)[italic(T)]~"("*.(round(datCredit_train[,.N]/1000))*"k)"),
             "c_Valid"=bquote(italic(C)[t]*": Validation set "*italic(D)[italic(V)]~"("*.(round(datCredit_valid[,.N]/1000))*"k)"))

# - Create graph
(g6 <- ggplot(port.aggr_def, aes(x=timeVar, y=Prop)) + theme_minimal() + 
    labs(x=bquote("Default spell cohorts (ccyymm): stop time "*italic(t[s])), y=bquote("Resolution rate (%) of type "*italic(kappa))) +
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # main line graph with overlaid points
    geom_line(aes(colour=Sample, linetype=Sample)) + 
    geom_point(aes(colour=Sample, shape=Sample), size=1) + 
    # facets
    facet_wrap(Facet~., labeller = label_parsed, scales = "free", nrow=length(resolDef_levels), strip.position="right") + 
    #annotations
    geom_text(data=dat_anno_def, aes(x=x, y=y, hjust=hjust, vjust=vjust, label = Label), family=chosenFont, size=3, parse=T) + 
    # scale options
    scale_colour_manual(name=bquote("Sample "*italic(bar(D))), values=col.v, labels=label.v) + 
    scale_shape_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + scale_linetype_discrete(name=bquote("Sample "*italic(bar(D))), labels=label.v) + 
    scale_y_continuous(breaks=pretty_breaks(), label=percent) + 
    scale_x_date(date_breaks=paste0(6, " month"), date_labels = "%b %Y"))

# - Save graph
ggsave(g6, file=paste0(genFigPath_Res_anc, "ResolutionRates_Def_ts_Subsample-", round(datCredit_smp[,.N]/1000),"k.png"), width=5000/(dpi*2.25), height=4000/(dpi*1.4), dpi=dpi, bg="white")

# - Cleanup
rm(dat_anno_def, ind, chosenFont, col.v, label.v, colnames, datGraph_Def, port.aggr_def, port.aggr_def2, maxDate, minDate)






# --- Cleanup
suppressWarnings(rm(g1, g2, g3, g4, g5_1, g5_2, g6, datCredit, datCredit_smp, stratifiers, clusVar, Counter, timeVar, End_Dte, Start_Dte, datGraph))
