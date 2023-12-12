# =============================== PROPORTIONAL HAZARD ASSUMPTION TEST ===============================
# Assessing the proportional hazards (PH) assumption using both a stratified approach (KM-estimators)
# and Schoenfeld residual plots.
# ---------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Marcel Muller
# --------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 4a.Data_Partition
# --------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Loading in required dataset
if (!exists('datCredit_real')) {  unpack.ffdf(paste0(genPath,"creditdata_analytics"), tempPath)}

# - Creating a dataset containing a single observation for each performance spell (used to create a survival object)
dat_Surv <- subset(datCredit_real, PerfSpell_Counter==1)

# Grapphing parameters
chosenFont <- "Cambria"
dpi <- 235




# ------ 2. Survival Function Comparison for Loan Type
# --- Survival objects
km_Default_lm_tpe1 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                              id = PerfSpell_Key, data=dat_Surv[LN_TPE=="WHL",])

km_Default_lm_tpe2 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                              id = PerfSpell_Key, data=dat_Surv[LN_TPE=="CHL",])

# --- Plotting dataset
dat_Plot <- rbind(cbind(surv_summary(km_Default_lm_tpe1), Loan_Type = "WHL"),
                  cbind(surv_summary(km_Default_lm_tpe2), Loan_Type = "CHL"))

# --- Plot
col.v <- brewer.pal(9, "Greens")[c(4,9)]
label.v <- c("WHL", "CHL")
chosenFont <- "Cambria"

(g_surv_comp_ln_tpe <- ggplot(dat_Plot, aes(x=log(time), y=log(-log(surv)), factor=Loan_Type)) +
                       theme_minimal() +
                       geom_line(aes(colour=Loan_Type, linetype=Loan_Type)) +
                       scale_colour_manual(name="Loan Type", values=col.v, labels=label.v) +
                       scale_linetype(name="Loan Type", labels=label.v) +
                       scale_x_continuous(name = bquote(italic(log(t))), limits = c(0,log(200))) +
                       scale_y_continuous(name = bquote(italic(log(-log(hat(S)(t)))))) +
                       theme(legend.position = "bottom", text = element_text(family = chosenFont)))

ggsave(file = "PH_Assumption_LN_TPE.png", plot = g_surv_comp_ln_tpe, device = 'png', path = genFigPath_ass)




# ------ 3. Survival Function Comparison for Payment Method
# --- Survival objects
km_Default_pmnt_meth1 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                              id = PerfSpell_Key, data=dat_Surv[slc_pmnt_method=="MISSING_DATA",])

km_Default_pmnt_meth2 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                              id = PerfSpell_Key, data=dat_Surv[slc_pmnt_method=="Salary",])

km_Default_pmnt_meth3 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                                 id = PerfSpell_Key, data=dat_Surv[slc_pmnt_method=="Statement",])

km_Default_pmnt_meth4 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                                 id = PerfSpell_Key, data=dat_Surv[slc_pmnt_method=="Debit Order FNB account",])

km_Default_pmnt_meth5 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                                 id = PerfSpell_Key, data=dat_Surv[slc_pmnt_method=="Debit Order other bank",])

km_Default_pmnt_meth6 <- survfit(Surv(time=TimeInPerfSpell-1, time2=PerfSpell_Age, PerfSpellResol_Type_Hist=="Defaulted", type="counting") ~ 1,
                                 id = PerfSpell_Key, data=dat_Surv[slc_pmnt_method=="Suspense",])


# --- Plotting dataset
dat_Plot2 <- rbind(cbind(surv_summary(km_Default_pmnt_meth1), PMNT_Method = "Missing"),
                   cbind(surv_summary(km_Default_pmnt_meth2), PMNT_Method = "Salary"),
                   cbind(surv_summary(km_Default_pmnt_meth3), PMNT_Method = "Statement"),
                   cbind(surv_summary(km_Default_pmnt_meth4), PMNT_Method = "Debit Order (FNB)"),
                   cbind(surv_summary(km_Default_pmnt_meth5), PMNT_Method = "Debit Order (Other Bank)"),
                   cbind(surv_summary(km_Default_pmnt_meth6), PMNT_Method = "Supsense"))

# --- Plot
col.v2 <- brewer.pal(9, "Greens")[c(3,4,5,7,8,9)]
label.v2 <- c("Missing", "Salary", "Statement", "Debit Order (FNB)", "Debit Order (Other Bank)", "Supsense")
chosenFont <- "Cambria"

(g_surv_comp_slc_pmnt_method <- ggplot(dat_Plot2, aes(x=log(time), y=log(-log(surv)), factor=PMNT_Method)) +
                                theme_minimal() +
                                geom_line(aes(colour=PMNT_Method, linetype=PMNT_Method)) +
                                scale_colour_manual(name="Payment Method", values=col.v2, labels=label.v2) +
                                scale_linetype(name="Payment Method", labels=label.v2) +
                                scale_x_continuous(name = bquote(log~italic((t))), limits = c(0,log(200))) +
                                scale_y_continuous(name = bquote(log(-log(italic(hat(S)(t)))))) +
                                theme(legend.position = "bottom", text = element_text(family = chosenFont)))

ggsave(file = "PH_Assumption_pmnt_method.png", plot = g_surv_comp_slc_pmnt_method, device = 'png', path = genFigPath_ass, dpi = 500)

# --- House keeping
suppressWarnings(rm(datCredit_real, dat_Surv, km_Default_lm_tpe1, km_Default_lm_tpe2, km_Default_pmnt_meth1, km_Default_pmnt_meth2, km_Default_pmnt_meth3,
                    km_Default_pmnt_meth4, km_Default_pmnt_meth5, km_Default_pmnt_meth5, dat_Plot, dat_Plot2, col.v, col.v2, label.v))




# ------ 4. Schoenfeld Residual Plots
# --- Remove some unessary variables to improve performance
colnames_datCredit_real <- colnames(datCredit_real)
datCredit_real <- datCredit_real %>% subset(select = -c(which(colnames_datCredit_real == "M_Inflation_Growth_1"):which(colnames_datCredit_real == "M_Repo_Rate_Vol_12")))


# --- Subsetting to only include observations with TimeInPerfSpell <= 300 and then grouping the data according to PerfSpell_Key and TimeInPerfSpell
datCredit_real <- datCredit_real[TimeInPerfSpell <= 300, ] %>% group_by(PerfSpell_Key, TimeInPerfSpell)


# --- Fitting a Cox PH model with multiple covariates including slc_pmnt_method
cph_Default_PH_test <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1) ~ 
                             Principal_wins + InterestRate_Margin +
                             LN_TPE + slc_pmnt_method,
                             data=datCredit_real, id=PerfSpell_Key)
# Principal_wins      : coef = -0.000000108801; exp(coef) = 0.999999891199; se(coef) = 0.000000006939; z = -15.679; p < 0.0000000000000002
# InterestRate_Margin : coef = 10.197476284966; exp(coef) = 26835.375702168672; se(coef) = 0.222489467262; z = 45.834; p < 0.0000000000000002
# LN_TPEWHL           : coef = -0.059004270950; exp(coef) = 0.942702742912; se(coef) = 0.013649668519; z = -4.323; p = 0.0000154
# slc_pmnt_methodDebit Order other bank: coef = 0.486377487299; exp(coef) = 1.626413830911; se(coef) = 0.012002908070; z = 40.522; p < 0.0000000000000002
# slc_pmnt_methodMISSING_DATA          : coef = 1.084559390140; exp(coef) = 2.958136147601; se(coef) = 0.010216908686; z = 106.153; p < 0.0000000000000002
# slc_pmnt_methodSalary                : coef = 1.023280319783; exp(coef) = 2.782306666418; se(coef) = 0.015606346742; z = 65.568; p < 0.0000000000000002
# slc_pmnt_methodStatement             : coef = 2.424936556182; exp(coef) = 11.301512385438; se(coef) = 0.009159261660; z = 264.752; p < 0.0000000000000002
# slc_pmnt_methodSuspense              : coef = 5.206898573136; exp(coef) = 182.527085030424; se(coef) = 0.011133146969 ; z = 467.693; p < 0.0000000000000002


# --- Testing the PH assumption using the PH test
# - Computing the (unscaled) Schoenfeld residuals using the custom functino
ass_cph_Default_PH_test <- cph_schoen(cph = cph_Default_PH_test, dat_train = datCredit_real, id = "PerfSpell_Key", time = "TimeInPerfSpell", status = "DefaultStatus1", verbose = F, max_time = 240)

ass_cph_Default_PH_test_data <- ass_cph_Default_PH_test$data
max_time <- 240

# Graphing parameters
chosenFont <- "Cambria"
dpi <- 360

# - [Principal]
g_Sch_Res_Principal <- ggplot(ass_cph_Default_PH_test_data[, list(Time, Sch_Res_Principal)], aes(x=Time, y=Sch_Res_Principal)) +
                       theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                       theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                       geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                       scale_x_continuous(limits = c(NA, max_time)) +
                       scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_Principal.png", g_Sch_Res_Principal, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# - [InterestRate_Margin]
g_Sch_Res_InterestRate_Margin <- ggplot(ass_cph_Default_PH_test_data[, list(Time, Sch_Res_InterestRate_Margin)], aes(x=Time, y=Sch_Res_InterestRate_Margin)) +
                                        theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                                        theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                                        geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                                        scale_x_continuous(limits = c(NA, max_time)) +
                                        scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_InterestRate_Margin.png", g_Sch_Res_InterestRate_Margin, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# - [LN_TPE]
# `WHL`
col.v <- brewer.pal(11, "Paired")[c(2,4,7,8)]
shape.v <- c(2,4)
label.v <- c("CHL" = "CHL",
             "WHL" = "WHL")

g_Sch_Res_LN_TPE <- ggplot(ass_cph_Default_PH_test_data, aes(x=Time, y=Sch_Res_LN_TPE_WHL)) +
                    theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                    theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                    geom_point(aes(col=Var_Val_LN_TPE_WHL, shape=Var_Val_LN_TPE_WHL)) + geom_smooth(aes(linetype="solid"), col=col.v[4], fill=col.v[3]) +
                    scale_x_continuous(limits = c(NA, max_time)) +
                    scale_color_manual(name = "Loan Type", values = col.v[1:2], labels = label.v) +
                    scale_shape_manual(name = "Loan Type", values = shape.v, labels = label.v) +
                    scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL) +
                    scale_fill_manual(name = "Loeass-smoother", values = col.v[3], labels = NULL)
                    
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_LN_TPE.png", g_Sch_Res_LN_TPE, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# - [slc_pmnt_method]
# 'Debit order other bank'
g_Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank <- ggplot(ass_cph_Default_PH_test_data[, list(Time, `Sch_Res_slc_pmnt_method_Debit Order other bank`)], aes(x=Time, y=`Sch_Res_slc_pmnt_method_Debit Order other bank`)) +
                                                    theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                                                    theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                                                    geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                                                    scale_x_continuous(limits = c(NA, max_time)) +
                                                    scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank.png", g_Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# 'MISSING_DATA'
g_Sch_Res_slc_pmnt_method_MISSING_DATA <- ggplot(ass_cph_Default_PH_test_data[, list(Time, Sch_Res_slc_pmnt_method_MISSING_DATA)], aes(x=Time, y=Sch_Res_slc_pmnt_method_MISSING_DATA)) +
                                          theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                                          theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                                          geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                                          scale_x_continuous(limits = c(NA, max_time)) +
                                          scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_slc_pmnt_method_MISSING_DATA.png", g_Sch_Res_slc_pmnt_method_MISSING_DATA, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# 'Salary'
g_Sch_Res_slc_pmnt_method_Salary <- ggplot(ass_cph_Default_PH_test_data[, list(Time, Sch_Res_slc_pmnt_method_Salary)], aes(x=Time, y=Sch_Res_slc_pmnt_method_Salary)) +
                                    theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                                    theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                                    geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                                    scale_x_continuous(limits = c(NA, max_time)) +
                                    scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_slc_pmnt_method_Salary.png", g_Sch_Res_slc_pmnt_method_Salary, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# 'Statement'
g_Sch_Res_slc_pmnt_method_Statement <- ggplot(ass_cph_Default_PH_test_data[, list(Time, Sch_Res_slc_pmnt_method_Statement)], aes(x=Time, y=Sch_Res_slc_pmnt_method_Statement)) +
                                       theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                                       theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                                       geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                                       scale_x_continuous(limits = c(NA, max_time)) +
                                       scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_slc_pmnt_method_Statement.png", g_Sch_Res_slc_pmnt_method_Statement, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# 'Suspense'
g_Sch_Res_slc_pmnt_method_Suspense <- ggplot(ass_cph_Default_PH_test_data[, list(Time, Sch_Res_slc_pmnt_method_Suspense)], aes(x=Time, y=Sch_Res_slc_pmnt_method_Suspense)) +
                                      theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                                      theme(text=element_text(family=chosenFont),legend.position = "bottom") +
                                      geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
                                      scale_x_continuous(limits = c(NA, max_time)) +
                                      scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
# Saving the plot
ggsave(file = "PH_Assumption_Sch_Res_slc_pmnt_method_Suspense.png", g_Sch_Res_slc_pmnt_method_Suspense, device = 'png', path = genFigPath_ass, height = 1800/dpi, width = 2200/dpi, dpi=dpi)

# - All levels of [slc_pmnt_method]
dat_slc_pmnt_method <- rbind(ass_cph_Default_PH_test_data[,list(ID = ID, Time = Time, Sch_Res = `Sch_Res_slc_pmnt_method_Debit Order other bank`, Level = "DebitOrderOtherBank")],
                             ass_cph_Default_PH_test_data[,list(ID = ID, Time = Time, Sch_Res = Sch_Res_slc_pmnt_method_MISSING_DATA, Level = "MISSING_DATA")],
                             ass_cph_Default_PH_test_data[,list(ID = ID, Time = Time, Sch_Res = Sch_Res_slc_pmnt_method_Salary, Level = "Salary")],
                             ass_cph_Default_PH_test_data[,list(ID = ID, Time = Time, Sch_Res = Sch_Res_slc_pmnt_method_Statement, Level = "Statement")],
                             ass_cph_Default_PH_test_data[,list(ID = ID, Time = Time, Sch_Res = Sch_Res_slc_pmnt_method_Suspense, Level = "Suspense")])
# Plot
col.v <- brewer.pal(11, "Spectral")[c(1,3,7,9,11)]
shape.v <- c(1,3,7,9,11)
linetype.v <- c(1,3,7,9,11)
label.v <- c("DebitOrderOtherBank" = "Debit Order Other Bank",
             "MISSING_DATA" = "Missing",
             "Salary" = "Salary",
             "Statement" = "Statement",
             "Suspense" = "Suspense")

g_Sch_Res_slc_pmnt_method <- ggplot(dat_slc_pmnt_method, aes(x=Time, y=Sch_Res)) +
                             theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
                             theme(text=element_text(family=chosenFont),legend.position = "bottom") + 
                             geom_point(aes(shape=Level, col=Level)) + geom_smooth(aes(col=Level, fill=Level, linetype=Level)) +
                             scale_x_continuous(limits = c(NA, max_time)) +
                             scale_color_manual(name = "SLC Payment Method", values = col.v, labels = label.v) +
                             scale_fill_manual(name = "SLC Payment Method", values = col.v, labels = label.v) +
                             scale_shape_manual(name = "SLC Payment Method", values = shape.v, labels = label.v) +
                             scale_linetype_manual(name = "SLC Payment Method", values = linetype.v, labels = label.v)

ggsave(file = "PH_Assumption_Sch_Res_slc_pmnt_method.png", g_Sch_Res_slc_pmnt_method, device = 'png', path = genFigPath_ass, height = 2000/dpi, width = 3000/dpi, dpi=dpi)


