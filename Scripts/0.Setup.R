# =================================== SETUP =============================================
# Setting up R environment, parameters, and function definitions
# ---------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha, Roelinde Bester, Marcel Muller, Roland Breedt

# DESCRIPTION: 
# This script installs and loads various libraries and packages, compiles all
# custom functions, and set requisite parameters.
# ---------------------------------------------------------------------------------------
# -- Inputs:
#   - DelinqM.R | Delinquency measures and related functions
# =======================================================================================



# ================ 0. Library setup

# ------ Install and load packages
# - data access and big data management
require(haven) # for SAS imports
require(ETLUtils)
require(ffbase)
require(ff)
tempPath <- "C:/TempData"; options("fftempdir"=tempPath)

# for data wrangling
require(tidyr)
require(dplyr)
require(data.table)
require(lubridate)
require(readr)
require(bit64) # for very big numeric values
require(stringr) # common string operations, e.g, str_pad
require(purrr) # mapping functions from tidyverse in working with matrices, lists

# for advanced looping functionality in simulation tasks
require(doBy)
require(foreach)
require(doParallel)

# for analyses
require(Hmisc)
require(moments) # for using skewness() function
require(regclass) # for VIF

# for modelling
require(survival) # for survival modelling
require(zoo)
require(car)
require(survivalROC) # for time-dependent ROC-analysis from Heagerty et al.
#require(survAUC) # for time-dependent ROC-analysis (alternative from Potapov et al.)
#require(tdROC) # for time-dependent ROC-analysis ([outdated?] alternative from Li et al.)
#require(timeROC) # for time-dependent ROC-analysis ([outdated?] alternative from Blanche)
require(pROC); require(ROCR) # both for cross-sectional ROC-analysis (main:pROC)
require(discSurv)
require(MASS)

#for plots
require(ggplot2)
require(ggpp) # Extensions to ggplot2, particularly geom_table
require(scales)
require(ggthemes)
require(RColorBrewer)
require(extrafont) #remotes::install_version("Rttf2pt1", version = "1.3.8"); Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.55.0/bin/gswin32c.exe"); font_import(); loadfonts(); loadfonts(device="win")
require(survminer)
require(gridExtra)
require(corrplot)
require(Metrics)



# ================ 1. Parametrisation

# - general R options
options(scipen=999) # Suppress showing scientific notation

# - Parameters used in calculating delinquency measures
sc.Thres <- 0.9; # repayment ratio - g1
d <- 3 # default threshold for g0/g1-measures of delinquency (payments in arrears)
k <- 6 # Probation period


# --- User specific paths
# - Getting the username of the current active user
Username <- Sys.getenv("USERNAME")

# - Setting the paths conditional on the current active user
if (Username == "WRQ"){ ### - Paths for AB
  # Custom path where R-scripts are saved
  path_cust <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Scripts/"
  
  # Common path for storing important R-objects as back-up
  genObjPath <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Objects/"
  genObjPath_Res <- "C:/Users/WRQ/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Objects/Resampling Schemes/" # Excel sheets for resampling tests
  
  # Common path for saving important analytics (e.g., sampling)
  genFigPath <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/" # General folder path
  genFigPath_CumHaz <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Baseline Hazard Comparison/" # Figures of baseline (cumulative) hazard comparison
  genFigPath_Res <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Main/" # Figures of resampling tests
  genFigPath_Res_anc <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Exploratory/" # Ancillary figures of resampling tests
  genFigPath_ass <- "C:/Users/WRQ/OneDrive - FRG/Analytix/Research/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Modelling Assumptions/" # Modelling assumptions path
} else if (Username == "R5532132"){ ### - Paths for MM
  # Custom path where R-scripts are saved
  path_cust <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Scripts/"
  
  # Common path for storing important R-objects as back-up
  genObjPath <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Objects/" # General objects"
  genObjPath_Res <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Objects/Resampling Schemes/" # Excel sheets for resampling tests
  
  # Common path for saving important analytics (e.g., sampling)
  genFigPath <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/" # General folder path
  genFigPath_CumHaz <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Baseline Hazard Comparison/" # Figures of baseline (cumulative) hazard comparison
  genFigPath_Res <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Main/" # Figures of resampling tests
  genFigPath_Res_anc <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Exploratory/" # Ancillary figures of resampling tests
  genFigPath_ass <- "C:/Users/R5532132/OneDrive - FRG/GCRM/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Modelling Assumptions/" # Modelling assumptions path
  path_cust <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Scripts/"
  
  # Common path for storing important R-objects as back-up
  genObjPath <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Objects/" # General objects"
  genObjPath_Res <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Objects/Resampling Schemes/" # Excel sheets for resampling tests
  
  # Common path for saving important analytics (e.g., sampling)
  genFigPath <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/" # General folder path
  genFigPath_CumHaz <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Baseline Hazard Comparison/" # Figures of baseline (cumulative) hazard comparison
  genFigPath_Res <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Main/" # Figures of resampling tests
  genFigPath_Res_anc <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Exploratory/" # Ancillary figures of resampling tests
  genFigPath_ass <- "C:/Users/R5532132/OneDrive - FRG/Default Survival Modelling/Term-Structure-Modelling-RetailMortgages/Figures/Modelling Assumptions/" # Modelling assumptions path
} else if (Username == "R5668395"){ ### - Paths for RB
  # Custom path where R-scripts are saved
  path_cust <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Scripts/"
  
  # Common path for storing important R-objects as back-up
  genObjPath <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Objects/" # General objects"
  genObjPath_Res <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Objects/Resampling Schemes/" # Excel sheets for resampling tests
  
  # Common path for saving important analytics (e.g., sampling)
  genFigPath <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Figures/" # General folder path
  genFigPath_CumHaz <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Figures/Baseline Hazard Comparison/" # Figures of baseline (cumulative) hazard comparison
  genFigPath_Res <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Main/" # Figures of resampling tests
  genFigPath_Res_anc <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Figures/Resampling Comparison/Exploratory/" # Ancillary figures of resampling tests
  genFigPath_ass <- "C:/Users/R5668395/Documents/GitHub TermStruc/Term-Structure-Modelling-RetailMortgages/Figures/Modelling Assumptions/" # Modelling assumptions path
} else { ### - If username not found in if statements, execution stops
  stop("User-specific paths not set for current user: ", Username, ". Please fix in setup script (0.Setup.R) before continuing.")
}


# --- General paths
# - Common path for saving big data objects
genPath <- "C:/Data/DefaultSurv_Data/"

# - Common path for importing raw data
genRawPath <- "C:/Data/"



# ================ 2. Custom functions

# ------ Custom function definitions
# - Load all custom functions defined in a separate R-script
source(paste0(path_cust,"0a(i).CustomFunctions.R"))

# - True End procedure functions defined in a separate R-script
source(paste0(path_cust,"0a(ii).TruEnd.R"))

# - Compile Delinquency Calculation Functions (CD, MD/DoD)
source(paste0(path_cust,'0a(iii).DelinqM.R'))

# - Survival functions
source(paste0(path_cust,'0a(iv).SurvFunc.R'))

