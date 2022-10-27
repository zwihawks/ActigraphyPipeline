# -----------------------------
# title: "Actigraphy pre- & post-processing demo"
# author: "Zoë Hawks"
# updated: "04/06/2022"
# Build info:
  # R version: 4.1.1 (2021-08-10) "Kick Things"
  # R Studio version: 1.4.1717
# -----------------------------

# -----------------------------
# Load libraries & define function inputs
# -----------------------------

path <- "" # path to directory containing current script & ExternalFuncCode folder
pat_ID <- 6 # participant ID
TMB_file <- T # do we have EMA data from TMB? (T for participants, F for some internal testing)

# Source functions
source(paste0(path, "ExternalFuncCode/Actigraphy_preFunc.R"))
source(paste0(path, "ExternalFuncCode/Actigraphy_postFunc.R"))
source(paste0(path, "ExternalFuncCode/wrangle_sleep.R"))
source(paste0(path, "ExternalFuncCode/parse_activity.R"))


# -----------------------------------
# EMA data
# -----------------------------------

if (TMB_file) {
  EMA_data <- read.csv(list.files(path = paste0(path, "Raw/ID_", as.character(pat_ID)), 
                                  pattern = "studies.*.csv$", full.names = TRUE))
}

# -----------------------------
# Main loop
# -----------------------------
# For each participant with TMB_file, run both with and without sleep log
for (sleepLog in c(TRUE, FALSE)) { 
  
  sleepLog <- as.logical(TMB_file*sleepLog) # sets sleepLog to FALSE if no TMB data
  print(paste0("Running with sleepLog = ", sleepLog))
  
  # -----------------------------------
  # Part 1: Wrangle sleep data
  # -----------------------------------
  # Create sleep log from EMA data, to be called in Actigraphy_preFunc.R
  # Debating whether to set replaceNA to TRUE vs. FALSE (or, to replace iff first night is NA)
  # With NAs, will only retain raw data for days (WW) with guiders (cf. save_ms5raw_without_invalid = TRUE in preFunc)
  
  if (sleepLog) {
    n_nights <- wrangle_sleep(path = path, patID = pat_ID, EMA_data = EMA_data, replaceNA = FALSE)
  } else {
    n_nights <- c()
  }
  
  # -----------------------------------
  # Part 2: Pre-processing
  # -----------------------------------
  # perform actigraphy pre-processing on raw bin file + EMA sleep log (iff log = TRUE) from single participant
  # Still evaluating whether parameters passed to myscript.R effectively estimate vigorous steps
  
  actigraphy_preP(path = path, 
                  patID = pat_ID, 
                  n_nights = n_nights,
                  log = sleepLog)
  
  # -----------------------------------
  # Part 3: Post-processing
  # -----------------------------------
  # outputs day summary with sleep and activity metrics, minute-level activity file, & optional QC plot
  actigraphy_postP(path = path, patID = pat_ID, output_plot = TRUE, log = sleepLog)
  
  # -----------------------------------
  # Part 4: EMA-level activity data
  # -----------------------------------
  # create log with step and/or activity data linked to EMA sitting ID
  if (TMB_file) {
    parse_activity(path = path, patID = pat_ID, EMA_data = EMA_data, 
                   direction = -1, offset_min = c(5, 20, 60, 180), 
                   steps = TRUE, activity = TRUE, log = sleepLog,
                   export = TRUE, return = FALSE)
  }
  
  if (!TMB_file) {
    print(paste0("No TMB file, breaking after one loop"))
    break
  }
}


