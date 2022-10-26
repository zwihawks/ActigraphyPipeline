# -----------------------------
# Actigraphy pre-processing
# -----------------------------

actigraphy_preP <- function (path, patID, n_nights = c(), log = TRUE) {
  #' perform actigraphy pre-processing on raw bin file + EMA sleep log from single participant
  #' 
  #' @details 
  #' Function reads in raw data from path + Raw/ID_patID folder
  #' Outputs contents to path + Processed/Pre folder
  #' Expects sleep log data & optimized for left wrist data from middle-aged adults
  #' Will need to modify part 3-4 if no sleep log
  #' Will need to modify part 5 for use in other populations
  #' Refer to: https://cran.r-project.org/web/packages/GGIR/vignettes/GGIR.html
  #' 
  #' @param path string indicating file pathway to folder containing Processed data folder
  #' @param patID participant ID, strictly numeric
  #' @param n_nights number of nights from start of actigraphy collection through final EMA
  #' @param log Boolean indicating whether to use EMA sleep/wake self-report data
  #' 
  if (!require("remotes")) install.packages("remotes"); require ("remotes")
  if (!require("GGIR")) remotes::install_github("wadpac/GGIR"); require ("GGIR")
  source(paste0(path, "ExternalFuncCode/myscript.R"))
  
  if (log) {
    outputdir <- paste0(path, "Processed/Pre/Log")
    log_loc <- c(paste0(path, "Raw/ID_", as.character(patID), "/sleepLog", as.character(patID), ".csv"))
    colid <- 1
    coln1 <- 2 
    do_visual <- TRUE
    outliers_only <- TRUE 
    criterror <- .5
  } else {
    outputdir <- paste0(path, "Processed/Pre/noLog")
    log_loc <- c()
    colid <- c()
    coln1 <- c()
    do_visual <- FALSE
    outliers_only <- FALSE 
    criterror <- c()
  }
  
  g.shell.GGIR(
    # ----------------
    # General parameters (parts 1-2)
    # ----------------
    datadir=paste0(path, "Raw/ID_", as.character(patID)),
    # e.g., specify data/processed/actigraphy >> creates "output_[input_subfolder]" in actigraphy folder
    outputdir=outputdir,
    overwrite = TRUE,
    do.cal = TRUE, # whether to apply auto-calibration to account for gravity
    idloc = 2, # where to find patID; 1 = file header, 2 = filename
    desiredtz = "America/New_York",
    do.enmoa = TRUE, # for activity estimation, below
    do.enmo = FALSE, # for activity estimation, below
    do.parallel = FALSE, # set to FALSE for step count algorithm
    strategy = 1,
    hrs.del.start = 0,
    hrs.del.end = 0,
    includedaycrit = 16, 
    # ----------------
    # Relevant to sleep estimation (parts 3-4)
    # ----------------
    ignorenonwear = TRUE, # ignore detected monitor non-wear periods in the detection of sustained inactivity bouts
    loglocation = log_loc, 
    do.visual = do_visual,
    outliers.only = outliers_only, # Relevant for do.visual == TRUE; if TRUE, used in conjunction with criterror (see below)
    criterror = criterror, 
    excludefirstlast = FALSE,
    colid = colid, 
    coln1 = coln1, 
    nnights = n_nights,
    sleepwindowType = "SPT", # Can also use "TimeInBed," which outputs sleep efficiency metrics
    def.noc.sleep = c(1), # how sleep period window should be estimated in absence of sleep guider
    # ----------------
    # Merging physical activity with sleep estimation (part 5)
    # ----------------
    # numbers below based on Eslinger et al. (2011) - n = 60, 40-65 y.o.
    # other activity cut-points for other body parts & ages (older/younger adults)
    # decision to have all participants wear on L wrist & use same cut-points
    # LEFT WRIST, ADULTS
    acc.metric = "ENMOa",
    threshold.lig = (217/(80*60)) * 1000,
    threshold.mod = (645/(80*60)) * 1000,
    threshold.vig = (1810/(80*60)) * 1000,
    mvpathreshold = (645/(80*60)) * 1000,
    save_ms5rawlevels = TRUE,
    part5_agg2_60seconds = TRUE,
    save_ms5raw_without_invalid = TRUE,
    # ----------------
    do.report = c(2, 4, 5),
    timewindow = c("WW"),
    visualreport = TRUE,
    dofirstpage = TRUE,
    myfun=vcs_fun # implement step count algorithm from external function
  )
}

