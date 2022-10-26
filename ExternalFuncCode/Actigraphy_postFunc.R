# -----------------------------
# Post-processing
# -----------------------------

actigraphy_postP <- function (path, patID, output_plot = FALSE, log = FALSE) {
  #' perform actigraphy post-processing
  #' 
  #' @details 
  #' Function reads in data from path + Processed/Pre/ID_patID folders
  #' Outputs contents to path + Processed/Post folder
  #' Outputs include: 
  #' 1. day summary with sleep and activity metrics; 
  #' 2. raw file with minute-level activity and step counts
  #' 3. optional QC plot. If plot looks wonky, adjust parameters in myscript.R, spec. indices 2 & 4 
  #' Code book: refer to https://cran.r-project.org/web/packages/GGIR/vignettes/GGIR.html
  #' 
  #' @param path string indicating pathway to folder that contains the Processed data folder
  #' @param patID participant ID, strictly numeric
  #' @param output_plot boolean indicating whether to save QC plot as TIFF
  #' @param log T/F indicating whether to use sleep data estimated using sleep log
  #' 
  if (!require("anytime")) install.packages("anytime"); require ("anytime")
  if (!require("remotes")) install.packages("remotes"); require ("remotes")
  if (!require("tidyverse")) remotes::install_github("tidyverse/tidyverse@v1.3.1"); require ("tidyverse")
  if (!require("lubridate")) remotes::install_github("tidyverse/lubridate@v1.8.0"); require ("lubridate")
  
  # Update home base
  path <- paste0(path, "Processed")
  if (log) {
    appendPre <- "/Pre/Log/output_ID_"
    appendPost <- "/Post/Log/output_ID_"
  } else {
    appendPre <- "/Pre/noLog/output_ID_"
    appendPost <- "/Post/noLog/output_ID_"
  }
  
  # read in and merge pre-processing files
  load(list.files(path = paste0(path, appendPre, as.character(pat_ID), "/meta/ms2.out"),
                  pattern = "*.RData", full.names = TRUE))
  files_to_read <- 
    list.files(path = paste0(path, appendPre, as.character(pat_ID), "/meta/ms5.outraw"), 
               pattern = "*.csv", recursive = TRUE, full.names = TRUE) 
  
  # determine tz
  tz1 <- ifelse(str_detect(IMP$metashort$timestamp[1], "-0500"), "Etc/GMT+5", "Etc/GMT+4")
  for (i in 1:length(IMP$metashort$timestamp)) {
    if(tz1 == "Etc/GMT+5" & !str_detect(IMP$metashort$timestamp[i], "-0500")) {
      stop(paste0("ERROR: Mismatched tz & offset in index ", i))
    } else if (tz1 == "Etc/GMT+4" & !str_detect(IMP$metashort$timestamp[i], "-0400")) {
      stop(paste0("ERROR: Mismatched tz & offset in index ", i))
    }
  }
  
  for (i in 1:length(files_to_read)) {
    if (str_detect(files_to_read[i], "code")) {
      codes <- read.csv(files_to_read[i])
    } else {
      raw <- read.csv(files_to_read[i]) %>%
        mutate(date_time = anytime(timenum, tz = tz1)) 
    }
    if (i > 2) {
      stop("ERROR: more files than expected")
    }
  }
  
  # Wrangle step data
  segmenting_steps <- 
    IMP$metashort %>% 
    separate(timestamp, into = c("day", "time"), sep = c("T")) %>% 
    separate(time, into = c("time", "omit"), sep = c("-")) %>% 
    select(-omit) %>% 
    unite("date_time", c("day", "time"), sep = " ", remove = FALSE) %>%
    mutate(date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%OS", tz = tz1)) %>%
    separate(time, c("hour", "minute", "second"), sep = ":", remove = FALSE) %>%
    group_by(day, hour, minute) %>%
    nest() %>%
    mutate(sum_steps_vig = map_dbl(.x = data, ~sum(.x$steps_vig)),
           sum_steps_nVig = map_dbl(.x = data, ~sum(.x$steps_nVig))) %>%
    unnest(cols = c(data)) %>%
    filter(second == "00") 
  
  # Merge with activity types
  # Code as vig if > 5 minute-bput of MVPA; otherwise, code as not vig
  comb <-
    raw %>%
    left_join(codes) %>% 
    left_join(segmenting_steps) %>%
    # sum if engaged in bout of PA for > 5 minutes | MVPA for > 1 minute
    mutate(sum_steps_final = if_else((str_detect(class_name, "LIG|MOD|MVPA") &
                                        str_detect(class_name, "bts_5|bts_10")) |
                                       str_detect(class_name, "MVPA_bts_1_5"), 
                                     sum_steps_vig, sum_steps_nVig))
  
  # Make sure merge worked as expected
  if (nrow(raw) != nrow(comb)) {
    stop("ERROR: nrows mismatched; double check left join")
  }
  
  # output raw file with step-counts & activity types broken down by minute
  minute_lvl_steps <-
    comb %>% select(date_time, day, time, class_name, anglez, ENMOa, contains("sum_steps")) 
  attr(minute_lvl_steps$date_time, "tzone") <-"America/New_York"
  minute_lvl_steps$time <- strftime(minute_lvl_steps$date_time, format="%T")
  
  if (!dir.exists(paste0(path, appendPost, as.character(pat_ID)))) {
    dir.create(paste0(path, appendPost, as.character(pat_ID)))
  } else {
    print("Warning: Directory already exists")
  }
  write.csv(minute_lvl_steps %>% select(-sum_steps_vig, -sum_steps_nVig) %>% mutate(TZ = tz(date_time)),
            paste0(path, appendPost, as.character(pat_ID), "/minute_lvl_activity.csv"),
            row.names = FALSE)
  
  # Day summaries (code book: https://cran.r-project.org/web/packages/GGIR/vignettes/GGIR.html, 4.3.1)
  read.csv(list.files(path = paste0(path, appendPre, as.character(pat_ID), "/results"),
                      pattern = "part5_daysummary*", full.names = TRUE)) %>%
    # select variables most likely to be of use in analysis
    select(ID:sleep_efficiency, contains("L5TIME"), contains("M5TIME"), daytype) %>%
    # join with step data
    left_join(comb %>%
                group_by(day) %>%
                summarise(sum_steps_final = sum(sum_steps_final)),
              by = c("calendar_date" = "day")) %>%
    mutate(actigraphy_TZ = tz1) %>%
    # output to CSV
    write.csv(., paste0(path, appendPost, as.character(pat_ID), "/person_lvl_summary.csv"),
              row.names = FALSE)
  
  # Summary stats & plots
  QC_plot <-
    comb %>%
    select(date_time, day, class_name, contains("sum_steps_"), ENMOa) %>%
    gather(key = key, value = value, -date_time, -class_name, -day) %>%
    ggplot(., aes(x = date_time, y = value)) +
    geom_point(alpha = 0.2) +
    geom_line(alpha=0.3) +
    geom_point(y = 0, aes(color = class_name), shape = 124) +
    facet_wrap(~key, scales = "free", nrow = 4) +
    theme_bw() +
    labs(color = "Class", y = "Steps", 
         title = "Evaluating actigraphy-derived step counts across different activity types",
         subtitle = "Note: y-axis scales differ across subplots") +
    guides(color = guide_legend(override.aes=list(shape = 15))) +
    scale_x_datetime(date_breaks = "1 day")
  
  if (output_plot) {
    ggsave(paste0(path, appendPost, as.character(pat_ID), "/QC_plot.tiff"), 
           QC_plot, width = 12, height = 6)
  }
}






