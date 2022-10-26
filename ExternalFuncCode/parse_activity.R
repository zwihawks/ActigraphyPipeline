# -----------------------------
# Minute-level activity data
# -----------------------------

parse_activity <- function (path, patID, EMA_data, direction = c(-1, 1), offset_min, 
                            steps = TRUE, activity = TRUE, log = FALSE,
                            export = FALSE, return = TRUE) {
  #' create log with step and/or activity data linked to EMA sitting ID
  #' 
  #' @details
  #' saves activity log as CSV to Post/output_ID_patID folder
  #' 
  #' @param path string indicating file pathway to folder containing Raw & Processed data folders
  #' @param patID participant ID, strictly numeric
  #' @param EMA_data data-frame; refer to PrePostProcessing_Demo.R for SQL query
  #' @param direction subset data before or after sitting time_stamp; use values -1 and 1, respectively
  #' @param offset_min array of length >= 1 containing duration in minutes to summarize data around EMA time_stamp
  #' @param steps T/F indicating whether to output step counts within designated window
  #' @param activity T/F indicating whether to output most common activity + percentage within designated window
  #' @param log T/F indicating whether to use sleep data estimated using sleep log
  #' @param export export activity log as CSV
  #' @param return return df as object
  #'
  if (!require("remotes")) install.packages("remotes"); require ("remotes")
  if (!require("tidyverse")) remotes::install_github("tidyverse/tidyverse@v1.3.1"); require ("tidyverse")
  
  if (log) {
    appendPost <- "Processed/Post/Log/output_ID_"
  } else {
    appendPost <- "Processed/Post/noLog/output_ID_"
  }
  
  minute_lvl_activity <-
    read.csv(paste0(path, appendPost, 
                    as.character(patID), "/minute_lvl_activity.csv")) %>%
    mutate(date_time = as.POSIXct(date_time, tz = "America/New_York"))
  
  start_times <- 
    EMA_data %>% 
    filter(!str_detect(battery_name, "Onboarding")) %>%
    select(user_id, sitting_id, start_time) %>%
    mutate(start_time = as.POSIXct(start_time, tz = "America/New_York")) %>%
    unique(.) %>% 
    # retain start times that have complete activity data for all offsets
    filter(start_time <= max(
      max(minute_lvl_activity$date_time),
      max((minute_lvl_activity$date_time + direction*max(offset_min)*60))))
  
  if ((max(minute_lvl_activity$date_time) < min(start_times$start_time)) | 
      (max(start_times$start_time) < min(minute_lvl_activity$date_time))) {
    stop("ERROR: EMA & binfile dates don't overlap")
  } 
  
  if (exists("tmp_comb_steps", envir=globalenv())) {rm(tmp_comb_steps, envir=globalenv())}
  if (exists("tmp_comb_activity", envir=globalenv())) {rm(tmp_comb_activity, envir=globalenv())}
  
  # loop through start times, binding one row to df at a time
  # there's a lot below that isn't efficient... still seems to be running fast enough, for now
  for (i in 1:nrow(start_times)) {
    time <- start_times$start_time[i]
    user_id = start_times$user_id[i]
    sitting_id = start_times$sitting_id[i]
    
    # prep subsetted_activity df, to be modified below
    subsetted_activity <- minute_lvl_activity
    
    # identify rows that fall within windows of interest
    for (j in offset_min) {
      val1 <- time
      val2 <- time + direction*j*60
      
      start <- min(val1, val2)
      end <- max(val1, val2)
      
      subsetted_activity <-
        subsetted_activity %>%
        mutate(tmp = if_else(date_time >= start & date_time < end, 1, 0))
      colnames(subsetted_activity)[colnames(subsetted_activity) == "tmp"] <- paste0("offset_min_", j)
    }
    
    # subset rows that fall within at least one of our offset windows
    subsetted_activity <- 
      subsetted_activity %>%
      filter_at(vars(contains("offset_min")), any_vars(. > 0))
    
    if (nrow(subsetted_activity) == 0) {
      # if no rows fall within offset windows, return to top of loop
      next
    } else {
      # sum steps
      if (steps) {
        
        tmp_steps <- 
          data.frame(
            user_id = user_id,
            sitting_id = sitting_id,
            EMA_time = time,
            direction = direction
          )
        
        # iterate over different offset windows
        for (j in offset_min) {
          sum <- 
            subsetted_activity %>%
            select(sum_steps_final, contains(as.character(j))) %>%
            filter_at(vars(contains(as.character(j))), any_vars(. > 0)) %>%
            summarise(sum = sum(sum_steps_final)) %>% pull(sum)
          
          # Add new variable to tmp_steps df + rename to specify offset
          tmp_steps$steps <- sum
          colnames(tmp_steps)[colnames(tmp_steps) == "steps"] <- paste0("steps_offset", j)
        }
        
        if (!exists("tmp_comb_steps")) {
          tmp_comb_steps <- tmp_steps
        } else {
          tmp_comb_steps <- tmp_comb_steps %>% bind_rows(tmp_steps)
        }
      }
      
      # compute percent time engaged in activities subtypes
      if (activity) {
        
        tmp_activity <- 
          data.frame(
            user_id = user_id,
            sitting_id = sitting_id,
            EMA_time = time,
            direction = direction)
        
        # iterate over different offset windows
        for (j in offset_min) {
          percents <-
            subsetted_activity %>% 
            select(class_name, contains(as.character(j))) %>%
            filter_at(vars(contains(as.character(j))), any_vars(. > 0)) %>%
            mutate(class_name = str_remove(class_name, "_bts_[0-9]*[_]*[0-9]*|_unbt")) %>%
            group_by(class_name) %>%
            summarise(count = n()) %>%
            mutate(sum = sum(count),
                   percent = count/sum*100) %>%
            ungroup() %>%
            slice_max(percent, n = 1, with_ties = FALSE) %>%
            select(class_name, percent)
          
          # Add new variables to tmp_activity df + rename to specify offset
          tmp_activity$top_class <- ifelse(nrow(percents) == 0, NA, percents$class_name)
          tmp_activity$top_class_val <- ifelse(nrow(percents) == 0, NA, percents$percent)
          colnames(tmp_activity)[colnames(tmp_activity) == "top_class"] <- paste0("topClass_offset", j)
          colnames(tmp_activity)[colnames(tmp_activity) == "top_class_val"] <- paste0("topClassPercent_offset", j)
        }
            
        if (!exists("tmp_comb_activity")) {
          tmp_comb_activity <- tmp_activity
        } else {
          tmp_comb_activity <- tmp_comb_activity %>% bind_rows(tmp_activity)
        }
      }
    }
  }
  
  # interested in steps, activity, or both?
  if (steps & activity) {
    x <- left_join(tmp_comb_steps, tmp_comb_activity)
    if ((nrow(x) != nrow(tmp_comb_steps)) | (nrow(x) != nrow(tmp_comb_activity))) {
      stop("ERROR: join failed.")
    }
  } else if (steps) {
    x <- tmp_comb_steps
  } else if (activity) {
    x <- tmp_comb_activity
  }
  
  # return results
  if (export) {
    write.csv(x, 
              paste0(path, appendPost, as.character(patID), "/EMA_lvl_activity.csv"),
              row.names = FALSE)
  }
  if (return) {
    return(x)
  }
}