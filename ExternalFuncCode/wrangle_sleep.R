# -----------------------------
# Actigraphy pre-processing: Wrangling sleep from EMA
# -----------------------------

wrangle_sleep <- function (path, patID, EMA_data, replaceNA = FALSE) {
  #' create sleep log from EMA data, to be called in Actigraphy_preFunc.R
  #' 
  #' @details
  #' saves sleep log as CSV to Raw/ID_patID folder
  #' returns list containing: 
  #' 1. # nights with overlapping EMA + actigraphy data 
  #' 2. initiate processing when EMA starts, measured in hours after start of actigraphy data collection  
  #' 3. terminate processing when EMA ends, measured in hours before the end of actigraphy data collection
  #' 
  #' @param path string indicating file pathway to folder containing Raw & Processed data folders
  #' @param patID participant ID, strictly numeric
  #' @param EMA_data data-frame; refer to PrePostProcessing_Demo.R for SQL query
  #' @param replaceNA whether to replace NAs in sleep log with widest empirical window, i.e., min bed to max wake
  #' 
  if (!require("readr")) install.packages("readr"); require ("readr")
  if (!require("chron")) install.packages("chron"); require ("chron")
  if (!require("remotes")) install.packages("remotes"); require ("remotes")
  if (!require("tidyverse")) remotes::install_github("tidyverse/tidyverse@v1.3.1"); require ("tidyverse")
  if (!require("GENEAread")) remotes::install_github("cran/GENEAread@v2.0.9"); require ("GENEAread")
  
  # identify min & max date range, including onboarding EMA
  dates_times_bounding_nights <- 
    EMA_data %>% 
    select(start_time) %>%
    mutate(start_time = as.POSIXct(start_time)) %>% 
    summarise(min = min(start_time),
              max = max(start_time)) %>%
    gather(key = key, value = value) %>% 
    pull(value) 
  dates_bounding_nights <- as.Date(dates_times_bounding_nights)
  
  # Extract self-report sleep & wake data from EMA morning battery 
  EMA_data_tmp <- 
    EMA_data %>%
    mutate(start_time = as.POSIXct(start_time, tz = "America/New_York")) %>% 
    filter(str_detect(study_name, "morn") & test_id == 220) %>%
    select(user_id, sitting_id, start_time, data) %>%
    mutate(data = str_remove_all(data, "[{}\"\\[\\]]")) %>%
    separate_rows(data, sep = ",") %>%
    separate(data, into = c("key", "value"), sep = ":", convert = TRUE) %>%
    filter(str_detect(key, "wake|sleep") & str_detect(key, "time$") & !str_detect(key, "timeOnPage")) %>%
    mutate(x = as.numeric(value),
           x = as.POSIXct(paste(floor(x), round((x-floor(x))*60), sep=":"), format="%H:%M")) %>%
    separate(x, into = c("omit", "SW_time"), sep = " ") %>%
    # Ensuring only 1 morning EMA per day
    separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>%
    arrange(start_time) %>%
    mutate(start_time_lag = lag(start_time),
           diff = as.numeric(difftime(start_time, start_time_lag, units = "hours")),
           track = if_else(key == "go_to_sleep_time" & diff < 12, 1, 0)) %>%
    ungroup() %>% group_by(sitting_id) %>% nest() %>%
    mutate(sum = map_dbl(.x = data, ~sum(.x$track))) %>%
    unnest(data) %>% ungroup() %>%
    filter(is.na(sum) | sum == 0) %>%
    select(-omit, -value, -time, -start_time, -start_time_lag, -diff, -track, -sum)
  
  # Expand grid to include nights where participant didn't provide EMA data
  study_dates_sleepWake <- 
    seq(from = dates_bounding_nights[1], to = dates_bounding_nights[2], by = "day") %>%
    data.frame(date = .) %>%
    expand_grid(., key = unique(EMA_data_tmp$key)) %>%
    left_join(EMA_data_tmp %>% mutate(date = as.Date(date)))
  
  # excerpt min & max dates from actigraphy bin file
  # code expects one bin file per raw folder
  binfile <- read.bin(list.files(path = paste0(path, "Raw/ID_", as.character(patID)),
                                 pattern = "*.bin", full.names = TRUE))
  bin_date_times <- c(min(binfile$page.timestamps), max(binfile$page.timestamps))
  attr(bin_date_times, "tzone") <-"America/New_York"
  bin_dates <- as.Date(bin_date_times)
  rm(binfile)
  
  # Align EMA & actigraphy dates
  if ((bin_date_times[2] < dates_times_bounding_nights[1]) | 
      (dates_times_bounding_nights[2] < bin_date_times[1])) {
    stop("ERROR: EMA & binfile dates don't overlap")
  } 
  # # No longer relevant - now, outputting sleep log for all days with actigrahpy data
  # if (bin_date_times[1] < dates_times_bounding_nights[1]) {
  #   # scenario where actigraphy starts recording prior to EMA onboarding
  #   hrs_del_start <- as.numeric(dates_times_bounding_nights[1] - bin_date_times[1], units = "hours")
  # } else {
  #   hrs_del_start <- 0
  # }
  # if (bin_date_times[2] > dates_times_bounding_nights[2]) {
  #   # scenario where EMA protocol stops while actigraphy is still recording
  #   hrs_del_end <- as.numeric(bin_date_times[2] - dates_times_bounding_nights[2], units = "hours")
  # } else {
  #   hrs_del_end <- 0
  # }
  
  # create df with two rows per day of actigraphy collection
  actigraphy_df <- 
    data.frame(date = seq(bin_dates[1], bin_dates[2], by = "day")) %>%
    mutate(var_1 = "go_to_sleep_time", var_2 = "wake_up_time") %>%
    gather(key = omit, value = key, -date) %>% arrange(date) %>% select(-omit)
  
  # create sleep log
  EMA_data_morn <-
    actigraphy_df %>% 
    left_join(study_dates_sleepWake, by = c("date", "key")) %>%
    mutate(ID = parse_number(user_id),
           key = str_remove_all(key, "_")) %>%
    group_by(date, sitting_id) %>% nest() %>% 
    ungroup() %>%
    mutate(Night = as.character(row_number())) %>%
    unnest(c(data)) %>%
    # mutate(Night = as.character(as.numeric(Night) - 1)) %>%
    filter(Night != 0) %>%
    {{. ->> tmp}} %>%
    unite(key, c("key", "Night"), sep = "_") %>%
    select(ID, key, SW_time) %>%
    mutate(ID = mean(ID,na.rm = TRUE)) %>%
    pivot_wider(names_from = key, values_from = SW_time)
  
  n_nights <- tmp %>% filter(is.na(user_id)) %>% 
    ungroup() %>% summarise(max = max(as.numeric(Night))) %>% pull(max)
  
  if (replaceNA) {
    extreme_vals <- 
      EMA_data_morn %>%
      gather(key = key, value = value, -ID) %>% 
      mutate(value = chron(times = value)) %>%
      mutate(key2 = if_else(str_detect(key, "wake"), "max_wake", "min_bed")) %>%
      group_by(key2) %>%
      summarise(min = min(value, na.rm = TRUE),
                max = max(value, na.rm = TRUE)) %>%
      mutate(val = if_else(str_detect(key2, "wake"), max, min)) %>%
      select(-min, -max)
    
    EMA_data_morn <- 
      EMA_data_morn %>% 
      gather(key = key, value = value, -ID) %>%
      mutate(key2 = if_else(str_detect(key, "wake"), "max_wake", "min_bed")) %>%
      left_join(extreme_vals) %>%
      mutate(value = if_else(is.na(value), as.character(val), value)) %>%
      select(-key2, -val) %>%
      pivot_wider(names_from = key, values_from = value)
  }
  
  if ((EMA_data_morn %>% pull(ID)) != patID) {
    stop("patID argument does not match EMA_data")
  } 
  
  # per GGIR Google group, want empty cell (rather than NA) in CSV 
  write.table(EMA_data_morn, 
              paste0(path, "Raw/ID_", as.character(pat_ID), "/sleepLog", as.character(pat_ID), ".csv"),
              na = "",
              row.names = FALSE,
              col.names = TRUE,
              append = FALSE,
              sep = ",")
  
  return(n_nights)
}

