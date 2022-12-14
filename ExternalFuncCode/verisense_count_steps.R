verisense_count_steps <- function(input_data=runif(500,min=-1.5,max=1.5), coeffs=c(0,0,0)) {
  # by Matthew R Patterson, mpatterson@shimmersensing.com
  # Find peaks of RMS acceleration signal according to Gu et al, 2017 method
  # This method is based off finding peaks in the summed and squared acceleration signal
  # and then using multiple thresholds to determine if each peak is a step or an artifact.
  # An additional magnitude threshold was added to the algorithm to prevent false positives 
  # in free living data. 
  #
  # returns sample location of each step
  # -------------
  # ZH mod:
  # rep 1 = shorter min period & more negative similarity threshold for vigorous/continuous activity
  # refer to Gu et al. 2017 & Verisense Toolbox validation paper on GitHub
  # -------------
  for (rep in c(1, 2)) {
    fs = 15 # temporary for now, this is manually set
    acc <- sqrt(input_data[,1]^2 + input_data[,2]^2 + input_data[,3]^2)
    
    if (sd(acc) < 0.025) {
      # acceleration too low, no steps
      num_seconds = round(length(acc) / fs)
      steps_per_sec = rep(0,num_seconds)
    } else {
      # Search for steps
      # Thresholds
      if (rep == 1){
        k <- coeffs[[1]]
        period_min <- coeffs[[2]]
        period_max <- coeffs[[3]]
        sim_thres <- coeffs[[4]]   # similarity threshold
        cont_win_size <- coeffs[[5]]  # continuity window size
        cont_thres <- coeffs[[6]]     # continuity threshold
        var_thres <- coeffs[[7]]  # variance threshold
        mag_thres <- coeffs[[8]]
      } else if (rep == 2) {
        k <- coeffs[[9]]
        period_min <- coeffs[[10]]
        period_max <- coeffs[[11]]
        sim_thres <- coeffs[[12]]   # similarity threshold
        cont_win_size <- coeffs[[13]]  # continuity window size
        cont_thres <- coeffs[[14]]     # continuity threshold
        var_thres <- coeffs[[15]]  # variance threshold
        mag_thres <- coeffs[[16]]
      } else {
        print("ERROR")
      }
      
      # find the peak rms value is every range of k
      half_k <- round(k/2)
      segments <- floor(length(acc) / k)
      peak_info <- matrix(NA,nrow=segments,ncol=5)
      
      # for each segment find the peak location
      for (i in 1:segments) {
        start_idx <- (i-1) * k + 1
        end_idx <- start_idx + (k-1)
        tmp_loc_a <- which.max(acc[start_idx:end_idx])
        tmp_loc_b <- (i-1) * k + tmp_loc_a
        # only save if this is a peak value in range of -k/2:+K/2
        start_idx_ctr <- tmp_loc_b - half_k
        if (start_idx_ctr < 1) {
          start_idx_ctr <- 1
        }
        end_idx_ctr <- tmp_loc_b + half_k
        if (end_idx_ctr > length(acc)) {
          end_idx_ctr <- length(acc)
        }
        check_loc <- which.max(acc[start_idx_ctr:end_idx_ctr])
        if (check_loc == (half_k + 1)) {
          peak_info[i,1] <- tmp_loc_b
          peak_info[i,2] <- max(acc[start_idx:end_idx])
        }
      }
      peak_info <- peak_info[is.na(peak_info[,1])!=TRUE,] # get rid of na rows
      
      # filter peak_info[,2] based on mag_thres
      peak_info <- peak_info[peak_info[,2] > mag_thres,]
      if (length(peak_info) > 10) {  # there must be at least two steps
        num_peaks <- length(peak_info[,1])
        
        no_steps = FALSE
        if (num_peaks > 2) {
          # Calculate Features (periodicity, similarity, continuity)
          peak_info[1:(num_peaks-1),3] <- diff(peak_info[,1]) # calculate periodicity
          peak_info <- peak_info[peak_info[,3] > period_min,] # filter peaks based on period_min
          peak_info <- peak_info[peak_info[,3] < period_max,]   # filter peaks based on period_max 
        } else {
          no_steps = TRUE
        }
      } else {
        no_steps = TRUE
      }
      
      if ( length(peak_info)==0 || length(peak_info) == sum(is.na(peak_info)) || no_steps == TRUE) {
        # no steps found
        num_seconds = round(length(acc) / fs)
        steps_per_sec = rep(0,num_seconds)
      } else {
        # calculate similarity
        num_peaks <- length(peak_info[,1])
        peak_info[1:(num_peaks-2),4] <- -abs(diff(peak_info[,2],2)) # calculate similarity
        peak_info <- peak_info[peak_info[,4] > sim_thres,]  # filter based on sim_thres
        peak_info <- peak_info[is.na(peak_info[,1])!=TRUE,] # previous statement can result in an NA in col-1
        
        # calculate continuity
        if (length(peak_info[,3]) > 5) {
          end_for <- length(peak_info[,3])-1
          for (i in cont_thres:end_for) {
            # for each bw peak period calculate acc var
            v_count <- 0 # count how many windows were over the variance threshold
            for (x in 1:cont_thres) {
              if (var(acc[peak_info[i-x+1,1]:peak_info[i-x+2,1]]) > var_thres) {
                v_count = v_count + 1
              }
            }
            if (v_count >= cont_win_size) {
              peak_info[i,5] <- 1 # set continuity to 1, otherwise, 0
            } else {
              peak_info[i,5] <- 0
            }
          }
        } 
        peak_info <- peak_info[peak_info[,5]==1,1] # continuity test - only keep locations after this
        peak_info <- peak_info[is.na(peak_info)!=TRUE] # previous statement can result in an NA in col-1
        
        if (length(peak_info)==0) {
          # no steps found
          num_seconds = round(length(acc) / fs)
          steps_per_sec = rep(0,num_seconds)
        } else {
          
          # for GGIR, output the number of steps in 1 second chunks
          start_idx_vec <- seq(from=1,to=length(acc),by=fs)
          steps_per_sec <- table(factor(findInterval(peak_info, start_idx_vec), levels = seq_along(start_idx_vec)))
          steps_per_sec <- as.numeric(steps_per_sec)
        }
      }
    }
    
    if (rep == 1) {
      steps_per_sec_vec <- steps_per_sec
    } else if (rep == 2) {
      steps_per_sec_vec <- cbind(steps_per_sec_vec, steps_per_sec)
    } else {
      print("ERROR")
    }
  }
  
  return(steps_per_sec_vec)
}