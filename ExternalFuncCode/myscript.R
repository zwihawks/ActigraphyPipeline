# update to local path to verisense_count_steps.R file
# Files from GitHub, algorithm based on Gu et al. (2017)
# https://github.com/ShimmerEngineering/Verisense-Toolbox/tree/master/Verisense_step_algorithm 
source("ExternalFuncCode/verisense_count_steps.R") 

# for windows laptop (after setting working directory to GGIR/)
# source("ExternalFuncCode/verisense_count_steps.R") 

vcs_fun =  list(FUN=verisense_count_steps,
              parameters= c(3, 4, 15, -5.0, 3, 4, 0.001, 1.2,
                            3, 5, 15, -0.5, 3, 4, 0.001, 1.2),
              expected_sample_rate= 15,
              expected_unit="g",
              colnames = c("steps_vig", "steps_nVig"),
              outputres = 1,
              minlength = 1,
              outputtype="numeric",
              aggfunction = sum,
              timestamp=F,
              reporttype="event")

