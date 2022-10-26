# ActigraphyPipeline
Repository contains code to process GeneActiv actigraphy data and merge output with ecological momentary assessment (EMA) data.   

**Description of folders & files:**  

1. PrePostProcessingDemo.R: Run pipeline      
2. ExternalFuncCode: Folder containing functions required to run pipeline  
   a. wrangle_sleep.R: create sleep log from EMA data  
   b. actigraph_preFunc.R: perform actigraphy pre-processing on raw bin file + EMA sleep log
   c. actigraph_postP.R: perform actigraphy post-processing
   d. parse_activity.R: create log with step and/or activity data linked to EMA sitting ID
   e. verisense_count_steps.R & myscript.R: compute step counts (called in Actigraphy_preFunc.R)  

**Additional notes:**  

1. Sleep and activity metrics derived using GGIR (Migueles JH, Rowlands AV, et al., 2019)  
2. Step counts estimated using open source code (Matthew Patterson, Verisense)  
3. Informal validity checks on step count estimates performed by computing MAPE relative to commerically available wearables (Apple Watch, Garmin Forerunner)  
