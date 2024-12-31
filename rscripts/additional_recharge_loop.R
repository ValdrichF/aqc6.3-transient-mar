# Add the additional recharge to all the stress periods in the RUN file.
library(stringr)
library(lubridate)
library(dplyr)

# setwd("D:/iMOD5/IMOD_USER/RUNFILES")
setwd('E:/AMIGO_needed')
# Read the Run file
original_run_path = "RUNFILES/avg_baseline.run"
rch_dir = "rch_transient"
save_dir = "RUNFILES/scenarios"

original_run = readLines(original_run_path)
rch_files = paste(getwd(),rch_dir, sep = "/")%>%
  list.files(".*idf", full.names = T, include.dirs = T)%>%
  str_replace_all("/","\\\\\\\\")

if(dir.exists(save_dir)) unlink(save_dir, recursive = T)
dir.create(save_dir)

# Create the character vector that I'd like to append to the run file
rch_line = c("  1,(RCH) Recharge package",
             "     1,   1.000000    ,   0.000000    ,RCH_PATH",
             "     1,   1.000000    ,   0.000000    ,0.000000")

# Add RCH package as an active module
original_run = c(original_run[1:19], 
                 "1,0 (RCH)",
                 original_run[20:length(original_run)])

worked = lapply(rch_files, function(scenario_rch_path){
  # Modify the RCH_PATH
  rch_line = rch_line%>%
    str_replace("RCH_PATH", scenario_rch_path)
  
  # Find the stress periods
  indexes = str_which(original_run, "20\\d{6},1,1,20\\d{6}$")%>%
    data.frame(indexes=.)
  
  # Check which are in the winter period
  indexes = indexes%>%
    mutate(date = original_run[indexes]%>%
             str_extract("20\\d{6}")%>%
             ymd())%>%
    mutate(month = month(date),
           year = year(date))%>%
    mutate(winter = month%in%c(1,2,10,11,12)&year%in%c(2011,2012))
  
  # Correct the indexes for getting shifted below cause of additional lines
  indexes = indexes%>%
    mutate(corrected_ind = indexes + seq_along(indexes)*2-1)
  
  run_tmp = character(length(original_run)+nrow(indexes)*2)
  run_tmp[indexes$corrected_ind] = rch_line[1]
  run_tmp[indexes$corrected_ind[indexes$winter]+1] = rch_line[2]
  run_tmp[indexes$corrected_ind[!indexes$winter]+1] = rch_line[3]
  run_tmp[c(indexes$corrected_ind, indexes$corrected_ind+1)*-1] = original_run
  
  # Change the save directory
  run_tmp[1] = paste0("'K:\\MODFLOW_transient_res\\", 
                      scenario_rch_path%>%basename()%>%tools::file_path_sans_ext(),
                      "'")
  
  filename = scenario_rch_path%>%
    basename()%>%
    tools::file_path_sans_ext()%>%
    paste0(".RUN")%>%
    paste(save_dir, ., sep = "\\\\")
  
  writeLines(run_tmp, filename)
  return(1)
  
})

file.copy(original_run_path, 
          original_run_path%>%
            basename()%>%
            paste(save_dir, ., sep = "/"))

