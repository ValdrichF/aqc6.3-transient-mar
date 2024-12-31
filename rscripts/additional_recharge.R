# Add the additional recharge to all the stress periods in the RUN file.
library(stringr)
library(lubridate)
library(dplyr)

# setwd("D:/iMOD5/IMOD_USER/RUNFILES")
setwd('E:/RUNFILES')
# Read the Run file
original_run_path = "scenario_weekly.run"
original_run = readLines(original_run_path)
scenario_rch_path = "D:\\\\Modeldatabase\\\\AMIGO31\\\\rch_transient\\\\rch_1_l1.idf"

# Create the character vector that I'd like to append to the run file
rch_line = c("  1,(RCH) Recharge package",
             "     1,   1.000000    ,   0.000000    ,RCH_PATH",
             "     1,   1.000000    ,   0.000000    ,0.000000")

# Add RCH package as an active module
original_run = c(original_run[1:19], 
                 "1,0 (RCH)",
                 original_run[20:length(original_run)])

# Modify the RCH_PATH
rch_line = rch_line%>%
  str_replace("RCH_PATH", scenario_rch_path)

# Find the stress periods
indexes = str_which(original_run, "20\\d{6},1,0,20\\d{6}$")%>%
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


# Loop through each rch scenario and create a new run file
run_tmp = character(length(original_run)+nrow(indexes)*2)
run_tmp[indexes$corrected_ind] = rch_line[1]
run_tmp[indexes$corrected_ind[indexes$winter]+1] = rch_line[2]
run_tmp[indexes$corrected_ind[!indexes$winter]+1] = rch_line[3]
run_tmp[c(indexes$corrected_ind, indexes$corrected_ind+1)*-1] = original_run

# Adjust the SHD to the closest head from the baseline scenario
heads_dir = "E:\\MODFLOW_transient_res\\baseline_avg_weekly\\head"
heads = list.files(heads_dir)
heads_dates = str_extract(heads, "\\d{8}")%>%
  unique()%>%
  ymd()
start_date = "20110328"
SHD_date = ymd(start_date) - heads_dates
SHD_date[SHD_date<0] = 9999 # ignore CHDs after the start date
head_last = heads_dates[which(SHD_date==min(SHD_date, na.rm = T))]%>%
  format("%Y%m%d")

SHD_ind = str_which(run_tmp, "015..SHD")+1
SHD_ind = seq(SHD_ind, length.out = 15)
run_tmp[SHD_ind] = str_replace(run_tmp[SHD_ind], "D:.*head\\\\", 
                               heads_dir%>%
                                 str_replace_all("\\\\","\\\\\\\\")%>%
                                 paste0("\\\\"))
run_tmp[SHD_ind] = str_replace(run_tmp[SHD_ind], "\\d{8}", head_last)

writeLines(run_tmp, original_run_path)
