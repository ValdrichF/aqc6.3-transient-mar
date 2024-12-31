# Split the RUN files into 11 groups and create a batch file to run each group sequentially
library(dplyr)
library(stringr)
setwd("E:/AMIGO_needed/RUNFILES/scenarios")
all_runs = list.files(pattern = ".*(RUN)|(run)$")

n_splits = 9
splits = rep(1:n_splits, length.out = length(all_runs))%>%
  sort%>%
  as.factor()

splits_dir = "splits"

folders = splits_dir%>%
  paste(1:n_splits, sep = "/")

for (folder in folders){
  if (dir.exists(folder)) unlink(folder, recursive = T)
  dir.create(folder, recursive = T)
}

imod = "D:\\iMOD5\\bin\\iMODFLOW-METASWAP_V5_5.exe"

for (i in 1:n_splits){
  this_split = all_runs[splits==i]
  sans_ext = tools::file_path_sans_ext(this_split)
  print(this_split)
  
  file.copy(this_split,
            paste0(splits_dir, "/", i, "/", this_split))
  if (FALSE){#(i <= 2){
    
    sapply(this_split, function(f){
      path = paste0(splits_dir, "/", i, "/", f)
      run = readLines(path)
      run[1] = str_replace(run[1], "G:", "E:")
      writeLines(run, path)
    })
  }
  
  bat = c(paste0(imod, " " , this_split, " > log_", sans_ext, ".txt"),
          "pause")
  
  bat_name = paste0(splits_dir, "/", i,  "/", "run_split.bat")
  writeLines(bat, bat_name)
}







