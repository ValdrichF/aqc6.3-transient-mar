# Script to read a big csv and save a copy with only the interesting cells
library(data.table)

big_grid = fread("ortho_grid.csv", select = c("x", "y"))

site_properties_files = list.files("C:/Users/valdf/MODFLOW_transient_res",
                                   "site_properties.*\\d+.csv",
                                   full.names = T)

sample = fread(site_properties_files[1], select = c("x", "y"))

# Align the x-y coordinates to the nearest cell
min_x = sample[,min(x)]
min_y = sample[,min(y)]
big_grid[,":="(
  x = round((x-min_x)/25)*25+min_x,
  y = round((y-min_y)/25)*25+min_y
)]

rm(sample)
gc()
# Loop through all the files & save the subset
setkey(big_grid, x, y)

lapply(site_properties_files,
       function(path){
         print(paste0("Working on: ", basename(path)))
         site_properties = fread(path)
         site_properties[, ":="(
           V1=NULL,
           dx = NULL,
           dy = NULL
           )]
         setkey(site_properties, x, y)
         
         site_properties = site_properties[big_grid]
         
         fwrite(site_properties, basename(path))
         return(1)
       })
