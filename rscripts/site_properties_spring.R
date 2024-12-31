# Script to process the site_properties
# Calculate the average during spring
library(data.table)

site_properties = fread("D:/Dropbox/WUR/Rscripts/paper_2/site_properties_spring_gaussian_blur_500.csv")
site_properties[, ":="(time=ymd(time))]

site_properties = site_properties[
  time>=ymd("20120225"),
  .(transmissivity = transmissivity[1],
    resistance = resistance[1],
    depth = mean(head(depth[8])),
    storage_coefficient = storage_coefficient[1],
    
    river_stage = mean(head(river_stage,8)),
    river_stage_ghg = river_stage_ghg[1],
    river_stage_zeros = mean(head(river_stage_zeros,8)),
    river_stage_ghg_zeros = river_stage_ghg_zeros[1],
    river_stage_depth = mean(head(river_stage_depth,8)),
    river_stage_maxdepth = mean(head(river_stage_maxdepth,8)),
    
    river_conductance = river_conductance[1],
    river_conductance_ghg = river_conductance_ghg[1],
    river_conductance_zeros = river_conductance_zeros[1],
    
    drain_binary = drain_binary[1],
    drain_level = mean(head(drain_level,8)),
    drain_level_maxdepth = mean(head(drain_level_maxdepth,8)),
    drain_level_depth = mean(head(drain_level_depth,8)),
    
    drain_conductivity = drain_conductivity[1],
    drain_conductivity_zeros = drain_conductivity_zeros[1],
    
    river_density = mean(head(river_density,8)),
    river_density_ghg = river_density_ghg[1]),
  .(x, y)
]
