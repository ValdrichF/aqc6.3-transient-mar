library(png)
library(recolorize)
library(dplyr)
# setwd("D:/Modeldatabase/Baaksebeek/Res_maps/decay_total/pngs")
dir = "C:/Users/valdf/MODFLOW_transient_res/"
fs = list.files(dir, pattern = "*.png$", full.names = T)

lapply(fs, function(i){
  paste0("working on: ", i)%>%
    print()
  f = readPNG(i)
  
  totals_y = apply(f, 2, sum)
  totals_x = apply(f, 1, sum)
  min_y = (max(totals_y)-totals_y)%>%
    cumsum()
  min_y = which(min_y==0)%>%
    max()
  min_x = (max(totals_x)-totals_x)%>%
    cumsum()
  min_x = which(min_x==0)%>%
    max()
  
  max_y = (max(totals_y)-totals_y)%>%
    rev()%>%
    cumsum()%>%
    rev()
  max_y = which(max_y==0)%>%
    min()
  max_x = (max(totals_x)-totals_x)%>%
    rev()%>%
    cumsum()%>%
    rev()
  max_x = which(max_x==0)%>%
    min()
  
  min_x = max(1, min_x-30)
  max_x = min(nrow(f), max_x+30)
  min_y = max(1, min_y-30)
  max_y = min(ncol(f), max_y+30)
  
  f_crop = f[min_x:max_x, min_y:max_y,]
  
  plotImageArray(f_crop)
  writePNG(f_crop, i)
  
})

  

