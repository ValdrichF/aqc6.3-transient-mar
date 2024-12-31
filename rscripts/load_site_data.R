# Script to load the data and libraries
# Made so that it can be called from another script using 'source()'

library(dplyr)
library(data.table)
library(lubridate)
library(ggplot2)
library(tagger)
library(cowplot)
library(tidyverse)
library(tidymodels)
library(purrr)
library(furrr)
library(DALEX)
source("functions.R")
big_grid = fread("ortho_grid.csv")
df = fread("response_sites_thresh0.csv")


setnames(big_grid, "area", "area_site")
df[, ":="(time=ymd(time),
          dx = NULL,
          dy = NULL,
          area = as.numeric(area))]

# Merge df with big_grid, to be used for grouping the responses for each recharge site
setkey(df, scenario, rch_site, time)
setkey(big_grid, scenario, rch_site)
df = big_grid[df]
decay_total = df%>%
  nest(.by = c(x, y, scenario, rch_site, recharge, side, area_site))%>%# side, recharge, area_site, 
  mutate(decay_model_total = map(data, fit_exp_decay_N0, y_var = "total"),
         decay_model_water = map(data, fit_exp_decay_N0, y_var = "total_water"))%>%
  mutate(slope_total = map(decay_model_total, ~coef_exp(.)[2]),
         slope_water = map(decay_model_water, ~coef_exp(.)[2]),
         
         initial_total = map(data, function(.dt){
           .dt[time==ymd("2012-02-25"), total]
         }),
         initial_water = map(data, function(.dt){
           .dt[time==ymd("2012-02-25"), total_water]
         }))%>%
  unnest(c(slope_total, slope_water, initial_total, initial_water))%>%
  mutate(slope_total = slope_total*-1,
         slope_water = slope_water*-1)

decay_total_sub = decay_total%>%
  select(-data, -decay_model_total, -decay_model_water)
decay_total_sub = as.data.table(decay_total_sub)
decay_total_sub[,":="(slope_total = log10(slope_total),
                      slope_water = log10(slope_water),
                      initial_total = log10(initial_total),
                      initial_water = log10(initial_water))]

rm(big_grid, df)


load_site_data = function(properties_path, df = decay_total_sub){
    site_properties = fread(properties_path)
    df = copy(df)
    
    # Align the x-y coordinates to the nearest cell
    # And merge with the relevant columns from decay_total
    min_x = site_properties[,min(x)]
    min_y = site_properties[,min(y)]
    df[,":="(
        x = round((x-min_x)/25)*25+min_x,
        y = round((y-min_y)/25)*25+min_y
    )]
    setkey(df, x, y)
    setkey(site_properties, x, y)
    df = site_properties[df]
    # df[,":="(total_recharge = recharge*area_site)]
    
    rm(site_properties)
    
    return(df)
    
}
