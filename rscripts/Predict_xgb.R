#########
# Estimate the intial response and the decay coefficient for each cell in the catchment
# Assume a recharge rate of 15 mm/day over 1 sq km
# Test it out using a XGBoost models and U-NET for the the initial response
# 



###########

library(stringr)
library(lubridate)
library(tidyverse)
library(ggplot2)
library(tagger)
library(data.table)
library(ggplot2)
library(DALEX)
library(stringr)
library(tagger)
library(shapviz)
library(readr)
library(patchwork)
source("../paper_2/functions.R")
common_theme = list(
    theme_light(),
    theme(text=element_text(size=12), 
          axis.text=element_text(size=8, color = 'black'),
          axis.title=element_text(size=12),
          plot.title=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12), 
          legend.position = 'bottom',
          strip.background = element_blank(),
          strip.text = element_blank(),
          tagger.panel.tag.text = element_text(size = 8, color = 'black', face = 'bold'),
          tagger.panel.tag.background = element_rect(fill = "white", color = 'grey')
    )
)

# Load the models
model_ini = read_rds("D:/Dropbox/WUR/Rscripts/paper_2/models_response_ini2/xgboost_response_ini2_.rds")
model_decay = read_rds("D:/Dropbox/WUR/Rscripts/paper_2/models_response_decay2/xgboost_response_decay2_.rds")
# Step 1: Load the new dataset from the CSV file
new_data_winter <- read_csv("C:/Users/valdf/MODFLOW_transient_res/site_properties_winter_gb_500.csv")
new_data_summer <- read_csv("C:/Users/valdf/MODFLOW_transient_res/site_properties_spring_gb_500.csv")

# Step 2: Preprocess the new data 
# select the columns that were used to train the model
# remove na drain_conductivity
# add recharge and area columns
new_data_winter = new_data_winter%>%
    select(any_of(c("x", "y", names(model_ini$training), "sto_baseline")))%>%
    replace_na(list(drain_conductivity_zeros = 0))%>%
    mutate(recharge = 15)%>% # same as paper 1: 15mm/day
    mutate(area_site = 10*1e4) # same as paper 1: 10 ha
new_data_summer = new_data_summer%>%
    select(any_of(c("x", "y", names(model_decay$training), "sto_baseline")))%>%
    replace_na(list(drain_conductivity_zeros = 0))%>%
    mutate(recharge = 15)%>% # same as paper 1: 15mm/day
    mutate(area_site = 10*1e4) # same as paper 1: 10 ha

# Step 3: Predict using the trained model
start_time  = Sys.time()
predictions_ini <- predict(model_ini, new_data_winter)
end_ini = Sys.time()
predictions_decay <- predict(model_decay, new_data_summer)
end_decay = Sys.time()

time_ini = end_ini - start_time # 15.21642 seconds
time_decay = end_decay - end_ini # 18.05426 seconds

# Step 4: Attach predictions to your new dataset if you want
predictions_ini_winter <- new_data_winter %>%
    mutate(predicted_ini = 10^predictions_ini$.pred)%>%
    as.data.table()
predictions_decay_summer <- new_data_winter %>%
    mutate(predicted_decay = 10^predictions_decay$.pred)%>%
    as.data.table()

# merge into a single dataset
setkey(predictions_decay_summer, x, y)
setkey(predictions_ini_winter, x, y)
predictions = predictions_decay_summer[predictions_ini_winter[, .(x, y, predicted_ini)]]

# save predictions to a CSV
write_csv(predictions_ini_winter, "C:/Users/valdf/MODFLOW_transient_res/predictions_ini_winter.csv")
write_csv(predictions_decay_summer, "C:/Users/valdf/MODFLOW_transient_res/predictions_decay_summer.csv")
write_csv(predictions, "C:/Users/valdf/MODFLOW_transient_res/predictions_output.csv")












