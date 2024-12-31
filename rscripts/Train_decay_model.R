### Train the model to predict the decay rate 
source("load_site_data.R")
sigma = 500
# geomline of log(total) vs time with abline showing the exponential decay coefficients and the trend line
# focus on 4 sites in one plot split by colour
a = copy(decay_total_sub)

setDT(a)
a[,depth_rank := rank(slope_rlm)]
# select 0.05, 0.33, 0.67, 0.95 percentile by depth to groundwater
a = a[depth_rank%in%round(c(0.05, 0.1, 0.51, 0.75)*.N), .(rch_site, scenario)]
setkey(a, rch_site, scenario)
setkey(decay_total_sub, rch_site, scenario)

decay_total_preds = decay_total%>%
  inner_join(a)%>%
  mutate(total_pred = mapply(fitted_exp,
                             .df = data,
                             decay_coefficient = slope_N0,
                             scale_0 = 1,
                             SIMPLIFY = F))%>%
  unnest(c(total_pred, data))%>%
  as.data.table()

decay_total_preds[, scenario:=factor(scenario,
                                     levels = .SD[order(slope_rlm), unique(scenario)])]


input_cols = c("recharge", "area_site", "transmissivity", "depth",
               "resistance", "river_conductance_zeros",
               "river_density_ghg") 


# fit the model
additional_input = c("river_stage", "storage_coefficient", "drain_conductivity_zeros")

for (in_put in additional_input) {
  
  paste0("working on: ", in_put)%>%
    print()
  
  cond = paste0("_decay_", 
                paste0(substr(in_put, 1,2), collapse = "_"))
  doParallel::registerDoParallel(12)
  model_decay = fit_xgb_cv(c(input_cols, in_put, "slope_N0"),
                           decay_total_sub%>%
                             replace_na(list(drain_level = 0,
                                             drain_conductivity_zeros = 0)),
                           slope_N0~.)
  doParallel::stopImplicitCluster()
  
  readr::write_rds(model_decay, 
                   paste0("compare_inputs/xgboost", cond, ".rds"))
  
}

