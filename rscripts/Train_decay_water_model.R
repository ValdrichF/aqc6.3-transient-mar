### Train the model to predict the decay rate 
source("load_site_data.R")

input_cols = c("recharge", "area_site", "transmissivity", "depth",
               "resistance", "river_conductance", "drain_conductivity_zeros",
               "river_density") 

# fit the model
additional_input = list(character(), "river_stage", "sto_baseline", "sto_unsat", "sto_por")

in_df = list(
  out_dir_name = c("models_response_decay3", "models_response_ini3",
                   "models_water_decay3", "models_water_ini3"),
  formulae = c(slope_total~., initial_total~.,
               slope_water~., initial_water~.)
)

for (i in 1:length(in_df$out_dir_name)) {
  
  if (dir.exists(in_df$out_dir_name[i])) next
  
  dir.create(in_df$out_dir_name[i])
  
  c("The inputs to this model are:",
    c(input_cols, y_var),
    "",
    "Notice that the river density and the river conductance are difference")%>%
    write_lines(paste0(in_df$out_dir_name[i], "/info.txt"))
  
  for (in_put in additional_input[c(1,2,3)]) {
    
    cond_prefix = in_df$out_dir_name[i]%>%
      str_remove("^models")
    
    cond = in_put%>%
      str_remove("_")%>%
      paste0(collapse = "_")%>%
      paste0(cond_prefix, "_", .)
    
    paste0("working on: ", cond)%>%
      print()
    
    y_var = in_df$formulae[[i]][[2]]%>%
      as.character()

    # Load the appropriate site properties, 
    # ie. winter for initial and spring for decay
    csv_path = y_var%>%
        str_detect("initial")%>%
        ifelse(
            "site_properties_winter_gb_500.csv",
            "site_properties_spring_gb_500.csv"
        )
    train_data = load_site_data(csv_path)
    
    doParallel::registerDoParallel(11)
    model_decay = fit_xgb_cv(c(input_cols, in_put, y_var),
                             train_data%>%
                               replace_na(list(drain_level = 0,
                                               drain_conductivity_zeros = 0)),
                             in_df$formulae[[i]])
    doParallel::stopImplicitCluster()
    
    readr::write_rds(model_decay, 
                     paste0(in_df$out_dir_name[i], "/xgboost", cond, ".rds"))
    
  }
  
}



