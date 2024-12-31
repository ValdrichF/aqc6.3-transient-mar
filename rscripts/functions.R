# Script containing all the functions that I'll reuse
library(data.table)
library(dplyr)
library(tidymodels)
library(DALEX)
library(ggplot2)
# library(stringr)
# library(lubridate)
# library(scales)
library(cowplot)
# library(RColorBrewer)
library(doFuture)

# Functions to help calculate the decay coefficients
calc_decay = function(arr, time){
  n_elements = which(arr<1e3)[1]-1
  n_elements = ifelse(is.na(n_elements), length(arr), n_elements) # all >1000
  if(n_elements==1) return(1)
  
  arr = arr[1:n_elements]/arr[1]
  log_arr = log(arr)
  n_days = as.numeric(time[1:n_elements]-time[1])
  mod = MASS::rlm(log_arr~n_days-1, maxit = 50)
  mod$coefficients[1]*-1
}

pred_exponential = function(arr, time){
  decay_rate = calc_decay(arr, time)
  n_days = as.numeric(time-time[1])
  arr[1]*exp(-decay_rate*n_days)
}

# Same as above but return the model instead of just the coefficients
fit_exp_decay_robust = function(.df, y_var = "total"){
  .df = .df[time>=ymd("2012-02-25")]
  arr = unlist(.df[,..y_var])
  time = .df$time
  
  n_elements = which(arr<1e2)[1]-1
  n_elements = ifelse(is.na(n_elements), length(arr), n_elements) # all >1000
  if(n_elements==1) return(1)
  arr = arr[1:n_elements]/arr[1]
  log_arr = log(arr)
  n_days = as.numeric(time[1:n_elements]-time[1])
  mod = MASS::rlm(log_arr~n_days, maxit = 50)
  if (mod$coefficients[1]>0){
    mod = MASS::rlm(log_arr~n_days-1, maxit = 50)
  }
  mod
}

fit_exp_decay_N0 = function(.df, y_var = "total", thresh = 1e2){
  .df = .df[time>=ymd("2012-02-25")]
  arr = unlist(.df[, ..y_var])
  time = .df$time
  
  n_elements = which(arr<thresh)[1]-1
  n_elements = ifelse(is.na(n_elements), length(arr), n_elements) # all >1000
  if(n_elements==1) return(1)
  arr = arr[1:n_elements]/arr[1]
  log_arr = log(arr)
  n_days = as.numeric(time[1:n_elements]-time[1])
  mod = MASS::rlm(log_arr~n_days-1, maxit = 50)
  mod
}

coef_exp = function(.mod){
  if(is.atomic(.mod)) return(c(0,-1))
  coef = c(0,coefficients(.mod))
  return(tail(coef, 2))
}

fitted_exp = function(.df, decay_coefficient, scale_0, y_var = "total") {
  res_len = .df[,.N]
  .df = .df[time>=ymd("2012-02-25")]
  arr0 = unlist(.df[1,..y_var])
  time = .df$time
  n_days = as.numeric(time-time[1])
  
  res = arr0*scale_0*exp(decay_coefficient*-1*n_days)
  res = replace(rep(NA, res_len),
                seq(to = res_len, length = length(res)),
                res)
  return(res)
}

## Main function that trains the model
fit_models = function(data, formula, n_hyperparameters = 500){
  
  wf = workflow()%>%
    add_formula(formula)
  
  set.seed(123)
  data_folds = vfold_cv(data)
  
  rf_spec = rand_forest(
    "regression",
    mtry = tune(),
    trees = tune(),
    min_n = tune()
  )%>%
    set_engine(
      "ranger",
      max.depth = tune()
    )
  rf_grid = grid_latin_hypercube(
    mtry(c(2,9)),
    trees(),
    min_n(),
    max.depth = tree_depth(),
    size = n_hyperparameters/2
  )
  rf_wf = wf%>%
    add_model(rf_spec)
  
  boost_spec = boost_tree(
    mode = "regression",
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    mtry = tune(),
    trees = 1000,
    learn_rate = tune()
  )%>%
    set_engine("xgboost", eval_metric = 'rmse')
  boost_grid = grid_latin_hypercube(
    tree_depth(),
    min_n(),
    loss_reduction(),
    sample_size = sample_prop(),
    mtry(c(2,9)),
    learn_rate(),
    size = n_hyperparameters
  )
  boost_wf = wf%>%
    add_model(boost_spec)
  
  # svm
  svm_spec = svm_rbf(
    mode = "regression",
    cost = tune(),
    rbf_sigma = tune(),
    margin = tune()
  )%>%
    set_engine("kernlab")
  svm_wf = wf%>%
    add_model(svm_spec)
  svm_grid = grid_latin_hypercube(
    cost(),
    rbf_sigma(),
    svm_margin(),
    size = n_hyperparameters
  )
  
  wf_list = list(rf_wf, boost_wf, svm_wf) 
  grid_list = list(rf_grid, boost_grid, svm_grid)
  met = metric_set(
    rmse, rsq, mape
  )
  
  doParallel::registerDoParallel(8)
  tune_res = mapply(
    tune_grid,
    object = wf_list,
    grid = grid_list,
    MoreArgs = list(resamples = data_folds,
                    metrics = met),
    SIMPLIFY = FALSE
  )
  best_fit = mapply(
    function(res, wf){
      hyperparams = select_best(res, metric = "rsq")
      wf%>%
        finalize_workflow(hyperparams)%>%
        fit(data)
    },
    res = tune_res,
    wf = wf_list,
    SIMPLIFY = FALSE
  )
  doParallel::stopImplicitCluster()
  best_fit
}

fit_xgb_cv = function(cols_to_include, df, formula){
  
  dat_compare = df%>%
    select(all_of(cols_to_include))
  set.seed(123)
  dat_split_tmp = dat_compare%>%
    initial_split(prop = 500/720)
  
  dat_training = training(dat_split_tmp)
  dat_testing = testing(dat_split_tmp)
  
  # norm_recipe = recipe(total~., data = dat_training)%>%
  #   step_unknown(new_level = -8)%>%
  #   step_normalize(all_predictors())%>%
  #   prep(training = dat_training)
  # dat_training = bake(norm_recipe, new_data = NULL)
  # dat_testing = bake(norm_recipe, new_data = dat_testing)
  
  set.seed(024)
  data_folds = vfold_cv(dat_training, v = 10, repeats = 5) # watch out!!!! Dramatically higher computation
  
  # hp-tune based on cv performance
  boost_spec = boost_tree(
    mode = "regression",
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    mtry = tune(),
    trees = 1000,
    learn_rate = tune()
  )%>%
    set_engine("xgboost", eval_metric = 'rmse')
  boost_grid = grid_latin_hypercube(
    tree_depth(),
    min_n(),
    loss_reduction(),
    sample_size = sample_prop(),
    mtry(c(2,9)),
    learn_rate(c(-4, -0.5)),
    size = 500
  )
  met = metric_set(
    rmse, rsq
  )
  boost_wf = workflow()%>%
    add_formula(formula)%>%
    add_model(boost_spec)
  tune_res = tune_grid(boost_wf,
                       grid = boost_grid,
                       resamples = data_folds,
                       metrics = met)
  hyperparams = select_best(tune_res, metric = 'rsq')
  
  # Train xgboost
  model = boost_wf%>%
    finalize_workflow(hyperparams)%>%
    fit(dat_training)
  # Add metadata to the model
  model$training = dat_training
  model$testing = dat_testing
  # save the cv performance
  setDT(hyperparams)
  performance = collect_metrics(tune_res)%>%
    as.data.table()
  on_cols = names(hyperparams)
  performance = performance[hyperparams, on = on_cols]
  
  model$cv_performance = performance
  model
}

fit_lightGBM = function(cols_to_include, df){
  require(bonsai)
  dat_compare = df%>%
    select(all_of(cols_to_include))
  set.seed(123)
  dat_split_tmp = dat_compare%>%
    initial_split(prop = 1-120/720)
  
  dat_training = training(dat_split_tmp)
  dat_testing = testing(dat_split_tmp)
  
  # norm_recipe = recipe(total~., data = dat_training)%>%
  #   step_unknown(new_level = -8)%>%
  #   step_normalize(all_predictors())%>%
  #   prep(training = dat_training)
  # dat_training = bake(norm_recipe, new_data = NULL)
  # dat_testing = bake(norm_recipe, new_data = dat_testing)
  
  set.seed(024)
  data_folds = vfold_cv(dat_training)
  
  # hp-tune based on cv performance
  boost_spec = boost_tree(
    mode = "regression",
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    # sample_size = 0.01,
    mtry = tune(),
    trees = 1000,
    learn_rate = tune()
  )%>%
    set_engine("lightgbm", eval_metric = 'rmse')
  boost_grid = grid_latin_hypercube(
    tree_depth(),
    min_n(),
    loss_reduction(),
    mtry(c(2,9)),
    learn_rate(),
    size = 500
  )
  met = metric_set(
    rmse, rsq, mape
  )
  boost_wf = workflow()%>%
    add_formula(total~.)%>%
    add_model(boost_spec)
  tune_res = tune_grid(boost_wf,
                       grid = boost_grid,
                       resamples = data_folds,
                       metrics = met)
  hyperparams = select_best(tune_res, metric = 'rsq')
  # Train xgboost
  model = boost_wf%>%
    finalize_workflow(hyperparams)%>%
    fit(dat_training)
  model$training = dat_training
  model$testing = dat_testing
  model
}

# Corr plot output is a png with WIDE margins
# this function crops the margins and returns the png as a matrix
# plot it using `recolorize::plotImageArray`
crop_blank_matrix = function(f){
  require(dplyr)
  totals_y = apply(f, 1, sum)
  totals_x = apply(f, 2, sum)
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
  
  min_x = max(0, min_x-30)
  max_x = min(max_x+30, ncol(f))
  min_y = max(0, min_y-30)
  max_y = min(max_y+30, nrow(f))
  
  f_crop = f[min_x:max_x, min_y:max_y,]
  f_crop
}

## Some helper functions
# Function to tell DALEX how to predict tidymodels models
custom_predict = function(model, data){
  predict(
    model,
    new_data = data
  )%>%
    unlist(use.names = FALSE)
}

# Estimates a 95% confidence interval of true population mean RMSE and MAE
error_range = function(truth, estimate){
  
  # temporal resolution of 7 days, so I could accept +-1/2 week
  errors = abs(truth - estimate)
  # using t-test to calculate the range for the mean
  errors_RMSE = t.test(errors^2)$conf.int%>%
    as.numeric()%>%
    sqrt()
  errors_MAE = t.test(errors)$conf.int%>%
    as.numeric()
  list(RMSE = errors_RMSE, MAE = errors_MAE)
}

# Function that plots 6 variables vs a y-variable
# Newer, cleaner version below. plot_variables is here only for old scripts
plot_variables = function(data, y_var, variables, color_col, color_breaks = waiver()){
  # apply through each
  plots = lapply(variables, function(c_name){
    data%>%
      # return the plot
      ggplot()+
      geom_point(aes(get(c_name), get(y_var), colour = get(color_col)), alpha = 0.4)+
      scale_color_steps(trans = "log10", breaks = color_breaks)+
      labs(x = c_name,
           y = "Decay rate [- per day]",
           colour = bquote("Recharge ["*m^3/day*"]"))+
      theme_minimal()
  })
  # plot it using plot_grid
  legend_common = get_legend(plots[[1]])
  plot_grid(plots[[1]]+theme(legend.position = "none"),
            plots[[2]]+theme(legend.position = "none",
                             axis.title.y = element_blank()),
            plots[[3]]+theme(legend.position = "none"),
            plots[[4]]+theme(legend.position = "none",
                             axis.title.y = element_blank()),
            plots[[5]]+theme(legend.position = "none"),
            plots[[6]]+theme(legend.position = "none",
                             axis.title.y = element_blank()),
            align = 'v',
            axis = "tb",
            hjust = -1,
            ncol = 2)%>%
    plot_grid(
      legend_common,
      ncol = 2,
      rel_widths = c(1, .3)
    )
}

plot_vars = function(col_names, df = decay_total_sub, y_var = "total", colour_var = "total_recharge", col_breaks = rch_breaks){
  
  rch_range = range(df$total_recharge)%>%
    log10()
  rch_breaks = (10^seq(rch_range[1], rch_range[2], length.out = 7))%>%
    round(digits = -2)
  
  dat = df%>%
    select(all_of(c(y_var, col_names, colour_var)))%>%
    pivot_longer(all_of(col_names))
  nvars = length(unique(dat$name))
  
  dat%>%
    ggplot()+
    geom_point(aes(value, get(y_var), colour = get(colour_var)), size = 0.7)+
    scale_color_steps(trans = "log10", breaks = rch_breaks)+
    facet_wrap(vars(name), scales = "free_x", strip.position = "bottom", ncol = 3)+
    labs(x = "",
         y = "Log Decay coefficient [per day]",
         colour = bquote("Recharge ["*m^3/day*"]"))+
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "bottom",
          legend.key.width = unit(2, 'cm'))+
    tag_facets()
}

# scatter plots of the 6 most important variables based on a dataframe 'vimp'
plot_imp_variables = function(vimp, y_var, data, must_include = "river_conductance"){
  # extract the important variables from vimp
  imp_vars = vimp%>%
    slice_max(importance, n = 6)%>%
    select(input_name)%>%
    unlist(use.names = F)
  # make sure river conductance is plotted
  if (!(must_include %in% imp_vars)){
    imp_vars[6] = must_include
  }
  plot_vars(imp_vars, data, y_var, color_col = "total_recharge", color_breaks = rch_breaks)
}

# Some extra functions to visualize and analyze the models

diagnostic_trans = function(df){
  df[, ":="(y = 10^y,
            y_hat = 10^y_hat)]
  df[, ":="(residuals = y_hat-y,
            abs_residuals = abs(y_hat-y))]
  df
}
custom_predict = function(model, data){
  predict(
    model,
    new_data = data
  )%>%
    unlist(use.names = FALSE)
}

model_profile_custom = function(ex, variables, groups = NULL, type='accumulated'){
  # accumulated profile if no groups
  if(is.null(groups)){
    cp_profiles = model_profile(
      ex, variables, type = type 
    )
    res = cp_profiles$agr_profiles
    setDT(res)
    setnames(res,
             c("_x_", "_yhat_"),
             c(variables, "yhat"))
    res[, type:=type]
    return(res)
  }
  # else...
  # If groups, use ingredients::ceteris_paribus and then group by equal intervaled bins
  cp_profiles = ingredients::ceteris_paribus(
    ex,
    new_observation = ex$data,
    variables = variables
  )
  setDT(cp_profiles)
  setnames(cp_profiles, 
           c(variables, "_yhat_", groups), 
           c('value', "yhat", "grouping_var"))
  
  cp_profiles[,":="(`_label_` =  variables)]
  
  groups_range = quantile(cp_profiles$grouping_var,
                          c(0.01, 0.99))
  groups_interval = abs(diff(groups_range)/4)
  n_digits = 1-floor(log10(groups_interval))
  groups_interval = round(groups_interval, n_digits)
  groups_min = floor(groups_range[1]/groups_interval)*groups_interval
  groups_breaks = 0:4*groups_interval+groups_min
  groups_labels = paste(groups_breaks[1:4], groups_breaks[2:5], sep = " to ")
  
  cp_profiles = cp_profiles[grouping_var>groups_breaks[1] & grouping_var<=groups_breaks[5]]
  cp_profiles = cp_profiles%>%
    mutate(grouping_quantile = case_when(
      grouping_var < groups_breaks[2] ~ groups_labels[1],
      grouping_var < groups_breaks[3] ~ groups_labels[2],
      grouping_var < groups_breaks[4] ~ groups_labels[3],
      TRUE ~ groups_labels[4]
    ))%>%
    mutate(grouping_quantile = factor(grouping_quantile, levels = groups_labels))
  setkey(cp_profiles, value, grouping_quantile)
  res = cp_profiles[, .(yhat = median(yhat)),
                    .(`_label_`, value, grouping_quantile)]
  setnames(res, 
           c("value", "grouping_quantile"),
           c(variables, groups))
  return(res)
}

analyze_model = function(mod, grouping_col = "river_conductance_zeros"){
  explainers = lapply(
    mod,
    function(model, ...){
      explain(
        model,
        label = model$col,
        data = model$testing%>%select(-total),
        y = model$testing$total,
        ...)},
    predict_function = custom_predict,
    verbose = FALSE
  )
  
  # y - yhat plot
  diagnostics = lapply(
    explainers, 
    function(ex){
      model_diagnostics(ex)%>%
        select(y, y_hat, label)
    }
  )%>%
    rbindlist()%>%
    diagnostic_trans()
  n_labels = length(unique(diagnostics$label))
  p1 = diagnostics%>%
    ggplot()+
    geom_point(aes(y, y_hat, color = label), show.legend = F)+
    geom_abline(slope = 1, intercept = 0, linetype = 2)+
    theme_default_dalex()+
    facet_wrap(vars(label))+
    tag_facets(tag_pool = letters[-(1:n_labels)],
               position = list(x = 0.8, y = 0.1))
  # Feature importance
  feat_imp = lapply(explainers, model_parts)
  p2 = do.call(plot, c(feat_imp, max_vars = 10))+
    labs(title = "Feature Importance", subtitle = "")+
    tag_facets(position = list(x = 0.8, y = 0.1))
  # model profile
  m_profile = lapply(explainers, model_profile_custom, groups = grouping_col)%>%
    rbindlist()
  p3 = m_profile%>%
    ggplot()+
    geom_line(aes(value, mean_yhat, colour = grouping_quantile))+
    facet_wrap(vars(`_label_`), scales = 'free_x')+
    labs(x = "", colour = grouping_col)+
    theme_default_dalex()+
    tag_facets(tag_pool = letters[-(1:(2*n_labels))],
               position = list(x = 0.8, y = 0.1))
  
  # Arrange and plot the figures
  p13 = plot_grid(p1, p3, ncol = 1)
  plot_grid(p2, p13, ncol = 2)
}
