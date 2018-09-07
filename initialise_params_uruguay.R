initialise_user_global_params <- function(){
  
  global_params = list()
  
  global_params$simulation_folder = paste0(path.expand('~'), '/offset_data/uruguay/')
  
  global_params$feature_raster_files = paste0(global_params$simulation_folder, 'uruguay_data/raster_tiff/ecosystem_services/', 
                                              list.files(path = paste0(global_params$simulation_folder, 'uruguay_data/raster_tiff/ecosystem_services/'), 
                                                         all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE,
                                                         include.dirs = FALSE, no.. = FALSE, pattern = '.tif'))
  
  global_params$planning_units_raster = '~/offset_data/uruguay/uruguay_data/raster_tiff/parcelas_uy.tif'
  global_params$number_of_cores = 'all'
  
  # The number of realizations to run
  global_params$realisation_num = 1
  
  global_params$save_simulation_outputs = TRUE
  
  global_params$overwrite_site_characteristics = TRUE
  global_params$store_zeros_as_sparse = TRUE
  global_params$run_from_simulated_data = FALSE
  global_params$save_simulation_outputs = TRUE
  global_params$overwrite_unregulated_probability_list = FALSE
  global_params$overwrite_dev_probability_list = FALSE
  global_params$overwrite_offset_probability_list = FALSE
  global_params$overwrite_management_dynamics = TRUE
  global_params$overwrite_feature_dynamics = TRUE
  global_params$overwrite_condition_classes = TRUE
  global_params$overwrite_site_features = TRUE
  
  return(global_params)
}


# define feature_layers dynamics by logistic curve

logistic_projection <- function(parcel_vals, min_eco_val, max_eco_val, current_dec_rate, time_vec){
  
  t_sh = -1/current_dec_rate * log( ((parcel_vals - min_eco_val)/(max_eco_val - parcel_vals)))
  
  # define logistic curve given logistic parameter set.
  eco_projected = min_eco_val + (max_eco_val - min_eco_val)/(1 + exp(-current_dec_rate*(time_vec - t_sh)))
  
  return(eco_projected)
}


create_dynamics_set <- function(logistic_params_set, condition_class_bounds, time_vec){
  
  dynamics_set = lapply(seq_along(logistic_params_set), 
                        function(i) lapply(seq_along(logistic_params_set[[i]]),
                                           function(j) lapply(seq_along(logistic_params_set[[i]][[j]]),
                                                              function(k) logistic_projection(parcel_vals = logistic_params_set[[i]][[j]][[k]][1], 
                                                                                              min_eco_val = condition_class_bounds[[i]][[j]][1], 
                                                                                              max_eco_val = condition_class_bounds[[i]][[j]][3], 
                                                                                              current_dec_rate = logistic_params_set[[i]][[j]][[k]][2], 
                                                                                              time_vec = time_vec))))
  dynamics_set = lapply(seq_along(logistic_params_set), 
                        function(i) lapply(seq_along(logistic_params_set[[i]]),
                                           function(j) setNames(dynamics_set[[i]][[j]], c('lower_bound', 'best_estimate', 'upper_bound'))))
  
  return(dynamics_set)
}


initialise_user_simulation_params <- function(){ 
  
  simulation_params = list()
  
  # what subset of features to use in the simulation
  simulation_params$features_to_use_in_simulation = 1
  
  simulation_params$transform_params = array(1, length(simulation_params$features_to_use_in_simulation))
  # The total number of layers to use in the offset calcuation (iterating from the start)
  simulation_params$features_to_use_in_offset_calc = simulation_params$features_to_use_in_simulation
  
  simulation_params$features_to_use_in_offset_intervention = simulation_params$features_to_use_in_simulation
  
  simulation_params$use_offset_metric = TRUE
  
  simulation_params$time_steps = 5
  simulation_params$intervention_num = 5
  
  # The maximum number of parcels can be selected to offset a single development
  
  simulation_params$max_offset_parcel_num = 10
  
  # Stops the offset from delivering any further gains once it has acheived the gains required
  simulation_params$limit_offset_restoration = TRUE
  
  # The probability per parcel of it being unregulatedly cleared, every parcel gets set to this number - set to zero to turn off
  simulation_params$unregulated_loss_prob = 0.001
  
  # Exclude parcels with less than this number of pixels.
  simulation_params$min_site_screen_size = 100
  # ignore parcels with size below this number of elements 
  simulation_params$max_site_screen_size_quantile = 0.95

  
  # when the interventions are set to take place, in this case force to occur once per year
  simulation_params$intervention_vec = build_stochastic_intervention(time_steps = simulation_params$time_steps, 
                                                                            intervention_start = 1, 
                                                                            intervention_end = simulation_params$time_steps, 
                                                                            intervention_num = simulation_params$intervention_num, 
                                                                            sd = 1)
  
  #   c('net_gains', 'restoration_gains', 'avoided_condition_decline', 'avoided_loss',
  #     'protected_condition', 'current_condition', 'restored_condition')
  
  simulation_params$offset_action_params = list(c('net_gains', 'restore'))
  
  # This is the equivalent of offset_calc_type for the dev parcel. Options
  # are: 'current_condition' - losses are calcuated relative to the value of
  # the site at the time of the intervention 
  # 'future_condition' - is the do nothing trjectory of the development site.
  simulation_params$dev_calc_type = 'future_condition'    #'future_condition', 'current_condition' 
  
  # Track accumulated credit from previous exchanges (eithger in current or
  # previous time step) and use them to allow developments to proceed if the
  # credit is large enough. FALSE means ignore any exces credit from offset exchanges
  simulation_params$allow_developments_from_credit = TRUE
  
  # How the development parcels are selected options are 'random' or
  # 'weighted'. Note tha weighted requires an additonal weighting layer. If
  # you are running on your own data you need to specify the weights file in
  # intialise_routines.R  (or put the files in simulation_inputs)
  simulation_params$development_selection_type = 'random'  
  
  # The time horizon in which the offset gains need to equal the devlopment impact
  simulation_params$offset_time_horizon = 30
  
  # Include future legal developments in calculating contribution of avoided
  # losses to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_potential_developments_in_offset_calc = list(FALSE)
  
  # Include future unregulated developments in calculating contribution of avoided losses
  # to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_unregulated_loss_in_offset_calc = list(FALSE)
  
  # Include unregulated clearing in the calculating the contribution of avoided
  # losses to the impact of the development. 
  simulation_params$include_unregulated_loss_in_dev_calc = simulation_params$include_unregulated_loss_in_offset_calc
  
  simulation_params$dev_counterfactual_adjustment = 'as_offset'
  # The development impacts is multiplied by this factor (irrespective of how
  # they were caluclated) and the offset impact then needs to match this
  # multiplied development impact
  simulation_params$offset_multiplier = 1
  
  return(simulation_params)
  
}


user_transform_function <- function(pool_vals, transform_params){
  scaled_scores <- lapply(seq_along(pool_vals), function(i) transform_params[i]/sum(transform_params)*100.68*(1 - exp(-5*( pool_vals[[i]]/transform_params[i] )^2.5) ))
  BAM_score <- sqrt(Reduce('+', scaled_scores[1:5]) * Reduce('+', scaled_scores[1:5]))
  return(BAM_score)
}



initialise_user_feature_params <- function(){
  
  feature_params = list()

  feature_params$background_dynamics_type = 'site_scale'
  feature_params$management_dynamics_type = 'site_scale'
  feature_params$scale_features = FALSE
  
  feature_params$site_sample_type = 'trunc_norm'
  feature_params$initial_site_sd = 0.05
  
  feature_params$initial_site_mean_sd = 0.2
  feature_params$dynamics_sample_type = 'by_initial_value' #'by_initial_value' 
  # Sample the restoration rates from a uniform distribution to they vary per parcel and per feature
  feature_params$management_dynamics_sample_type = 'by_distribution'
  
  feature_params$project_by_mean = TRUE
  
  feature_params$management_update_dynamics_by_differential = TRUE
  feature_params$background_update_dynamics_by_differential = TRUE
  
  feature_params$perform_management_dynamics_time_shift = FALSE
  feature_params$perform_background_dynamics_time_shift = FALSE
  
  feature_params$update_offset_dynamics_by_time_shift = TRUE
  
  feature_params$sample_management_dynamics = TRUE
  
  # Sample the background dynamics from a uniform distribution to they vary per site and per feature
  feature_params$sample_background_dynamics = TRUE
  
  feature_params$simulated_time_vec = 0:80
  
  feature_params$condition_class_bounds = rep(list(list(c(0, 0.5, 1))), 6)
  
  mean_decline_rate = -0.02
  mean_restoration_rate = 0.04
  
  background_logistic_params_set = rep(list(list(list(c(0, mean_decline_rate), c(0.5, mean_decline_rate), c(1, mean_decline_rate)))), 6)
  
  management_logistic_params_set = rep(list(list(list(c(0.01, 0.04), c(0.01, 0.05), c(0.01, 0.06)))), 6)
  
  feature_params$simulated_time_vec = 0:200
  
  feature_params$background_dynamics_bounds <- create_dynamics_set(background_logistic_params_set, 
                                                                   feature_params$condition_class_bounds,
                                                                   feature_params$simulated_time_vec)
  
  
  feature_params$management_dynamics_bounds <- create_dynamics_set(management_logistic_params_set, 
                                                                   feature_params$condition_class_bounds,
                                                                   feature_params$simulated_time_vec)
  
  feature_params$initial_condition_class_bounds = lapply(seq_along(feature_params$background_dynamics_bounds), 
                                                         function(i) lapply(seq_along(feature_params$background_dynamics_bounds[[i]]), 
                                                                            function(j) c(feature_params$background_dynamics_bounds[[i]][[j]]$lower_bound[1], 
                                                                                          feature_params$background_dynamics_bounds[[i]][[j]]$best_estimate[1], 
                                                                                          feature_params$background_dynamics_bounds[[i]][[j]]$upper_bound[1])))
  
  return(feature_params)
}

setup_sub_plots <- function(nx, ny, x_space, y_space){
  par(mfrow = c(ny, nx))
  par(cex = 0.6)
  par(mar = c(x_space, y_space, 1, 0), oma = c(2, 4, 2.5, 0.5))
  
  par(tcl = -0.25)
  par(mgp = c(2, 0.3, 0))
  
}


initialise_user_output_params <- function(){
  output_params = list()
  output_params$output_plot_folder = vector()
  output_params$plot_type = 'impacts' # can be 'outcomes'  or 'impacts' or 'none'
  output_params$realisation_num = 'all' # 'all' or number to plot
  output_params$write_pdf = TRUE
  
  output_params$plot_site = TRUE
  output_params$plot_program = TRUE
  output_params$plot_landscape = TRUE
  output_params$plot_offset_metric = TRUE
  
  output_params$scenario_vec = 'all' #c(1,4,7,10, 8, 2,3,5,6,9,11,12 ) #1:12

  output_params$plot_subset_type = 'all' #c('offset_action_type') # 'offset_calc_type', 'offset_action_type', offset_time_horizon'
  output_params$plot_subset_param = 'all' #c('maintain') # 'net_gains', 'restore', 15
  output_params$features_to_plot = 1:6
  output_params$print_dev_offset_sites = FALSE
  output_params$sets_to_plot = 10
  output_params$site_outcome_plot_lims_set = list(c(0, 1e2))
  output_params$program_outcome_plot_lims_set = list(c(0e6, 1e4))
  output_params$landscape_outcome_plot_lims_set = list(c(0, 2e4))
  output_params$nx = 3 
  output_params$ny = 4
  output_params$site_impact_plot_lims_set = list(c(-1e2, 1e2), c(-1e2, 1e2), c(-1e2, 1e2), c(-1e3, 1e3), c(-1e3, 1e3), c(-1e3, 1e3))
  output_params$program_impact_plot_lims_set = list(c(-1e4, 1e4), c(-2e4, 2e4), c(-2e4, 2e4), c(-2e4, 2e4), c(-2e4, 2e4), c(-2e4, 2e4)) 
  output_params$landscape_impact_plot_lims_set = list(c(-1e4, 1e4), c(-1e4, 1e4), c(-1e4, 1e4), c(-1e5, 1e5), c(-1e5, 1e5), c(-1e5, 1e5))
  
  return(output_params)
}

