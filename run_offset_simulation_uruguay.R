library(offsetsim)
source('initialise_params_uruguay.R')

user_simulation_params = initialise_user_simulation_params()
user_global_params = initialise_user_global_params()
user_feature_params = initialise_user_feature_params()
user_output_params <- initialise_user_output_params()

osim.run(user_global_params, user_simulation_params, user_feature_params, loglevel = 'TRACE')
current_simulation_folder = find_current_run_folder(base_folder = user_global_params$simulation_folder)
osim.output(user_output_params, current_simulation_folder, loglevel = 'TRACE')
