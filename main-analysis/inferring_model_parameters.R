#library(devtools)
devtools::load_all(path = "C:/Users/putne/OneDrive/Documents/git-repos/malclimsim")
#detach("package:malclimsim", unload = TRUE)
#devtools::install_github("https://github.com/SwissTPH/malclimsim", ref = "development")
library(malclimsim)
library(mcstate)

#############################################
## Climate data downloading and processing  #
#############################################
# Years used for the analysis
years_clim <- 2017:2023
start_date <- ymd("2018-01-01")  # Start date of the analysis
end_date <- ymd("2023-12-31")  # End date of the analysis
years_analysis <- year(start_date) : year(end_date)

# Latitude and longitude where climate data (rainfall and temperature) is to be saved
# Rainfall data is from CHIRPS and temperature data is from ERA5
lat <- 8.34
lon <- 17.77

path_to_chad <- "C:/Users/putne/switchdrive/Chad/" # Replace this with location of chad switchdrive folder

# Defining names of files to be saved
rain_file_name <- "Data/climate-data/chirps_moissala.rds"
temp_file_name <- "Data/climate-data/era5_moissala.nc"

# Saving the rainfall data locally based on specified path
#save_climate_data(lon = lon, lat = lat, years = years_clim,
#                  path_to_data = path_to_chad,
#                  rain_file_name = rain_file_name, temp_file_name = temp_file_name)

# Reading in climate data saved by `save_climate_data`
temp_path <- paste0(path_to_chad, "data/climate-data/era5_moissala.nc")
rain_path <- paste0(path_to_chad, "data/climate-data/chirps_moissala.rds")

met_365 <- process_climate_data(lon, lat, years = years_clim, temp_path = temp_path,
                                rain_path = rain_path, months_30_days = FALSE)

#####################################################################
## Seasonal malaria chemoprevention data downloading and processing #
#####################################################################
path_to_SMC <- file.path(path_to_chad,"data/data-raw/CPS/CPS_coverage_imput_2018_with_2023.xlsx")

# Load and clean SMC data
SMC_clean <- load_clean_smc_data(path_to_excel = path_to_SMC)

# Convert to schedule
smc_schedule <- smc_schedule_from_data(
  smc_cov = SMC_clean,
  months_30_days = FALSE,
  years = years_analysis
)

#############################
## Defining population size #
#############################
district_pops <- readxl::read_excel(paste0(path_to_chad, "data/model-inputs/district_population_estimates.xlsx")) # 2023 population estimates for select districts of Mandoul
N <- district_pops %>% filter(District == "Moissala") %>% dplyr::select(Population) %>% unlist() %>% as.numeric()

########################################
## Ensuring alignment of model inputs  #
########################################
# It is necessary that the climate and SMC inputs begin on the same date
# To ensure this, we have a helper function that checks this condition/trims them accordingly
# Furthermore, as we have a lag on climate, we need the climate vectors to have
# data before the date when SMC starts. The following code helps to ensure everything is aligned.

met_365 <- impute_climate_to_end_date(met_365, max(smc_schedule$dates)) # climate data only available until December 2023

# Specify lag
clim_smc_lag <- 180

# Validate continuity, span, and alignment
# validate_smc_climate_alignment(smc_schedule, met_365, clim_smc_lag)

# Trim & align
aligned <- lag_and_trim_smc_climate(smc_schedule, met_365, clim_smc_lag)
smc_schedule    <- aligned$smc
met_365 <- aligned$climate

validate_smc_climate_alignment(smc_schedule, met_365, clim_smc_lag)

########################################
## Define population growth data frame #
########################################
r_daily_df <- readRDS(paste0(path_to_chad, "data/model-inputs/daily_growth_rates_chad.rds"))

###########################
## Defining malaria model #
###########################
malaria_model <- load_model("model_new_R_with_FOI")  # Load the deterministic climate model

###########################################
## Combining model inputs into one object #
###########################################
# Extract rainfall and temperature data
rain <- met_365$anom  # Standardized rolling mean of rainfall
temp <- met_365$temp  # Temperature data

# Extract key SMC schedule information
SMC <- smc_schedule$SMC  # Indicator for days when an SMC round started (1s and 0s)
decay <- smc_schedule$decay  # Efficacy decay of SMC over time
cov <- smc_schedule$cov  # SMC coverage over time

# Define parameter inputs for the malaria model simulation
param_inputs <- list(
  # Rate and Demographic parameters
  mu_TS = 1/30, mu_IR = 1/5, mu_RS = 1/195,
  mu_EI = 1/10, delta_b = 47/(1000*365), delta_d = 47/(1000*365),
  delta_a = 1/(5 * 365), N = N, percAdult = 0.81,

  # Immunity and Reporting Parameters
  pi_s_1 = 0.75, c_s = 0.12,
  qR = 0.24,

  # Immunity and Reporting Parameters
  fT_C = 0.28,

  # Population Scaling factor
  s = 9,

  # Growth rates
  r_C_0 = r_daily_df$r_daily_u5,
  r_A_0 = r_daily_df$r_daily_o5,

  # Time-varying Inputs
  decay = decay,
  SMC = as.numeric(SMC),
  cov_SMC = cov,
  c_R_D = rain, temp = temp,

  # Climate Parameters
  b = -log(0.92),
  R_opt = 0,
  alpha = 1,
  sigma_LT = 4, sigma_RT = 4,
  k1 = 1.5, lag_R = 30, lag_T = 15,

  # Likelihood function parameters
  size_1 = 15,

  # SMC parameters
  eff_SMC = 0.8,

  # Misc
  clim_SMC_lag = clim_smc_lag
)

##########################
## Loading observed data #
##########################
# Reading in observed data
obs_cases <- readRDS(file.path(path_to_chad,"data/cases-data/MSF-data/opd_district_MOISSALA_2018_2023.rds"))
head(obs_cases)

############################
## Weighting observed data #
############################
# SMC was deployed each year except for 2019
# Means maximizing posterior probability globally gives more weight to non-SMC years
# This helps to minimize the influence of imbalanced treatment-control group sizes
obs_cases$treatment <- ifelse(year(obs_cases$date_ymd) != 2019, 1, 0)
control_weight <- mean(obs_cases$treatment) / (1 - mean(obs_cases$treatment))
obs_cases$obs_weight <- ifelse(year(obs_cases$date_ymd) == 2019, control_weight, 1)

##########################################
## Defining model parameters to estimate #
##########################################
# Listing the model parameters to be estimated
params_to_estimate <- c(lag_R = "lag_R", lag_T = "lag_T",
                        sigma_LT = "sigma_LT", sigma_RT = "sigma_RT",
                        k1 = "k1", size_1 = "size_1",
                        eff_SMC = "eff_SMC", s = "s",
                        b = "b")


###########################################################################
## Defining parameters of the (adaptive) Random-walk Metropolis algorithm #
###########################################################################
params_default <- create_mcmc_params(stage = "stage1")
adaptive_params_1 <- params_default$adaptive_params
control_params_1 <- params_default$control_params
control_params_1$n_steps = 5000

# Defining starting values for chains
start_values_1 <- create_start_values(params_to_estimate, control_params_1,
                                      model = malaria_model, param_inputs = param_inputs,
                                      random = TRUE, seed = NULL)
# Proposal matrix
proposal_matrix_1 <- create_proposal_matrix(params_to_estimate = params_to_estimate,
                                            model = malaria_model, param_inputs = param_inputs)

param_names <- rownames(proposal_matrix_1)

####################################################################
## Defining aspects of the inference procedure not related to MCMC #
####################################################################
# Defining if using monthly or weekly data, age groups, population growth, etc
inf_config <- make_obs_config(
  use_monthly = FALSE,
  age_group = "u5",
  include_pop_growth = TRUE,
)

# Small wrapper around inf_run() + diagnostics + saving
run_stage <- function(proposal_matrix, start_values, control_params,
                      adaptive_params, out_dir, out_file, rerun_n = 1000,
                      rerun_random = TRUE, seed = 24) {

  # Infer with all the shared args baked in
  results <- inf_run(model = malaria_model, param_inputs = param_inputs,
                     control_params = control_params,
                     params_to_estimate = params_to_estimate,
                     proposal_matrix = proposal_matrix,
                     adaptive_params = adaptive_params,
                     start_values = start_values,
                     noise = FALSE, seed = seed,
                     month_unequal_days = FALSE,
                     dates = c(start_date, end_date),
                     synthetic = FALSE, incidence_df = obs_cases,
                     save_trajectories = FALSE, rerun_n = rerun_n,
                     rerun_random = rerun_random, param_priors = NULL,
                     n_years_warmup = 3, obs_config = inf_config)

  # Save to disk
  saveRDS(results, file = paste0(out_dir, out_file))

  # Return the results object
  invisible(results)
}


######################
## Running inference #
######################
mcmc_results_out_dir = paste0(path_to_chad, "mcmc-results/")


### Stage 1 ###
mcmc_stage_1 <- run_stage(proposal_matrix_1, start_values_1,
                          control_params = control_params_1,
                          adaptive_params = adaptive_params_1,
                          out_dir = mcmc_results_out_dir,
                          out_file = "mcmc_results_stage_1.rds")

### Stage 2 ###
mcmc_stage_1 <- readRDS(paste0(mcmc_results_out_dir, "mcmc_results_stage_1.rds"))
inf_params_stage_2 <- update_inf_stage(results_obj = mcmc_stage_1,
                                       proposal_matrix = proposal_matrix_1,
                                       param_names = rownames(proposal_matrix_1),
                                       S_prev = 3000, draw_n = ncol(start_values_1),
                                       shrink = 0.95, stage = "stage2", n_steps = 1000)
rm(mcmc_stage_1)
mcmc_stage_2 <- run_stage(proposal_matrix = inf_params_stage_2$proposal_matrix,
                          start_values = inf_params_stage_2$start_values,
                          control_params = inf_params_stage_2$control_params,
                          adaptive_params = inf_params_stage_2$adaptive_params,
                          out_dir = mcmc_results_out_dir,
                          out_file = "mcmc_results_stage_2.rds")

#trace <- MCMC_diag(mcmc_stage_2, params = "trace")
#rm(trace)

### Stage 3 ###
mcmc_stage_2 <- readRDS(paste0(mcmc_results_out_dir, "mcmc_results_stage_2.rds"))
inf_params_stage_3 <- update_inf_stage(results_obj = mcmc_stage_2,
                                       proposal_matrix = inf_params_stage_2$proposal_matrix,
                                       param_names = rownames(proposal_matrix_1),
                                       S_prev = 3000, draw_n = 4,
                                       shrink = 0.8, stage = "noadapt", n_steps = 1000)
inf_params_stage_3$control_params$n_burnin <- 0
rm(mcmc_stage_2)
mcmc_stage_3 <- run_stage(proposal_matrix = inf_params_stage_3$proposal_matrix,
                          start_values = inf_params_stage_3$start_values,
                          control_params = inf_params_stage_3$control_params,
                          adaptive_params = inf_params_stage_3$adaptive_params,
                          out_dir = mcmc_results_out_dir,
                          out_file = "mcmc_results_stage_3.rds")

#trace <- MCMC_diag(mcmc_stage_3, params = "trace")

### Stage 4 ###
mcmc_stage_3 <- readRDS(paste0(mcmc_results_out_dir, "mcmc_results_stage_3.rds"))
inf_params_stage_4 <- update_inf_stage(results_obj = mcmc_stage_3,
                                       proposal_matrix = inf_params_stage_3$proposal_matrix,
                                       param_names = rownames(proposal_matrix_1),
                                       S_prev = 3000, draw_n = 4,
                                       shrink = 0.8, stage = "noadapt2", n_steps = 30000)
rm(mcmc_stage_3)
mcmc_results_final <- run_stage(proposal_matrix = inf_params_stage_4$proposal_matrix,
                          start_values = inf_params_stage_4$start_values,
                          control_params = inf_params_stage_4$control_params,
                          adaptive_params = inf_params_stage_4$adaptive_params,
                          out_dir = mcmc_results_out_dir,
                          out_file = "mcmc_results_final.rds")


saveRDS(mcmc_results_final, paste0(mcmc_results_out_dir, "mcmc_results_final_weighting2019.rds"))
