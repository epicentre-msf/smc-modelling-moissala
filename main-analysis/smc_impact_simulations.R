##########
## Setup #
##########
detach("package:malclimsim", unload = TRUE)
devtools::load_all(path = "C:\\Users\\putnni\\Documents\\git-hub-repositories\\malclimsim\\")
#library(malclimsim)
path_to_chad <- "C:/Users/putnni/switchdrive/Chad/"

start_date <- ymd("2018-01-01")
end_date   <- ymd("2023-12-31")

# Initialize or load sim_results
# This is so progress isn't lost if R crashes (which happens often)
output_dir <- file.path(path_to_chad, "model-outputs/")
sim_results_path <- file.path(output_dir, "weighting-noSMC/sim_results.rds")
sim_results <- if (file.exists(sim_results_path)) {
  readRDS(sim_results_path)
} else {
  list()
}

######################
## Loading SMC data  #
######################
path_to_SMC  <- file.path(path_to_chad, "Data/data-raw/CPS/CPS_coverage_imput_2018_with_2023.xlsx")
smc_data <- load_clean_smc_data(path_to_excel = path_to_SMC)

##############################
## Defining helper functions #
##############################
# Helper to run + save result if missing
run_and_save <- function(name, expr) {
  # Always read latest state from disk before modifying
  sim_results <- if (file.exists(sim_results_path)) {
    readRDS(sim_results_path)
  } else {
    list()
  }

  if (!name %in% names(sim_results)) {
    message("Running simulation: ", name)
    result <- expr  # <- run before assigning
    sim_results[[name]] <- result
    saveRDS(sim_results, sim_results_path)
    message("Saved: ", name)
  } else {
    message("Skipping ", name, " (already exists)")
  }
}

###########################################
##  Load MCMC Results and Observations    #
###########################################
mcmc_results <- readRDS(file.path(path_to_chad, "mcmc-results/mcmc_results_final.rds"))
mcmc_results <- readRDS(file.path(path_to_chad, "mcmc-results/weighting-noSMC/mcmc_results_final_weighting2019.rds"))
malaria_model <- load_model("model_new_R_with_FOI")
obs_cases <- mcmc_results$incidence_df

########################################################
## Sampling parameter sets from posterior distribution #
########################################################
n_samples      <- 500
param_samples  <- sample_mcmc_steps(mcmc_results$coda_pars, n_samples)

#####################################
## Defining population growth input #
#####################################
r_df <- get_population_scaling(
  n             = nrow(obs_cases),
  month         = FALSE,
  growth_rate_C = 1.000071,
  growth_rate_A = 1.000092
)

pop_growth_df <- data.frame(
  date_ymd = obs_cases$date_ymd,
  r_C      = r_df$r_C
)

mu_transform_C <- function(inc_df, param_inputs) {
  inc_df$inc_C * inc_df$r_C
}

#######################################
## Simulating using true SMC schedule #
#######################################
run_and_save("true_SMC", run_simulations_from_samples(
  model = malaria_model,
  param_inputs = mcmc_results$param_inputs,
  param_samples = param_samples,
  start_date = start_date,
  end_date = end_date,
  prewarm_years = 2,
  mu_transform_C = mu_transform_C,
  covariate_matrix = pop_growth_df,
  noise = FALSE,
  month = FALSE
))

######################################################
## Simulating using no SMC across entire time period #
######################################################
param_inputs_noSMC <- mcmc_results$param_inputs
param_inputs_noSMC$eff_SMC <- 0
param_samples[,"eff_SMC"] <- 0

run_and_save("no_SMC", run_simulations_from_samples(
  model = malaria_model,
  param_inputs = param_inputs_noSMC,
  param_samples = param_samples,
  start_date = start_date,
  end_date = end_date,
  prewarm_years = 2,
  mu_transform_C = mu_transform_C,
  covariate_matrix = pop_growth_df,
  noise = FALSE,
  month = FALSE
))


sim_results$no_SMC[[1]]$inc_C %>% plot()
+##############################
## Simulating no SMC in 2019 #
##############################
# Load SMC coverage data for 2019
SMC_clean_2019 <- load_clean_smc_data(path_to_rds = file.path(path_to_chad, "data/smc-inputs/SMC_in_2019.rds"))
SMC_clean_2019 <- SMC_clean_2019 %>%
  mutate(date_start = as.character(as.Date(date_start))) %>%
  select(date_start, coverage)

smc_in_2019 <- smc_schedule_from_data(
  smc_cov        = SMC_clean_2019,
  months_30_days = FALSE,
  years          = 2018:2023
) %>% filter(dates >= start_date, dates <= end_date)

param_inputs_smc2019 <- mcmc_results$param_inputs
param_inputs_smc2019$SMC     <- smc_in_2019$SMC
param_inputs_smc2019$cov_SMC <- smc_in_2019$cov
param_inputs_smc2019$decay   <- smc_in_2019$decay

# Running simulation
run_and_save("SMC_in_2019", run_simulations_from_samples(
  model            = malaria_model,
  param_inputs     = param_inputs_smc2019,
  param_samples    = param_samples,
  start_date       = start_date,
  end_date         = end_date,
  prewarm_years    = 2,
  mu_transform_C   = mu_transform_C,
  covariate_matrix = pop_growth_df,
  noise            = FALSE,
  month            = FALSE
))

################################
## Simulating 4 rounds 2021-22 #
################################
# Load SMC coverage data for July 2021/22 counterfactual
SMC_clean_2021_2022 <- load_clean_smc_data(path_to_rds = file.path(path_to_chad, "data/smc-inputs/SMC_2021_2022_4rounds.rds"))
SMC_clean_2021_2022 <- SMC_clean_2021_2022 %>%
  mutate(date_start = as.character(as.Date(date_start))) %>%
  select(date_start, coverage)

smc_in_2021_2022 <- smc_schedule_from_data(
  smc_cov        = SMC_clean_2021_2022,
  months_30_days = FALSE,
  years          = 2018:2023
)

param_inputs_smc2021_2022 <- mcmc_results$param_inputs
param_inputs_smc2021_2022$SMC     <- smc_in_2021_2022$SMC
param_inputs_smc2021_2022$cov_SMC <- smc_in_2021_2022$cov
param_inputs_smc2021_2022$decay   <- smc_in_2021_2022$decay

# Running simulation
run_and_save("SMC_2021_2022", run_simulations_from_samples(
  model            = malaria_model,
  param_inputs     = param_inputs_smc2021_2022,
  param_samples    = param_samples,
  start_date       = start_date,
  end_date         = end_date,
  prewarm_years    = 2,
  mu_transform_C   = mu_transform_C,
  covariate_matrix = pop_growth_df,
  noise            = FALSE,
  month            = FALSE
))

###############################
## Simulating July start 2023 #
###############################
# Load SMC coverage data for 2023 counterfactual
SMC_clean_2023 <- load_clean_smc_data(path_to_rds = file.path(path_to_chad, "Data/smc-inputs/SMC_2023_July.rds"))
SMC_clean_2023 <- SMC_clean_2023 %>%
  mutate(date_start = as.character(as.Date(date_start))) %>%
  select(date_start, coverage)

smc_in_2023 <- smc_schedule_from_data(
  smc_cov        = SMC_clean_2023,
  months_30_days = FALSE,
  years          = 2018:2023
)

param_inputs_smc2023 <- mcmc_results$param_inputs
param_inputs_smc2023$SMC     <- smc_in_2023$SMC
param_inputs_smc2023$cov_SMC <- smc_in_2023$cov
param_inputs_smc2023$decay   <- smc_in_2023$decay

# Running simulation
run_and_save("SMC_2023", run_simulations_from_samples(
  model            = malaria_model,
  param_inputs     = param_inputs_smc2023,
  param_samples    = param_samples,
  start_date       = start_date,
  end_date         = end_date,
  prewarm_years    = 2,
  mu_transform_C   = mu_transform_C,
  covariate_matrix = pop_growth_df,
  noise            = FALSE,
  month            = FALSE
))

