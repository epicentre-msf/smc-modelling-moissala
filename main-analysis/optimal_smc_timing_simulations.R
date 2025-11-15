##########
## Setup #
##########
devtools::load_all(path = "C:\\Users\\putnni\\Documents\\git-hub-repositories\\malclimsim\\")
library(malclimsim)
path_to_chad <- "C:/Users/putnni/switchdrive/Chad/"
mcmc_results <- readRDS(file.path(path_to_chad, "mcmc-results/mcmc_results_final.rds"))
mcmc_results <- readRDS(file.path(path_to_chad, "mcmc-results/weighting-noSMC/mcmc_results_final_weighting2019.rds"))

malaria_model <- load_model("model_new_R_with_FOI")
start_date <- ymd("2018-01-01")
end_date   <- ymd("2023-12-31")

##########################################
## Defining average monthly SMC coverage #
##########################################
path_to_SMC  <- file.path(path_to_chad, "Data/data-raw/CPS/CPS_coverage_imput_2018_with_2023.xlsx")
smc_data <- readxl::read_excel(path_to_SMC)

avg_cov <- smc_data %>%
  filter(smc_couv_tot > 0) %>%
  summarize(smc_couv_tot = mean(smc_couv_tot)) %>%
  unlist()

month_cov <- rep(avg_cov, 12)
names(month_cov) <- as.character(1:12)

#################################################
##  Defining different potential SMC strategies #
#################################################
n_samples <- 500
param_samples <- sample_mcmc_steps(mcmc_results$coda_pars, n_samples)

# Define basic monthly SMC patterns used across scenarios
pattern_4_rounds_june <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
pattern_4_rounds_july <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0)
pattern_5_rounds_june <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0)
pattern_5_rounds_july <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0)

######################################
## Simulating comparator time series #
######################################
years <- 2018:2023
pattern_matrix <- matrix(rep(pattern_4_rounds_july, length(years)), nrow = length(years), byrow = TRUE)

# Generate SMC schedule for this pattern
schedule_4july <- gen_smc_schedule(
  start_date       = start_date,
  end_date         = end_date,
  years            = years,
  months_30_days   = FALSE,
  months_active    = pattern_matrix,
  coverage         = month_cov,
  smc_day_of_month = 15
)

# Apply the schedule to parameters
param_inputs_4july <- mcmc_results$param_inputs
param_inputs_4july$SMC     <- schedule_4july$SMC
param_inputs_4july$cov_SMC <- schedule_4july$cov
param_inputs_4july$decay   <- schedule_4july$decay

sim_4july <- run_simulations_from_samples(
  model            = malaria_model,
  param_inputs     = param_inputs_4july,
  param_samples    = param_samples,
  start_date       = start_date,
  end_date         = end_date,
  prewarm_years    = 2,
  mu_transform_C   = NULL,
  covariate_matrix = NULL,
  noise            = FALSE,
  month            = FALSE
)

#####################################
## Define strategies to be compared #
#####################################
# Define basic scenario labels with corresponding monthly patterns
patterns <- list(
  "4 rounds (July start)" = pattern_4_rounds_july,
  "4 rounds (June start)" = pattern_4_rounds_june,
  "5 rounds (June start)" = pattern_5_rounds_june,
  "5 rounds (July start)" = pattern_5_rounds_july
)

# Define target years of simulation
target_years <- 2018:2023

# Expand each pattern to apply uniformly across all target years
expanded_patterns <- lapply(patterns, function(pat) {
  setNames(replicate(length(target_years), pat, simplify = FALSE), as.character(target_years))
})

#######################################
## Run strategy comparison simulation #
#######################################
cases_averted_fn <- function(y1, y0) {
  sum(y1$inc_C) - sum(y0$inc_C)
}

# Run all SMC timing scenarios and collect posterior estimates + time series summaries
strategy_comparison_results_15 <- evaluate_multiple_scenarios(
  patterns             = expanded_patterns,
  smc_day_of_month     = 15,
  model                = malaria_model,
  param_inputs         = mcmc_results$param_inputs,
  param_samples        = param_samples,
  start_date           = start_date,
  end_date             = end_date,
  avg_cov              = month_cov,
  years                = target_years,
  mu_transform_C       = NULL,
  o1                   = sim_4july,
  outcome_fn           = cases_averted_fn,
  out_dir              = NULL,
  month                = FALSE,
  apply_decay          = TRUE,
  exclude_years        = NULL,
  use_SMC_as_covariate = FALSE,
  noise                = FALSE
)

# Run all SMC timing scenarios and collect posterior estimates + time series summaries
strategy_comparison_results_1 <- evaluate_multiple_scenarios(
  patterns             = expanded_patterns,
  smc_day_of_month     = 1,
  model                = malaria_model,
  param_inputs         = mcmc_results$param_inputs,
  param_samples        = param_samples,
  start_date           = start_date,
  end_date             = end_date,
  avg_cov              = month_cov,
  years                = target_years,
  mu_transform_C       = NULL,
  o1                   = sim_4july,
  outcome_fn           = cases_averted_fn,
  out_dir              = NULL,
  month                = FALSE,
  apply_decay          = TRUE,
  exclude_years        = NULL,
  use_SMC_as_covariate = FALSE,
  noise                = FALSE
)

sapply(strategy_comparison_results_15$estimates, quantile, c(0.025, 0.5, 0.975)) / 6
sapply(strategy_comparison_results_1$estimates, quantile, c(0.025, 0.5, 0.975)) / 6


# save the timing‐scenario results for plotting and tables later
output_dir <- file.path(path_to_chad, "model-outputs/")

saveRDS(
  strategy_comparison_results_15,
  file = file.path(output_dir, "weighting-noSMC/strategy_comparison_results_4JulyBaseline_15_weighting_noSMC.rds")
)

saveRDS(
  strategy_comparison_results_1,
  file = file.path(output_dir, "weighting-noSMC/strategy_comparison_results_4JulyBaseline_1_weighting_noSMC.rds")
)
