##############################
##  Environment and Packages ##
##############################
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)       # for rollmean
library(devtools)

load_all(path = "C:\\Users\\putnni\\Documents\\git-hub-repositories\\malclimsim\\")
####################################
##  File Paths and Core Inputs    ##
####################################
path_to_repo <- "C:/Users/putnni/Documents/git-hub-repositories/moissala-smc/"
path_to_chad    <- "C:/Users/putnni/switchdrive/Chad/"

# Source auxiliary functions
source(paste0(path_to_repo, "paper/aux-functions/aux_functions.R"))

start_date <- ymd("2018-01-01")
end_date   <- ymd("2023-12-31")

#############################################
##  Load Simulation Results & Observations ##
#############################################
sim_comp_1 <- readRDS(
  file.path(path_to_chad, "model-outputs/strategy_comparison_results_4JulyBaseline_1.rds")
)

sim_comp_15 <- readRDS(
  file.path(path_to_chad, "model-outputs/strategy_comparison_results_4JulyBaseline_15.rds")
)

sim_comp_1 <- readRDS(
  file.path(path_to_chad, "model-outputs/weighting-noSMC/strategy_comparison_results_4JulyBaseline_1_weighting_noSMC.rds")
)

sim_comp_15 <- readRDS(
  file.path(path_to_chad, "model-outputs/weighting-noSMC/strategy_comparison_results_4JulyBaseline_15_weighting_noSMC.rds")
)

##################################
## Defining necessary functions ##
##################################
calc_perc_red <- function(i, df_list, base_list) {
  base_cases   <- base_list[[i]]$inc_C %>% sum()
  July4R_cases <- df_list$`4 rounds (July start)`[[i]]$inc_C %>% sum()
  June4R_cases <- df_list$`4 rounds (June start)`[[i]]$inc_C %>% sum()
  July5R_cases <- df_list$`5 rounds (July start)`[[i]]$inc_C %>% sum()
  June5R_cases <- df_list$`5 rounds (June start)`[[i]]$inc_C %>% sum()

  tibble(
    index = i,
    July4R = July4R_cases / base_cases,
    June4R = June4R_cases / base_cases,
    July5R = July5R_cases / base_cases,
    June5R = June5R_cases / base_cases
  )
}

calc_case_diff <- function(i, df_list, base_list) {
  base_cases   <- base_list[[i]]$inc_C %>% sum()
  July4R_cases <- df_list$`4 rounds (July start)`[[i]]$inc_C %>% sum()
  June4R_cases <- df_list$`4 rounds (June start)`[[i]]$inc_C %>% sum()
  July5R_cases <- df_list$`5 rounds (July start)`[[i]]$inc_C %>% sum()
  June5R_cases <- df_list$`5 rounds (June start)`[[i]]$inc_C %>% sum()

  tibble(
    index = i,
    July4R = (July4R_cases - base_cases) / 6,
    June4R = (June4R_cases - base_cases) / 6,
    July5R = (July5R_cases - base_cases) / 6,
    June5R = (June5R_cases - base_cases) / 6
  )
}


#######################################################################
## Cases difference and % different for 15th of month SMC start date ##
#######################################################################
base_list <- sim_comp_15$o1$`4 rounds (July start)`
n_samples <- sim_comp_15$o1$`4 rounds (July start)` %>% length()

df_list_15 <- sim_comp_15$o2

results_perc_15 <- purrr::map_dfr(1:n_samples, ~ calc_perc_red(.x, df_list_15, base_list))
results_abs_15 <- purrr::map_dfr(1:n_samples, ~ calc_case_diff(.x, df_list_15, base_list))

quantile(results_perc_15$July4R, c(0.025, 0.5, 0.975))
june4_15_rel <- quantile(results_perc_15$June4R, c(0.025, 0.5, 0.975))
july5_15_rel <- quantile(results_perc_15$July5R, c(0.025, 0.5, 0.975))
june5_15_rel <- quantile(results_perc_15$June5R, c(0.025, 0.5, 0.975))

quantile(results_abs_15$July4R, c(0.025, 0.5, 0.975))
june4_15_abs <- quantile(results_abs_15$June4R, c(0.025, 0.5, 0.975))
july5_15_abs <- quantile(results_abs_15$July5R, c(0.025, 0.5, 0.975))
june5_15_abs <- quantile(results_abs_15$June5R, c(0.025, 0.5, 0.975))

######################################################################
## Cases difference and % different for 1st of month SMC start date ##
######################################################################
df_list_1 <- sim_comp_1$o2

results_perc_1 <- purrr::map_dfr(1:n_samples, ~ calc_perc_red(.x, df_list_1, base_list))
results_abs_1 <- purrr::map_dfr(1:n_samples, ~ calc_case_diff(.x, df_list_1, base_list))

july4_1_rel <- quantile(results_perc_1$July4R, c(0.025, 0.5, 0.975))
june4_1_rel <- quantile(results_perc_1$June4R, c(0.025, 0.5, 0.975))
july5_1_rel <- quantile(results_perc_1$July5R, c(0.025, 0.5, 0.975))
june5_1_rel <- quantile(results_perc_1$June5R, c(0.025, 0.5, 0.975))

july4_1_abs <- quantile(results_abs_1$July4R, c(0.025, 0.5, 0.975))
june4_1_abs <- quantile(results_abs_1$June4R, c(0.025, 0.5, 0.975))
july5_1_abs <- quantile(results_abs_1$July5R, c(0.025, 0.5, 0.975))
june5_1_abs <- quantile(results_abs_1$June5R, c(0.025, 0.5, 0.975))

######################################################################
## Constructing a table with results ##
######################################################################
make_summary_df <- function(scenario, label, rel, abs) {
  tibble(
    scenario = scenario,
    label = label,
    rel_lower = rel[[1]] %>% round(3),
    rel_median = rel[[2]] %>% round(3),
    rel_upper = rel[[3]] %>% round(3),
    abs_lower = abs[[1]] %>% round(0),
    abs_median = abs[[2]] %>% round(0),
    abs_upper = abs[[3]] %>% round(0)
  )
}

summary_15 <- bind_rows(
  make_summary_df("July 4R", "15th of month", quantile(results_perc_15$July4R, c(0.025, 0.5, 0.975)),
                  quantile(results_abs_15$July4R, c(0.025, 0.5, 0.975))),
  make_summary_df("June 4R", "15th of month", june4_15_rel, june4_15_abs),
  make_summary_df("July 5R", "15th of month", july5_15_rel, july5_15_abs),
  make_summary_df("June 5R", "15th of month", june5_15_rel, june5_15_abs)
)

summary_1 <- bind_rows(
  make_summary_df("July 4R", "1st of month", quantile(results_perc_1$July4R, c(0.025, 0.5, 0.975)),
                  quantile(results_abs_1$July4R, c(0.025, 0.5, 0.975))),
  make_summary_df("June 4R", "1st of month", june4_1_rel, june4_1_abs),
  make_summary_df("July 5R", "1st of month", july5_1_rel, july5_1_abs),
  make_summary_df("June 5R", "1st of month", june5_1_rel, june5_1_abs)
)

# Combine into one
summary_all <- bind_rows(summary_15, summary_1)
summary_all

# Write df to a .csv

readr::write_csv(summary_all, paste0(path_to_repo, "paper/paper/plots/tables/weighting-noSMC/optimal_smc_timing_table_weighting_noSMC.csv"))
readr::write_csv(summary_all, paste0(path_to_repo, "paper/paper/plots/tables/optimal_smc_timing_table.csv"))







