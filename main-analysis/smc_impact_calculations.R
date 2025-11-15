##############################
##  Environment and Packages ##
##############################
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)       # for rollmean
library(devtools)

load_all(path = "C:/Users/putne/OneDrive/Documents/git-repos/malclimsim")
####################################
##  File Paths and Core Inputs    ##
####################################
path_to_repo <- "C:/Users/putne/OneDrive/Documents/git-repos/moissala-smc/"
path_to_chad <- "C:/Users/putne/switchdrive/Chad/"

# Source auxiliary functions
source(paste0(path_to_repo, "paper-code/aux-functions/aux_functions.R"))

start_date <- ymd("2018-01-01")
end_date   <- ymd("2023-12-31")

#############################################
##  Load Simulation Results & Observations ##
#############################################
sim_results <- readRDS(
  file.path(path_to_chad, "model-outputs/sim_results.rds")
)

sim_results <- readRDS(
  file.path(path_to_chad, "model-outputs/weighting-noSMC/sim_results.rds")
)

#############################################
##  Defining Necessary Functions ##
#############################################
# Keep original filtering logic
filter_dates_2019 <- function(mat) {
  mat$date_ymd <- as.Date(mat$date_ymd)
  mat[year(mat$date_ymd) == 2019 & month(mat$date_ymd) %in% 6:12, ]
}
filter_dates_2122 <- function(mat) {
  mat$date_ymd <- as.Date(mat$date_ymd)
  mat[year(mat$date_ymd) %in% c(2021, 2022) & month(mat$date_ymd) %in% 6:12, ]
}
filter_dates_2023 <- function(mat) {
  mat$date_ymd <- as.Date(mat$date_ymd)
  mat[year(mat$date_ymd) == 2023 & month(mat$date_ymd) %in% 6:12, ]
}

#############################################
##  Percent Case Change + Cases Averted ##
#############################################

###################
# 2019 Comparison #
###################
true_smc_2019 <- lapply(sim_results$true_SMC, filter_dates_2019)
smc_2019      <- lapply(sim_results$SMC_in_2019, filter_dates_2019)

true_2019_total <- sapply(true_smc_2019, function(x) sum(x$inc_C_transformed))
smc_2019_total <- sapply(smc_2019, function(x) sum(x$inc_C_transformed))

# Total cases each strategy
quantile(true_2019_total , c(0.025, 0.5, 0.975))  %>% round()
quantile(smc_2019_total , c(0.025, 0.5, 0.975))  %>% round()

# By stopping SMC in 2019, we have about ...% more cases than if SMC had been kept...
rel_2019 <- quantile((true_2019_total - smc_2019_total) / smc_2019_total , c(0.025, 0.5, 0.975)) %>% round(2)

# By stopping SMC in 2019, we have about ... more cases than if SMC had been kept...
abs_2019 <- quantile(true_2019_total - smc_2019_total, c(0.025, 0.5, 0.975)) %>% round()

rel_2019
abs_2019
###########################
# 2021–2022 5R vs 4R July #
###########################
true_smc_2122   <- lapply(sim_results$true_SMC, filter_dates_2122)
four_rounds_2122 <- lapply(sim_results$SMC_2021_2022, filter_dates_2122)

true_2122_total <- sapply(true_smc_2122, function(x) sum(x$inc_C_transformed))
four_round_july_total <- sapply(four_rounds_2122, function(x) sum(x$inc_C_transformed))

# Total cases each strategy
quantile((true_2122_total / 2), c(0.025, 0.5, 0.975)) %>% round()
quantile((four_round_july_total / 2), c(0.025, 0.5, 0.975)) %>% round()

# By moving to 5 rounds of SMC, we have about ...% more cases than if 4 rounds of SMC had been kept...
rel_2122 <- quantile(((true_2122_total - four_round_july_total) / four_round_july_total), c(0.025, 0.5, 0.975))

# By moving to 5 rounds of SMC, we have about less cases (per year) than if 4 rounds of SMC had been kept...
abs_2122 <- (quantile((true_2122_total - four_round_july_total), c(0.025, 0.5, 0.975)) / 2) %>% round()

rel_2122
abs_2122
#####################
# 2023 June vs July #
#####################
true_smc_2023 <- lapply(sim_results$true_SMC, filter_dates_2023)
july_smc_2023 <- lapply(sim_results$SMC_2023, filter_dates_2023)

true_2023_total <- sapply(true_smc_2023, function(x) sum(x$inc_C_transformed))
july_smc_total <- sapply(july_smc_2023, function(x) sum(x$inc_C_transformed))

# Total cases each strategy
quantile((true_2023_total), c(0.025, 0.5, 0.975)) %>% round()
quantile((july_smc_total), c(0.025, 0.5, 0.975)) %>% round()

# By moving to June start date, we have about ...% more cases than if July start date had been kept...
rel_23 <- quantile(((true_2023_total - july_smc_total) / july_smc_total), c(0.025, 0.5, 0.975))

# By moving to June start date, we have about less cases (per year) than if July start date had been kept...
abs_23 <- (quantile((true_2023_total - july_smc_total), c(0.025, 0.5, 0.975)) / 2) %>% round()

rel_23
abs_23
########################################
# SMC Effectiveness (June to December) #
########################################
# Raw effectiveness
obs_cases <- readRDS(file.path(path_to_chad, "Data/cases-data/MSF-data/opd_district_MOISSALA_2012_2023.rds")) %>%
  filter(date_ymd %in% as.Date((start_date : end_date)))

obs_cases_in_season <- obs_cases %>% filter(month(date_ymd) %in% 6:12)
raw_cases_noSMC <- obs_cases_in_season %>% filter(year(date_ymd) == 2019)
raw_cases_SMC <- obs_cases_in_season %>% filter(year(date_ymd) != 2019)

rr_raw <- (sum(raw_cases_SMC$inc_C) / 5) / sum(raw_cases_noSMC$inc_C)

rr_raw


# Adjusted effectiveness
filter_in_season <- function(df){
  df <- df %>% filter(month(date_ymd) %in% 6:12)
}

noSMC <- lapply(sim_results$no_SMC, filter_in_season)
SMC <- lapply(sim_results$SMC_in_2019, filter_in_season)

no_smc_total <- sapply(noSMC, function(x) sum(x$inc_C_transformed))
smc_total <- sapply(SMC, function(x) sum(x$inc_C_transformed))

rr_adj <- quantile(smc_total / no_smc_total, c(0.025, 0.5, 0.975))
rr_adj %>% round (2)

# Adjusted differences in cases
ca_adj <- -quantile(smc_total - no_smc_total, c(0.025, 0.5, 0.975)) / 6
ca_adj %>% signif(3)

###########################################################
# SMC Effectiveness (June to December, not includng 2019) #
###########################################################
# Adjusted effectiveness
filter_no_2019_in_season <- function(df){
  df <- df %>% filter(month(date_ymd) %in% 6:12,
                      !(year(date_ymd) %in% 2019))
}

noSMC_2 <- lapply(sim_results$no_SMC, filter_no_2019_in_season)
SMC_2 <- lapply(sim_results$SMC_in_2019, filter_no_2019_in_season)

no_smc_total_2 <- sapply(noSMC_2, function(x) sum(x$inc_C_transformed))
smc_total_2 <- sapply(SMC_2, function(x) sum(x$inc_C_transformed))

rr_adj_2 <- quantile(smc_total_2 / no_smc_total_2, c(0.025, 0.5, 0.975))
rr_adj_2 %>% round (2)

# Adjusted differences in cases
ca_adj_2 <- -quantile(smc_total_2 - no_smc_total_2, c(0.025, 0.5, 0.975)) / 5
ca_adj_2 %>% signif(3)

##########################################
# SMC Effectiveness (August to November) #
##########################################
obs_cases_in_season <- obs_cases %>% filter(month(date_ymd) %in% 8:11)
raw_cases_noSMC <- obs_cases_in_season %>% filter(year(date_ymd) == 2019)
raw_cases_SMC <- obs_cases_in_season %>% filter(year(date_ymd) != 2019)

rr_raw <- (sum(raw_cases_SMC$inc_C) / 5) / sum(raw_cases_noSMC$inc_C)

rr_raw

# Adjusted effectiveness
filter_in_season <- function(df){
  df <- df %>% filter(month(date_ymd) %in% 8:11)
}

noSMC <- lapply(sim_results$no_SMC, filter_in_season)
SMC <- lapply(sim_results$SMC_in_2019, filter_in_season)

no_smc_total <- sapply(noSMC, function(x) sum(x$inc_C_transformed))
smc_total <- sapply(SMC, function(x) sum(x$inc_C_transformed))

rr_adj <- quantile(smc_total / no_smc_total, c(0.025, 0.5, 0.975))
rr_adj %>% round (2)

# Adjusted differences in cases
ca_adj <- -quantile(smc_total - no_smc_total, c(0.025, 0.5, 0.975)) / 6
ca_adj %>% signif(3)


###############################
# Reduction in FOI due to SMC #
###############################
# Model parameter effectiveness
mcmc_results <- readRDS(file.path(path_to_chad, "mcmc-results/weighting-noSMC/mcmc_results_final_weighting2019.rds"))
eff_SMC_quants <- quantile(mcmc_results$coda_pars[,"eff_SMC"], c(0.025, 0.5, 0.975)) %>% round (2)


##########################################
## Combining all values in a data frame  #
##########################################
impact_df <- data.frame(Year = c("2019", "2021-22", "2023"),
           Change = c("No SMC", "5th SMC Round", "June SMC Start Date"),
           `RelCasesChange2.5` = c(rel_2019[1], rel_2122[1], rel_23[1]),
           `RelCasesChange50` = c(rel_2019[2], rel_2122[2], rel_23[2]),
           `RelCasesChange97.5` = c(rel_2019[3], rel_2122[3], rel_23[3]),
           `AbsCasesChange2.5` = c(abs_2019[1], abs_2122[1], abs_23[1]),
           `AbsCasesChange50` = c(abs_2019[2], abs_2122[2], abs_23[2]),
           `AbsCasesChange97.5` = c(abs_2019[3], abs_2122[3], abs_23[3]))

impact_df[3:5] <- impact_df[3:5] %>% round(3)
impact_df[6:8] <- impact_df[6:8] %>% round(0)

impact_df
write.csv(impact_df, paste0(path_to_repo, "paper/paper/plots/tables/weighting-noSMC/impact_calculations_weighting_noSMC.csv"), row.names = FALSE)

