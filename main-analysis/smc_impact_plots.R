###############################
##  Environment and Packages ##
###############################
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(cowplot)
library(zoo)       # for rollmean
library(devtools)
library(patchwork)

#install.packages("odin.dust", repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
#install.packages("mcstate", repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))

detach("package:malclimsim", unload = TRUE)
devtools::load_all(path = "C:/Users/putne/OneDrive/Documents/git-repos/malclimsim")
#library(malclimsim)

# MASTER COLOURS FOR ALL SCENARIOS
SMC_COLORS <- c(
  "SMC in 2019"               = "#c79a70",
  "4 rounds in 2021–2022" = "#0072b2",
  "July start in 2023" = "#009e73",
  "No SMC" = "#333333",
  "SMC" = "#c79a70",
  "Implemented SMC schedule"         = "#d55e00"
)

SMC_LINETYPES <- c(
  "SMC in 2019"               = "dotted",
  "4 rounds in 2021–2022" = "dashed",
  "July start in 2023" = "dotdash",
  "No SMC" = "twodash",
  "SMC" = "dotted",
  "Implemented SMC schedule"         = "solid"
)

####################################
##  File Paths and Core Inputs    ##
####################################
path_to_repo <- "C:/Users/putne/OneDrive/Documents/git-repos/moissala-smc/"
path_to_chad <- "C:/Users/putne/switchdrive/Chad/"
path_to_figures <- file.path(path_to_repo, "paper-plots")

start_date <- ymd("2018-01-01")
end_date   <- ymd("2023-12-31")

source(paste0(path_to_repo, "paper-code/aux-functions/aux_functions.R"))
source(paste0(path_to_repo, "paper-code/aux-functions/plotting_functions.R"))

######################################
## Making SMC Patterns for Plotting ##
######################################
# Define basic monthly SMC patterns used across scenarios
pattern_4_rounds_june <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
pattern_4_rounds_july <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0)
pattern_5_rounds_june <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0)
pattern_5_rounds_july <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0)
pattern_none          <- rep(0, 12)

# Define the actual (historical) SMC timing schedule by year
true_schedule <- list(
  "2018" = pattern_4_rounds_july,
  "2019" = pattern_none,
  "2020" = pattern_4_rounds_july,
  "2021" = pattern_5_rounds_july,
  "2022" = pattern_5_rounds_july,
  "2023" = pattern_5_rounds_june
)

smc_patterns <- list()

# Add the Implemented SMC Schedule to the list of evaluated patterns
smc_patterns[["Implemented SMC schedule"]] <- true_schedule

# Counterfactual 1: use 4 rounds in July for 2021 and 2022 (rest as observed)
cf_4rounds_july_21_22 <- true_schedule
cf_4rounds_july_21_22[["2021"]] <- pattern_4_rounds_july
cf_4rounds_july_21_22[["2022"]] <- pattern_4_rounds_july
smc_patterns[["4 rounds in 2021–2022"]] <- cf_4rounds_july_21_22

# Counterfactual 2: use 5 rounds in July for 2023 (instead of June)
cf_5rounds_july_2023 <- true_schedule
cf_5rounds_july_2023[["2023"]] <- pattern_5_rounds_july
smc_patterns[["July start in 2023"]] <- cf_5rounds_july_2023

smc_in_2019 <- list(
  "2018" = pattern_4_rounds_july,
  "2019" = pattern_4_rounds_july,
  "2020" = pattern_4_rounds_july,
  "2021" = pattern_5_rounds_july,
  "2022" = pattern_5_rounds_july,
  "2023" = pattern_5_rounds_june
)

smc_patterns[["SMC in 2019"]] <- smc_in_2019

# Counterfactual: no SMC the entire time period
smc_no <- list(
  "2018" = pattern_none,
  "2019" = pattern_none,
  "2020" = pattern_none,
  "2021" = pattern_none,
  "2022" = pattern_none,
  "2023" = pattern_none
)

smc_patterns[["No SMC"]] <- smc_no

######################
## Loading SMC data  #
######################
path_to_SMC  <- file.path(path_to_chad, "Data/data-raw/CPS/CPS_coverage_imput_2018_with_2023.xlsx")
smc_data <- load_clean_smc_data(path_to_excel = path_to_SMC)

# Defining SMC start dates for plotting
years <- 2018:2023
smc_day_mat_true <- build_smc_day_matrix(smc_data, years)

# 2019 starting SMC days
smc_day_mat_2019 <- smc_day_mat_true
smc_day_mat_2019["2019", 7:10] <- smc_day_mat_true["2018", 7:10]

# 2021-22 starting SMC days
smc_day_mat_2122 <- smc_day_mat_true
smc_day_mat_2122["2021",11] <- NA

# 2023 starting SMC days
smc_day_mat_23 <- smc_day_mat_true
smc_day_mat_23["2023", 6] <- NA
smc_day_mat_23["2023", 11] <- 7

smc_day_mat_list <- list(
  "Implemented SMC schedule" = smc_day_mat_true,
  "SMC in 2019" = smc_day_mat_2019,
  "4 rounds in 2021–2022" = smc_day_mat_2122,
  "July start in 2023" = smc_day_mat_23
)

############################################
##  Load Simulation Results & Observations ##
############################################
# For results using no weighting
sim_results <- readRDS(
  file.path(path_to_chad, "model-outputs/sim_results.rds")
)

mcmc_results <- readRDS(
  file.path(path_to_chad, "mcmc-results/mcmc_results_final.rds")
)

# For results using likelihood weighting only for no SMC vs SMC years
sim_results <- readRDS(
  file.path(path_to_chad, "model-outputs/weighting-noSMC/sim_results.rds")
)

mcmc_results <- readRDS(
  file.path(path_to_chad, "mcmc-results/weighting-noSMC/mcmc_results_final_weighting2019.rds")
)

obs_cases <- mcmc_results$incidence_df

################################################
##  Summarize Simulation Output for Plotting   ##
################################################
true_summary <- summarize_simulation_ci(sim_results$true_SMC, variables = "inc_C_transformed") %>% mutate(label = "Implemented SMC schedule")
in_2019_summary <- summarize_simulation_ci(sim_results$SMC_in_2019, variables = "inc_C_transformed") %>% mutate(label = "SMC in 2019")
four_round_2021_22_summary <- summarize_simulation_ci(sim_results$SMC_2021_2022, variables = "inc_C_transformed") %>% mutate(label = "4 rounds in 2021–2022")
five_round_july_23_summary <- summarize_simulation_ci(sim_results$SMC_2023, variables = "inc_C_transformed") %>% mutate(label = "July start in 2023")
no_SMC_summary <- summarize_simulation_ci(sim_results$no_SMC, variables = "inc_C_transformed") %>% mutate(label = "No SMC")


# Combine all (excluding full_SMC)
case_summaries <- list(
  `Implemented SMC schedule`      = true_summary,
  `SMC in 2019`   = in_2019_summary,
  `4 rounds in 2021–2022` = four_round_2021_22_summary,
  `July start in 2023` = five_round_july_23_summary,
  `No SMC` = no_SMC_summary
)

###############################################
## What if SMC had been implemented in 2019? ##
###############################################
# ---- Plot Settings ----
x_axis_label_size <- 12.5
y_axis_title_size <- 14.5
#slide_2019_levels <- c("Implemented SMC Schedule", "SMC in 2019")
slide_2019_levels <- c("SMC in 2019", "Implemented SMC schedule")


width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 3.2
height_1yr    <- 5.2

# ---- Full Time Series Plots ----
true_df_full     <- case_summaries[["Implemented SMC schedule"]]
in_2019_df_full  <- case_summaries[["SMC in 2019"]]
compare_df_full  <- bind_rows(true_df_full, in_2019_df_full)

xlim_2019   <- range(true_df_full$date_ymd)
ylim_true   <- c(0, max(true_df_full$upper, na.rm = TRUE) * 1.2)
ylim_both   <- c(0, max(compare_df_full$upper, na.rm = TRUE) * 1.2)

# Plot: True SMC + SMC in 2019
p_2019_1 <- create_time_series_plot(summary_df = compare_df_full,
                        obs_df         = NULL,
                        smc_day_of_month_list = smc_day_mat_list,
                        xlim           = xlim_2019,
                        ylim           = ylim_both,
                        xlab           = NULL,
                        ylab           = "Weekly Cases (<5 yrs)",
                        scenario_levels = slide_2019_levels,
                        scenario_cols = SMC_COLORS[slide_2019_levels],
                        x_axis_label_size = 12.5,
                        y_axis_title_size = 14.5,
                        legend_text_size = 13,
                        axis_text_size = 10,
                        legend_position = "None",
                        add_x_lab = FALSE)


# ---- Zoomed 2019 Plots ----
true_df_2019 <- true_df_full %>%
  filter(year(date_ymd) == 2019, month(date_ymd) %in% 5:12)

in_2019_df_2019 <- in_2019_df_full %>%
  filter(year(date_ymd) == 2019, month(date_ymd) %in% 5:12)

compare_df_2019 <- bind_rows(true_df_2019, in_2019_df_2019)

xlim_2019_zoom   <- range(true_df_2019$date_ymd)
ylim_2019_true   <- c(0, max(true_df_2019$upper, na.rm = TRUE) * 1.2)
ylim_2019_both   <- c(0, max(compare_df_2019$upper, na.rm = TRUE) * 1.2)

p_2019_2 <- create_time_series_plot(summary_df = compare_df_2019,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat_list,
                                    xlim           = xlim_2019_zoom,
                                    ylim           = ylim_2019_true,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = slide_2019_levels,
                                    scenario_cols = SMC_COLORS[slide_2019_levels],
                                    scenario_linetypes = SMC_LINETYPES[slide_2019_levels],
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "none",
                                    add_x_lab = TRUE,
                                    plot_title = "B",
                                    title_size = 14)

####################################################
## What if 4 rounds of SMC had been used in 21-22? #
####################################################
# ---- Plot Settings ----
x_axis_label_size <- 12.5
y_axis_title_size <- 14.5

width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 3.2
height_1yr    <- 5.2

# ---- Scenario Labels ----
#compare_labels <- c("Implemented SMC Schedule", "4 rounds in 2021–2022")
compare_labels <- c("4 rounds in 2021–2022", "Implemented SMC schedule")

# ---- Full Time Series ----
true_df_full <- case_summaries[["Implemented SMC schedule"]]
cf_df_full   <- case_summaries[["4 rounds in 2021–2022"]]

compare_df_full <- bind_rows(true_df_full, cf_df_full)

xlim_full <- range(true_df_full$date_ymd)
ylim_full <- c(0, max(compare_df_full$upper, na.rm = TRUE) * 1.2)

# True SMC + 4 Rounds July
p_2021_1 <- create_time_series_plot(summary_df = compare_df_full,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat_list,
                                    xlim           = xlim_full,
                                    ylim           = ylim_full,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels,
                                    scenario_cols = SMC_COLORS[compare_labels],
                                    scenario_linetypes = SMC_LINETYPES[compare_labels],
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "none",
                                    add_x_lab = FALSE)

# 2021–2022 Only
true_df_2122 <- true_df_full %>%
  filter(year(date_ymd) %in% c(2021, 2022), month(date_ymd) %in% 1:12)

cf_df_2122 <- cf_df_full %>%
  filter(year(date_ymd) %in% c(2021, 2022), month(date_ymd) %in% 1:12)

compare_df_2122 <- bind_rows(true_df_2122, cf_df_2122)

xlim_2122 <- range(true_df_2122$date_ymd)
ylim_2122 <- c(0, max(compare_df_2122$upper, na.rm = TRUE) * 1.2)

p_2021_2 <- create_time_series_plot(summary_df = compare_df_2122,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat_list,
                                    xlim           = xlim_2122,
                                    ylim           = ylim_2122,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels,
                                    scenario_cols = SMC_COLORS[compare_labels],
                                    scenario_linetypes = SMC_LINETYPES[compare_labels],
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "none",
                                    add_x_lab = TRUE,
                                    plot_title = "2021-2022",
                                    title_size = 14)


# 2021 Only (June to December)
true_df_2021 <- true_df_full %>%
  filter(year(date_ymd) == 2021, month(date_ymd) %in% 6:12)

cf_df_2021 <- cf_df_full %>%
  filter(year(date_ymd) == 2021, month(date_ymd) %in% 6:12)

compare_df_2021 <- bind_rows(true_df_2021, cf_df_2021)

xlim_2021 <- range(true_df_2021$date_ymd)
ylim_2021 <- c(0, max(compare_df_2021$upper, na.rm = TRUE) * 1.2)

# 2021 Comparison
p_2021_3 <- create_time_series_plot(summary_df = compare_df_2021,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat_list,
                                    xlim           = xlim_2021,
                                    ylim           = ylim_2021,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels,
                                    scenario_cols = SMC_COLORS[compare_labels],
                                    scenario_linetypes = SMC_LINETYPES[compare_labels],
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "none",
                                    add_x_lab = TRUE,
                                    plot_title = "C",
                                    title_size = 14)

############################################################
## Impact of SMC starting in June not July in 2023 ##
############################################################
# ---- Plot Settings ----
x_axis_label_size <- 12.5
y_axis_title_size <- 14.5

width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 2.45 / 1.2
height_1yr    <- 6.2 / 1.2

# ---- Plot Settings ----
x_axis_label_size <- 12.5
y_axis_title_size <- 14.5

width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 3.2
height_1yr    <- 5.2

# Scenario Labels
# compare_labels <- c("Implemented SMC Schedule", "July Start in 2023")
compare_labels <- c("July start in 2023", "Implemented SMC schedule")

# Full Time Series Data
true_df_full <- case_summaries[["Implemented SMC schedule"]]
cf_df_full   <- case_summaries[["July start in 2023"]]

compare_df_full <- bind_rows(true_df_full, cf_df_full)

xlim_full <- range(compare_df_full$date_ymd)
ylim_full <- c(0, max(compare_df_full$upper, na.rm = TRUE) * 1.2)

# True SMC + 5 rounds July
p_2023_1 <- create_time_series_plot(summary_df = compare_df_full,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat_list,
                                    xlim           = xlim_full,
                                    ylim           = ylim_full,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels,
                                    scenario_cols = SMC_COLORS[compare_labels],
                                    scenario_linetypes = SMC_LINETYPES[compare_labels],
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "none",
                                    add_x_lab = FALSE)


# 2023 Only
true_df_2023 <- true_df_full %>%
  filter(year(date_ymd) == 2023, month(date_ymd) %in% 5:12)

cf_df_2023 <- cf_df_full %>%
  filter(year(date_ymd) == 2023, month(date_ymd) %in% 5:12)

compare_df_2023 <- bind_rows(true_df_2023, cf_df_2023)

xlim_2023 <- range(true_df_2023$date_ymd)
ylim_2023 <- c(0, max(compare_df_2023$upper, na.rm = TRUE) * 1.2)


p_2023_2 <- create_time_series_plot(summary_df = compare_df_2023,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat_list,
                                    xlim           = xlim_2021,
                                    ylim           = ylim_2021,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels,
                                    scenario_cols = SMC_COLORS[compare_labels],
                                    scenario_linetypes = SMC_LINETYPES[compare_labels],
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "none",
                                    add_x_lab = TRUE,
                                    plot_title = "D",
                                    title_size = 14)

#######################################################
## Plotting full time series all strategies one plot ##
#######################################################
# Full Time Series Data
true_summary <- summarize_simulation_ci(sim_results$true_SMC, variables = "inc_C_transformed") %>% mutate(label = "Implemented SMC schedule")
in_2019_summary <- summarize_simulation_ci(sim_results$SMC_in_2019, variables = "inc_C_transformed") %>% mutate(label = "SMC in 2019")
four_round_2021_22_summary <- summarize_simulation_ci(sim_results$SMC_2021_2022, variables = "inc_C_transformed") %>% mutate(label = "4 rounds in 2021–2022")
five_round_july_23_summary <- summarize_simulation_ci(sim_results$SMC_2023, variables = "inc_C_transformed") %>% mutate(label = "July start in 2023")
no_SMC_summary <- summarize_simulation_ci(sim_results$no_SMC, variables = "inc_C_transformed") %>% mutate(label = "No SMC")

compare_df_full <- bind_rows(true_summary, in_2019_summary, four_round_2021_22_summary, five_round_july_23_summary)


compare_labels <- c(
  "SMC in 2019",
  "4 rounds in 2021–2022",
  "July start in 2023",
  "Implemented SMC schedule"
)

legend_order <- c(
  "SMC in 2019",
  "4 rounds in 2021–2022",
  "July start in 2023",
  "Implemented SMC schedule"
)



p_all <- create_time_series_plot(summary_df = compare_df_full,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = NULL,
                                    xlim           = xlim_full,
                                    ylim           = ylim_full,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels,
                                    scenario_cols = SMC_COLORS[compare_labels],
                                 scenario_linetypes =  SMC_LINETYPES[compare_labels],
                                 legend_order = legend_order,
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "right",
                                 plot_title = "A",
                                    add_x_lab = FALSE)

###################################################################
## Plotting full time series all strategies one plot with no SMC ##
###################################################################
compare_df_full_with_no <- bind_rows(true_summary, in_2019_summary, four_round_2021_22_summary, five_round_july_23_summary, no_SMC_summary)

legend_order <- c(
  "Implemented SMC schedule",
  "SMC in 2019",
  "4 rounds in 2021–2022",
  "July start in 2023",
  "No SMC"
)

compare_all_labels <- c(
  "SMC in 2019",
  "4 rounds in 2021–2022",
  "July start in 2023",
  "No SMC",
  "Implemented SMC schedule"
)


p_all_with_no <- create_time_series_plot(summary_df = compare_df_full_with_no,
                                 obs_df         = NULL,
                                 smc_day_of_month_list = NULL,
                                 xlim           = xlim_full,
                                 ylim           = ylim_full,
                                 xlab           = NULL,
                                 ylab           = "Weekly Cases (<5 yrs)",
                                 scenario_levels = compare_all_labels,
                                 scenario_cols = SMC_COLORS[compare_all_labels],
                                 scenario_linetypes =  SMC_LINETYPES[compare_all_labels],
                                 legend_order = legend_order,
                                 x_axis_label_size = 12.5,
                                 y_axis_title_size = 14.5,
                                 axis_text_size = 10,
                                 legend_position = "right",
                                 plot_title = "A",
                                 add_x_lab = FALSE)

###################################################
## Plotting only all SMC vs no SMC entire period ##
###################################################
SMC_COLORS["SMC"] <- SMC_COLORS["SMC in 2019"]
SMC_LINETYPES["SMC"] <- SMC_LINETYPES["SMC in 2019"]

# compare_labels <- c("Implemented SMC Schedule", "July Start in 2023")
compare_labels_eff <- c("SMC", "No SMC")

legend_order_eff <- c(
  "SMC",
  "No SMC"
)

# Full Time Series Data
smc_df_full <- case_summaries[["SMC in 2019"]] %>% mutate(label = "SMC")
no_smc_df_full   <- case_summaries[["No SMC"]]

compare_df_eff <- bind_rows(smc_df_full, no_smc_df_full)

xlim_full <- range(compare_df_full$date_ymd)
ylim_full <- c(0, max(compare_df_full$upper, na.rm = TRUE) * 1.2)

# SMC vs no SMC full time series
p_no_vs_all <- create_time_series_plot(summary_df = compare_df_eff,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = NULL,
                                    xlim           = xlim_full,
                                    ylim           = ylim_full,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = compare_labels_eff,
                                    scenario_cols = SMC_COLORS[compare_labels_eff],
                                    scenario_linetypes = SMC_LINETYPES[compare_labels_eff],
                                    legend_order = legend_order_eff,
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    axis_text_size = 10,
                                    legend_position = "top",
                                    add_x_lab = FALSE)


######################################################
## Plotting implemented SMC vs no SMC entire period ##
#####################################################
# compare_labels <- c("Implemented SMC Schedule", "July Start in 2023")
compare_labels_eff <- c("Implemented SMC schedule", "No SMC")

legend_order_eff <- c("Implemented SMC schedule", "No SMC")

# Full Time Series Data
smc_df_full <- case_summaries[["Implemented SMC schedule"]]
no_smc_df_full   <- case_summaries[["No SMC"]]

compare_df_eff <- bind_rows(smc_df_full, no_smc_df_full)

xlim_full <- range(compare_df_full$date_ymd)
ylim_full <- c(0, max(compare_df_full$upper, na.rm = TRUE) * 1.2)

# SMC vs no SMC full time series
p_no_vs_true <- create_time_series_plot(summary_df = compare_df_eff,
                                       obs_df         = NULL,
                                       smc_day_of_month_list = NULL,
                                       xlim           = xlim_full,
                                       ylim           = ylim_full,
                                       xlab           = NULL,
                                       ylab           = "Weekly Cases (<5 yrs)",
                                       scenario_levels = compare_labels_eff,
                                       scenario_cols = SMC_COLORS[compare_labels_eff],
                                       scenario_linetypes = SMC_LINETYPES[compare_labels_eff],
                                       legend_order = legend_order_eff,
                                       x_axis_label_size = 12.5,
                                       y_axis_title_size = 14.5,
                                       axis_text_size = 10,
                                       legend_position = "top",
                                       add_x_lab = FALSE)

#####################################################################
## Plotting only all SMC vs no SMC entire period + rainfall mirror ##
#####################################################################
rain_path <- paste0(path_to_chad, "data/climate-data/chirps_moissala.rds")
rain_df <- readRDS(rain_path) %>% filter(year(date) %in% 2018:2023)
colnames(rain_df) <- c("dates", "rainfall")

xlim_shared <- range(
  min(compare_df_eff$date_ymd, na.rm = TRUE),
  max(compare_df_eff$date_ymd, na.rm = TRUE)
)

# Recalculate rainfall if needed
monthly_rain_df <- rain_df %>%
  mutate(year = year(dates),
         month = month(dates)) %>%
  group_by(year, month) %>%
  summarise(rainfall = sum(rainfall, na.rm = TRUE), .groups = "drop") %>%
  mutate(dates = as.Date(sprintf("%04d-%02d-01", year, month))) %>%
  select(dates, rainfall) %>%
  arrange(dates)

# Updated rainfall plot with flipped axis (mirror effect)
p_rain <- ggplot(monthly_rain_df, aes(x = dates, y = rainfall)) +
  geom_col(fill = "steelblue", width = 28) +
  scale_y_reverse() +
  labs(y = "Monthly Cumulative Rainfall (mm)", x = NULL) +
  theme_pubr() + scale_x_date(
    breaks      = "6 months",
    date_labels = "%b\n%Y",
    limits      = xlim_shared,
    expand      = expansion(mult = c(0.01, 0.01))
  ) +
  theme(
    axis.title.y = element_text(size = 13),
    axis.text.y  = element_text(size = 12)
  )

# Modify top plot: remove x-axis title and ticks for alignment
p_no_vs_all <- p_no_vs_all  + theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5, r = 5, b = -30, l = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)# tight bottom
  )

# Combine with patchwork
p_smc_no_smc_with_rain <- (p_no_vs_all / p_rain) +
  plot_layout(heights = c(2, 1.3))  # adjust height ratio as needed

p_smc_no_smc_with_rain


#####################################################################
## Plotting true SMC vs no SMC entire period + rainfall mirror ##
#####################################################################
p_true_no_smc <- p_no_vs_true + theme_pubr() + theme(legend.title = element_blank(),
legend.text = element_text(size = 14))# tight bottom

p_no_vs_true_for_rain <- p_no_vs_true  + theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5, r = 5, b = -30, l = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)# tight bottom
  )

# Combine with patchwork
p_true_no_smc_with_rain <- (p_no_vs_true_for_rain / p_rain) +
  plot_layout(heights = c(2, 1.3))  # adjust height ratio as needed

p_true_no_smc_with_rain

#####################
## Saving Plots ##
#####################
# Location for plots to be saved
path_to_figures <- paste0(path_to_repo, "paper-plots/retrospective-impact-plots/weighting-noSMC/")

# SMC effectiveness plot with rain mirror
ggsave(filename = paste0(path_to_figures, "p_smc_no_smc_with_rain.pdf"),
       plot = p_smc_no_smc_with_rain,
       height = 7,
       width = 11,
       dpi = 500)

# SMC effectiveness plot with rain mirror (with true SMC schedule not always SMC)
ggsave(filename = paste0(path_to_figures, "p_true_no_smc.pdf"),
       plot = p_true_no_smc,
       height = 5,
       width = 11,
       dpi = 500)

# SMC effectiveness plot with rain mirror (with true SMC schedule not always SMC)
ggsave(filename = paste0(path_to_figures, "p_true_no_smc_no_rain.pdf"),
       plot = p_,
       height = 7,
       width = 11,
       dpi = 500)

# Combined plot (Full TS + zoomed without no SMC line)
p_all <- p_all + theme_pubr() + theme(legend.title = element_blank(),
                                      legend.text = element_text(size = 14),
                                      legend.key.size = unit(1.2, "cm"),
                                      title = element_text(face = "bold", size = 14))
p_2019_2 <- p_2019_2 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))
p_2021_2 <- p_2021_2 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))
p_2021_3 <- p_2021_3 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))
p_2023_2 <- p_2023_2 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))

p_zoom <- p_2019_2 + p_2021_3 + p_2023_2 + plot_layout(axis_titles = "collect")

p_combined_plot <- p_all / p_zoom + plot_layout(axis_titles = "collect", heights = c(0.4, 1))
p_combined_plot

ggsave(filename = paste0(path_to_figures, "retro_impact_combined_plot.pdf"),
       plot = p_combined_plot,
       height = 10,
       width = 11,
       dpi = 500)

# Combined plot (Full TS + zoomed with no SMC line)
p_all_with_no <- p_all_with_no + theme_pubr() + theme(legend.title = element_blank(),
                                      legend.text = element_text(size = 14),
                                      legend.key.size = unit(1.2, "cm"),
                                      title = element_text(face = "bold", size = 14))
p_2019_2 <- p_2019_2 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))
p_2021_2 <- p_2021_2 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))
p_2021_3 <- p_2021_3 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))
p_2023_2 <- p_2023_2 + theme_pubr() + theme(legend.position = "none", title = element_text(face = "bold", size = 14))

p_zoom <- p_2019_2 + p_2021_3 + p_2023_2 + plot_layout(axis_titles = "collect")

p_combined_plot <- p_all_with_no / p_zoom + plot_layout(axis_titles = "collect", heights = c(0.4, 1))
p_combined_plot

ggsave(filename = paste0(path_to_figures, "retro_impact_combined_plot_with_no_SMC_line.pdf"),
       plot = p_combined_plot,
       height = 10,
       width = 11,
       dpi = 500)

