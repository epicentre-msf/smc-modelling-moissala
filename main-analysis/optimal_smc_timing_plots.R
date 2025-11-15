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
#library(malclimsim)
library(patchwork)
library(ggpubr)
library(forcats)
library(ggpattern)
library(lubridate)
devtools::load_all(path = "C:/Users/putne/OneDrive/Documents/git-repos/malclimsim")

# MASTER COLOURS FOR ALL SCENARIOS
SMC_COLORS <- c(
  "4R July 15 start (comparator)"     = "#e9f542",
  "4R June 15 start"     = "#c79a70",
  "5R June 15 start"     = "#0072b2",
  "5R July 15 start"     = "#009e73",
  "4R July 1 start"     = "#d55e00",
  "4R June 1 start"     = "#c79a70",
  "5R June 1 start"     = "#0072b2",
  "5R July 1 start"     = "#009e73"
)

SMC_LINETYPES <- c(
  "4R July 15 start (comparator)"               = "solid",
  "4R June 15 start" = "dotdash",
  "5R June 15 start" = "dashed",
  "5R July 15 start"         = "dotted",
  "4R July 1 start"               = "twodash",
  "4R June 1 start" = "dotdash",
  "5R June 1 start" = "dashed",
  "5R July 1 start"         = "dotted"
)

####################################
##  File Paths and Core Inputs    ##
####################################
path_to_repo <- "C:/Users/putne/OneDrive/Documents/git-repos/moissala-smc/"
path_to_chad <- "C:/Users/putne/switchdrive/Chad/"
path_to_figures <- file.path(path_to_repo, "paper-plots")

start_date <- ymd("2018-01-01")
end_date   <- ymd("2023-12-31")

source(file.path(path_to_repo, "paper-code/aux-functions/aux_functions.R"))
source(file.path(path_to_repo, "paper-code/aux-functions/plotting_functions.R"))

# read in the timing‐scenario results from the simulation step
strategy_comparison_results_1 <- readRDS(
  file.path(path_to_chad, "model-outputs/strategy_comparison_results_4JulyBaseline_1.rds")
)

strategy_comparison_results_1 <- readRDS(
  file.path(path_to_chad, "model-outputs/weighting-noSMC/strategy_comparison_results_4JulyBaseline_1_weighting_noSMC.rds")
)

# read in the timing‐scenario results from the simulation step
strategy_comparison_results_15 <- readRDS(
  file.path(path_to_chad, "model-outputs/strategy_comparison_results_4JulyBaseline_15.rds")
)

strategy_comparison_results_15 <- readRDS(
  file.path(path_to_chad, "model-outputs/weighting-noSMC/strategy_comparison_results_4JulyBaseline_15_weighting_noSMC.rds")
)

######################################
## Making SMC Patterns for Plotting ##
######################################
# Define target years of simulation
target_years <- 2018:2023

# For when SMC starts 15th of month
# Define basic monthly SMC patterns used across scenarios
pattern_4_rounds_july_15 <- c(NA, NA, NA, NA, NA, NA, 15, 15, 15, 15, NA, NA)
pattern_4_rounds_june_15 <- c(NA, NA, NA, NA, NA, 15, 15, 15, 15, NA, NA, NA)
pattern_5_rounds_june_15 <- c(NA, NA, NA, NA, NA, 15, 15, 15, 15, 15, NA, NA)
pattern_5_rounds_july_15 <- c(NA, NA, NA, NA, NA, NA, 15, 15, 15, 15, 15, NA)


smc_day_mat_july4R_15 <- matrix(data = rep(pattern_4_rounds_july_15, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(2018:2023, sprintf("%02d", 1:12)))

smc_day_mat_june4R_15 <- matrix(data = rep(pattern_4_rounds_june_15, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))

smc_day_mat_july5R_15 <- matrix(data = rep(pattern_5_rounds_july_15, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))

smc_day_mat_june5R_15 <- matrix(data = rep(pattern_5_rounds_june_15, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))

# For when SMC starts 1st of month
pattern_4_rounds_july_1 <- c(NA, NA, NA, NA, NA, NA, 1, 1, 1, 1, NA, NA)
pattern_4_rounds_june_1 <- c(NA, NA, NA, NA, NA, 1, 1, 1, 1, NA, NA, NA)
pattern_5_rounds_june_1 <- c(NA, NA, NA, NA, NA, 1, 1, 1, 1, 1, NA, NA)
pattern_5_rounds_july_1 <- c(NA, NA, NA, NA, NA, NA, 1, 1, 1, 1, 1, NA)

smc_day_mat_july4R_1 <- matrix(data = rep(pattern_4_rounds_july_1, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))

 smc_day_mat_june4R_1 <- matrix(data = rep(pattern_4_rounds_june_1, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))

smc_day_mat_july5R_1 <- matrix(data = rep(pattern_5_rounds_july_1, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))

smc_day_mat_june5R_1 <- matrix(data = rep(pattern_5_rounds_july_1, length(target_years)),
                                nrow = length(target_years), byrow = TRUE,
                                dimnames = list(as.character(target_years), sprintf("%02d", 1:12)))


smc_day_mat <- list(
  "4R July 15 start (comparator)" = smc_day_mat_july4R_15,
  "4R June 15 start" = smc_day_mat_june4R_15,
  "5R July 15 start" = smc_day_mat_july5R_15,
  "5R June 15 start" = smc_day_mat_june5R_15,
  "4R July 1 start" = smc_day_mat_july4R_1,
  "4R June 1 start" = smc_day_mat_june4R_1,
  "5R July 1 start" = smc_day_mat_july5R_1,
  "5R June 1 start" = smc_day_mat_june5R_1

)

################################################
##  Summarize Simulation Output for Plotting   ##
################################################
# 15th month
four_july_df_15 <- strategy_comparison_results_15$summaries[["4 rounds (July start)"]] %>%
  mutate(label = "4R July 15 start (comparator)")

four_june_df_15 <- strategy_comparison_results_15$summaries[["4 rounds (June start)"]] %>%
  mutate(label = "4R June 15 start")

five_july_df_15 <- strategy_comparison_results_15$summaries[["5 rounds (July start)"]] %>%
  mutate(label = "5R July 15 start")

five_june_df_15 <- strategy_comparison_results_15$summaries[["5 rounds (June start)"]] %>%
  mutate(label = "5R June 15 start")

# 1st month
four_july_df_1 <- strategy_comparison_results_1$summaries[["4 rounds (July start)"]] %>%
  mutate(label = "4R July 1 start")

four_june_df_1 <- strategy_comparison_results_1$summaries[["4 rounds (June start)"]] %>%
  mutate(label = "4R June 1 start")

five_july_df_1 <- strategy_comparison_results_1$summaries[["5 rounds (July start)"]] %>%
  mutate(label = "5R July 1 start")

five_june_df_1 <- strategy_comparison_results_1$summaries[["5 rounds (June start)"]] %>%
  mutate(label = "5R June 1 start")

#########################################
## Best strategy vs. baseline strategy ##
#########################################
# ---- Plot Settings ----
x_axis_label_size <- 12.5
y_axis_title_size <- 14.5
levels_1 <- c("4R July 15 start (comparator)", "5R June 15 start")

legend_order <- c("4R July 15 start (comparator)", "5R June 15 start")

width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 3.2
height_1yr    <- 5.2

summary_df_1 <- bind_rows(four_july_df_15, five_june_df_15)

xlim_1   <- range(summary_df_1$date_ymd)
ylim_1   <- c(0, max(summary_df_1$upper, na.rm = TRUE) * 1.2)

p_baseline_best <- create_time_series_plot(summary_df = summary_df_1,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = smc_day_mat,
                                    xlim           = xlim_1,
                                    ylim           = ylim_1,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = levels_1,
                                    scenario_cols = SMC_COLORS[levels_1],
                                    scenario_linetypes = SMC_LINETYPES[levels_1],
                                    legend_order = legend_order,
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    legend_text_size = 13,
                                    axis_text_size = 10,
                                    legend_position = "top",
                                    add_x_lab = FALSE,
                                    smc_rect_shape = "bar",
                                    smc_bar_gap = 0,
                                    smc_bar_height = 0.03,
                                    smc_alpha = 1)

p_baseline_best + theme_pubr() + theme(legend.title = element_blank(),
                                       legend.text = element_text(size = 14))

###########################################
## All 15th month start day on same plot ##
###########################################
levels_2 <- c("4R July 15 start (comparator)",
              "4R June 15 start",
              "5R July 15 start",
              "5R June 15 start")

legend_order_2 <- c("4R July 15 start (comparator)",
                    "4R June 15 start",
                    "5R July 15 start",
                    "5R June 15 start")

width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 3.2
height_1yr    <- 5.2

summary_df_2 <- bind_rows(four_july_df_15, four_june_df_15,
                          five_july_df_15, five_june_df_15)

xlim_2   <- range(summary_df_2$date_ymd)
ylim_2   <- c(0, max(summary_df_2$upper, na.rm = TRUE) * 1.2)

p_15_all <- create_time_series_plot(summary_df = summary_df_2,
                                           obs_df         = NULL,
                                           smc_day_of_month_list = NULL,
                                           xlim           = xlim_1,
                                           ylim           = ylim_1,
                                           xlab           = NULL,
                                           ylab           = "Weekly Cases (<5 yrs)",
                                           scenario_levels = levels_2,
                                           scenario_cols = SMC_COLORS[levels_2],
                                           scenario_linetypes = SMC_LINETYPES[levels_2],
                                           legend_order = legend_order_2,
                                           x_axis_label_size = 12.5,
                                           y_axis_title_size = 14.5,
                                           legend_text_size = 13,
                                           axis_text_size = 10,
                                           legend_position = "top",
                                           add_x_lab = FALSE)

p_15_all + theme_pubr() + theme(legend.title = element_blank(),
                                       legend.text = element_text(size = 14))

###########################################
## All 1st month start day on same plot ##
###########################################
levels_3 <- c("4R July 15 start (comparator)",
              "4R July 1 start",
              "4R June 1 start",
              "5R July 1 start",
              "5R June 1 start")

legend_order_3 <- c("4R July 15 start (comparator)",
                    "4R July 1 start",
                    "4R June 1 start",
                    "5R July 1 start",
                    "5R June 1 start")

width_full_ts <- 7.5
height_full_ts <- 4.5
width_1yr     <- 3.2
height_1yr    <- 5.2

summary_df_3 <- bind_rows(four_july_df_15, four_july_df_1, four_june_df_1,
                          five_july_df_1, five_june_df_1)

xlim_3   <- range(summary_df_3$date_ymd)
ylim_3   <- c(0, max(summary_df_3$upper, na.rm = TRUE) * 1.2)

p_1_all <- create_time_series_plot(summary_df = summary_df_3,
                                    obs_df         = NULL,
                                    smc_day_of_month_list = NULL,
                                    xlim           = xlim_1,
                                    ylim           = ylim_1,
                                    xlab           = NULL,
                                    ylab           = "Weekly Cases (<5 yrs)",
                                    scenario_levels = levels_3,
                                    scenario_cols = SMC_COLORS[levels_3],
                                    scenario_linetypes = SMC_LINETYPES[levels_3],
                                    legend_order = legend_order_3,
                                    x_axis_label_size = 12.5,
                                    y_axis_title_size = 14.5,
                                    legend_text_size = 13,
                                    axis_text_size = 10,
                                    legend_position = "top",
                                    add_x_lab = FALSE)

p_1_all + theme_pubr() + theme(legend.title = element_blank(),
                                legend.text = element_text(size = 14))


#############
## Barplot ##
#############
cases_averted_list_15 <- strategy_comparison_results_15$estimates[c(
  "4 rounds (July start)",
  "4 rounds (June start)",
  "5 rounds (July start)",
  "5 rounds (June start)"
)]

names(cases_averted_list_15) <- c("4R July 15 start",
                                     "4R June 15 start",
                                     "5R July 15 start",
                                     "5R June 15 start")

cases_averted_list_1 <- strategy_comparison_results_1$estimates[c(
  "4 rounds (July start)",
  "4 rounds (June start)",
  "5 rounds (July start)",
  "5 rounds (June start)"
)]

names(cases_averted_list_1) <- c("4R July 1 start",
                                 "4R June 1 start",
                                 "5R July 1 start",
                                 "5R June 1 start")

cases_averted_list <- c(cases_averted_list_15, cases_averted_list_1)


# Convert input into data frame
barplot_df <- purrr::imap_dfr(cases_averted_list, function(est, label) {
  if (is.numeric(est) && !is.list(est)) {
    quant <- quantile(est, probs = c(0.025, 0.975), na.rm = TRUE)
    tibble::tibble(
      strategy       = label,
      cases_averted = mean(est, na.rm = TRUE),
      lower         = quant[1],
      upper         = quant[2]
    )
  } else if (is.list(est) && all(c("mean", "lower", "upper") %in% names(est))) {
    tibble::tibble(
      strategy       = label,
      cases_averted = est$mean,
      lower         = est$lower,
      upper         = est$upper
    )
  } else {
    stop("Each element of 'estimates' must be a numeric vector or a list with 'mean', 'lower', and 'upper'.")
  }
})

# Define categories for use in plotting
barplot_df$smc_start_day <- c(rep(15, 4), rep(1,4))
barplot_df$smc_month_rounds <- c("4R July start", "4R June start", "5R July start", "5R June start",
                                 "4R July start", "4R June start", "5R July start", "5R June start")

barplot_df$cases_averted <- barplot_df$cases_averted / 6
barplot_df$lower <- barplot_df$lower / 6
barplot_df$upper <- barplot_df$upper / 6

# Define plotting function inputs
df <- barplot_df
orientation = "horizontal"

fill_colors_barplot <- c(
  "4R July start" = "#e9f542",
  "4R June start"     = "#c79a70",
  "5R June start"     = "#0072b2",
  "5R July start"     = "#009e73"
)


# Make plot
p_barplot <- plot_smc_barplot(df = barplot_df,
                 fill_colors = fill_colors_barplot,
                 orientation = "horizontal",
                 plot_title = "",
                 x_axis_lab = NULL,
                 y_axis_lab = "Yearly Reduction in Cases")

##################
## Saving plots ##
##################
# Location for plots to be saved
path_to_figures <- paste0(path_to_repo, "paper-plots/optimal-strategy-plots/weighting-noSMC/")

# Main analysis figure
p_baseline_best_save <- p_baseline_best + theme_pubr() +
  ggtitle("A") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18),
        axis.title.y = element_text(size = 16, face = "bold"))
p_baseline_best_save

p_barplot_save <- p_barplot + ggtitle("B") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold"))

opt_strategy_main_figure <- (p_baseline_best_save / p_barplot_save) + plot_layout(heights = c(1, 0.8))
#((p_barplot_save + ggtitle("A")) / (p_baseline_best_save + ggtitle("B"))) + plot_layout(heights = c(1, 0.8))

opt_strategy_main_figure

ggsave(filename = paste0(path_to_figures, "opt_strategy_main_plot.pdf"),
       plot = opt_strategy_main_figure,
       height = 10,
       width = 13,
       dpi = 500)

# Supplementary figure
p_15_all_save <- p_15_all + ggtitle("A") + theme_pubr() + theme(
  plot.title = element_text(face = "bold", size = 18),
  legend.text = element_text(size = 16),
  legend.key.size = unit(1.2, "cm"),
  legend.title = element_blank(),
  axis.title.y = element_text(size = 16, face = "bold"))

p_15_all_save

p_1_all_save <- p_1_all + ggtitle("B") + theme_pubr() + theme(
  plot.title = element_text(face = "bold", size = 18),
  legend.text = element_text(size = 16),
  legend.key.size = unit(1.2, "cm"),
  legend.title = element_blank(),
  axis.title.y = element_text(size = 16, face = "bold"))

opt_strategy_supp_figure <- p_15_all_save / p_1_all_save

ggsave(filename = paste0(path_to_figures, "opt_strategy_supp_plot.pdf"),
       plot = opt_strategy_supp_figure,
       height = 10,
       width = 14,
       dpi = 500)

