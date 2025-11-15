##########################
## Inference diagnostics #
##########################
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(scales)
library(malclimsim)
library(mcstate)
devtools::load_all(path = "C:/Users/putne/OneDrive/Documents/git-repos/malclimsim")

# Define file paths
mcmc_results_out_dir = paste0(path_to_chad, "mcmc-results/")
path_to_chad <- "C:/Users/putne/switchdrive/Chad/" # Replace this with location of chad switchdrive folder
path_to_repo <- "C:/Users/putne/OneDrive/Documents/git-repos/moissala-smc/"

# Load mcmc results saved from "inferring_model_parameters.R"
mcmc_results_final <- readRDS(paste0(mcmc_results_out_dir, "weighting-noSMC/mcmc_results_final_weighting2019.rds"))

malaria_model <- load_model("model_new_R_with_FOI")  # Load the deterministic climate model

# Define paramerers
params_to_estimate <- setdiff(colnames(mcmc_results_final$coda_pars),
                              c("log_prior", "log_likelihood", "log_posterior"))

# Observed data
obs_cases <- mcmc_results_final$incidence_df

# Start and end date
start_date <- ymd("2018-01-01")  # Start date of the analysis
end_date <- ymd("2023-12-31")  # End date of the analysis

# Some helpful plotting functions
source(paste0(path_to_repo, "paper-code/aux-functions/plotting_functions.R"))

# Location for plots to be saved
path_to_figures <- paste0(path_to_repo, "paper-plots/mcmc-diagnostic-plots/")
output_dir <- path_to_figures

# Examining trace plots
mcmc_trace <- save_mcmc_trace(mcmc_results_final,
                              output_dir = path_to_figures,
                              file_prefix = "mcmc",
                              save = FALSE,
                              save_format = "pdf",
                              param_exclude = c("log_likelihood", "log_prior", "log_posterior"),
                              trace_colors = NULL,
                              plot_width = 10,
                              plot_height = 2,
                              base_size = 12,
                              line_size = 0.8,
                              line_alpha = 0.5,
                              ncol_traces = 3,
                              param_order = params_to_estimate,
                              ggtheme = ggplot2::theme_classic(base_size = 12))

# Gelman-rubin test
gelman_diag <- MCMC_diag(mcmc_results_final, params = "gelman", save = TRUE,
                         output_dir = path_to_figures)

corr_plot <- save_corr_plot(
  results = mcmc_results_final,
  params        = params_to_estimate,
  sample_size   = 5000,
  corr_method   = "pearson",
  smooth_method = "loess",
  hist_fill     = "grey80",
  point_alpha   = 0.3,
  smooth_size   = 0.5,
  cor_size      = 3.5,
  strip_size    = 4,
  ggtheme       = theme_bw(base_size = 12),
  save          = FALSE,
  output_dir    = path_to_figures,
  file_prefix   = "mcmc_corr",
  save_format   = "pdf",
  width         = 11,
  height        = 11,
  dpi           = 300)

ggsave(paste0(path_to_figures, "mcmc_corr.pdf"), height = 11, width = 11,
       dpi = 300, plot = corr_plot)

# Marginal posterior distributions
post_plots <- save_posterior_prior_plot(
  results_list = list(mcmc_results_final),
  params = params_to_estimate,
  nrow = 3,
  ncol = 3,
  show_true = FALSE,
  true_vals = NULL,
  show_prior = TRUE,
  prior_n = 1000,
  plot_type = "density",
  title = NULL,
  fill_alpha = 0.5,
  line_size = 1,
  base_size = 12,
  show_legend = FALSE,
  prior_fill = "grey50",
  prior_alpha = 0.2,
  palette = RColorBrewer::brewer.pal(8, "Set2"),
  #save_path = paste0(path_to_figures, "marginal_posteriors_plot.pdf"),
  save_path = NULL,
  save_width = 10,
  save_height = 6,
  save_dpi = 1000
)

########################
## Assessing model fit #
########################
# Sampling from the posterior distribution
n_samples <- 500
param_samples <- sample_mcmc_steps(mcmc_results_final$coda_pars, n_samples)

# Population growth data
r_daily_df <- readRDS(paste0(path_to_chad, "data/model-inputs/daily_growth_rates_chad.rds"))

# Create covariate matrix (population growth)
r_df <- get_population_scaling(
  n             = nrow(obs_cases),
  month         = FALSE,
  growth_rate_C = r_daily_df$r_daily_u5,
)

pop_growth_df <- data.frame(
  date_ymd = obs_cases$date_ymd,
  r_C = r_df$r_C
)

# Define transformation for SMC-adjusted incidence
mu_transform_C <- function(inc_df, param_inputs) {
  inc_df$inc_C * inc_df$r_C
}

# Simulating cases using the best parameters
max_posterior_params <- extract_max_posterior_params(mcmc_results_final)
param_inputs_best <- update_param_list(mcmc_results_final$param_inputs, max_posterior_params)

sim_best <- data_sim(
  model         = malaria_model,
  param_inputs  = param_inputs_best,
  start_date    = start_date,
  end_date      = end_date,
  noise         = FALSE,
  mu_transform_C = mu_transform_C,
  mu_transform_A = NULL,
  covariate_matrix = pop_growth_df,
  month = FALSE,
  save = FALSE
)

sim_best %>% group_by(year(date_ymd)) %>%
  summarize(inc_C = sum(inc_C))

obs_cases %>% group_by(year(date_ymd)) %>%
  summarize(inc_C = sum(inc_C))

sim_best_long <- tidyr::pivot_longer(sim_best, -date_ymd, names_to = "variable", values_to = "value")

# Constructing the uncertainty (prediction) interval
simulations <- run_simulations_from_samples(
  model         = malaria_model,
  param_inputs  = param_inputs_best,
  param_samples = param_samples,
  start_date    = start_date,
  end_date      = end_date,
  prewarm_years = 2,
  mu_transform_C = mu_transform_C,
  mu_transform_A = NULL,
  covariate_matrix = pop_growth_df,
  noise         = TRUE,
  month         = FALSE
)

# Summarize the simulations
ci_data <- summarize_simulation_ci(simulations, variables = "inc_C_transformed", ci_level = 0.95)

# Format the simulated and observed data for plotting
ppc_data <- malclimsim::prepare_ppc_data(
  ci_data = ci_data,
  best_df = sim_best_long,
  obs_df  = mcmc_results_final$incidence_df,
  sim_var = "inc_C_transformed",
  obs_var = "inc_C"
)

# Plot the posterior predictive check
path_to_SMC  <- file.path(path_to_chad, "Data/data-raw/CPS/CPS_coverage_imput_2018_with_2023.xlsx")
smc_data <- load_clean_smc_data(path_to_excel = path_to_SMC)

# Defining SMC start dates for plotting
source(paste0(path_to_repo, "paper-code/aux-functions/aux_functions.R"))
years <- 2018:2023
smc_day_mat_true <- build_smc_day_matrix(smc_data, years)

xlim_shared <- range(
  min(ppc_data$date_ymd, na.rm = TRUE),
  max(ppc_data$date_ymd, na.rm = TRUE)
)

ppc_plot <- plot_ppc(ppc_data, ci_data,
                     smc_day_of_month_list = smc_day_mat_true,        # names must match ppc_data$variable
                     smc_rect_shape = "bar",
                     smc_bar_height = 0.03,
                     smc_bar_gap    = 0.012,
                     smc_fill       = "grey10",
                     smc_alpha = 1,
                     title       = "",
                     xlab        = "",
                     ylab        = "Weekly Cases (<5 yrs)",
                     show_obs    = TRUE,
                     show_best   = TRUE,
                     show_ribbon = TRUE
) + theme(legend.title = element_blank(),
          legend.text = element_text(size = 14)) + scale_x_date(
            breaks      = "6 months",
            date_labels = "%b\n%Y",
            limits      = xlim_shared,
            expand      = expansion(mult = c(0.01, 0.01))
          )

#############################
## Model Fit + Climate Plot #
#############################
## Config / Inputs
years_keep   <- 2018:2023
rain_path    <- paste0(path_to_chad, "data/climate-data/chirps_moissala.rds")

# xlim_shared <- c(as.Date("2018-01-01"), as.Date("2023-12-31"))

## Load & Prepare Climate Data (Monthly)
# Temperature (monthly mean)
temp_monthly <- met_365 %>%
  transmute(dates = as.Date(dates), temp) %>%
  filter(year(dates) %in% years_keep) %>%
  mutate(year = year(dates), month = month(dates)) %>%
  group_by(year, month) %>%
  summarise(dates = as.Date(sprintf("%04d-%02d-01", first(year), first(month))),
            temp  = mean(temp, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(dates)

# Rainfall (monthly sum)
monthly_rain <- readRDS(rain_path) %>%
  transmute(dates = as.Date(date), rainfall = rainfall) %>%
  filter(year(dates) %in% years_keep) %>%
  mutate(year = year(dates), month = month(dates)) %>%
  group_by(year, month) %>%
  summarise(dates    = as.Date(sprintf("%04d-%02d-01", first(year), first(month))),
            rainfall = sum(rainfall, na.rm = TRUE),
            .groups  = "drop") %>%
  arrange(dates)

# Join to shared monthly grid
clim_monthly <- monthly_rain %>%
  left_join(temp_monthly, by = "dates")

# x-axis limits (fallback if not provided)
if (!exists("xlim_shared")) {
  xlim_shared <- range(clim_monthly$dates, na.rm = TRUE)
}

## Axis Mapping (Rainfall ↔ Temperature)
r_min <- 0
r_max <- max(clim_monthly$rainfall, na.rm = TRUE)
t_min <- min(clim_monthly$temp,     na.rm = TRUE)
t_max <- max(clim_monthly$temp,     na.rm = TRUE)

temp_to_rain <- function(t) {
  if (isTRUE(t_max == t_min)) return(rep((r_max + r_min)/2, length(t)))
  (t - t_min) / (t_max - t_min) * (r_max - r_min) + r_min
}
rain_to_temp <- function(r) {
  if (isTRUE(r_max == r_min)) return(rep((t_max + t_min)/2, length(r)))
  (r - r_min) / (r_max - r_min) * (t_max - t_min) + t_min
}

clim_monthly <- clim_monthly %>%
  mutate(temp_scaled = temp_to_rain(temp))

## Figure: Rainfall (bars, left) + Temperature (red line, right)
p_rain <- ggplot(clim_monthly, aes(x = dates)) +
  geom_col(aes(y = rainfall), fill = "steelblue", width = 28) +
  geom_line(aes(y = temp_scaled), color = "red", linewidth = 1.4) +
  scale_y_continuous(
    name    = "Monthly Cumulative Rainfall (mm)",
    trans   = "reverse",  # mirror effect
    sec.axis = sec_axis(~ rain_to_temp(.), name = "Avg Monthly Temperature (°C)")
  ) +
  scale_x_date(
    breaks      = "6 months",
    date_labels = "%b\n%Y",
    limits      = xlim_shared,
    expand      = expansion(mult = c(0.01, 0.01))
  ) +
  theme_pubr() +
  theme(
    axis.title.y.left  = element_text(size = 13),
    axis.text.y.left   = element_text(size = 12),
    axis.title.y.right = element_text(size = 13),
    axis.text.y.right  = element_text(size = 12)
  )

## Top PPC Plot (cosmetics only)
p_ppc_final <- ppc_plot +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5, r = 5, b = -30, l = 5),
    legend.title = element_blank(),
    legend.text  = element_text(size = 14)
  )

## Combine Panels
p_model_fit_with_rain_smc <- (p_ppc_final / p_rain) +
  plot_layout(heights = c(2, 1.3))

p_model_fit_with_rain_smc

#################
## Saving Plots #
#################
# Correlation plot
ggsave(paste0(path_to_figures, "mcmc_corr_plot.pdf"), height = 11, width = 11,
       dpi = 1000, plot = corr_plot)

# Posteriors + Trace
mcmc_trace <- mcmc_trace +
  plot_annotation(
    title = "A",
    theme = theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.02),
      plot.subtitle = element_blank(),      # in case you had a subtitle slot
      plot.caption  = element_blank()       # ditto for caption
    )
  )

p_post_plots <- post_plots +
  plot_annotation(
    title = "B",
    theme = theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.02),
      plot.subtitle = element_blank(),      # in case you had a subtitle slot
      plot.caption  = element_blank()       # ditto for caption
    )
  )

mcmc_trace_post_row <- plot_grid(mcmc_trace, p_post_plots, nrow = 2)
mcmc_trace_post_col <- plot_grid(mcmc_trace, p_post_plots, ncol = 2)

ggsave(paste0(path_to_figures, "weighting-noSMC/mcmc_trace_post_row_weighting.pdf"), height = 15, width = 12,
       dpi = 600, plot = mcmc_trace_post_row)
ggsave(paste0(path_to_figures, "weighting-noSMC/mcmc_trace_post_col.pdf"), height = 8, width = 12,
       dpi = 600, plot = mcmc_trace_post_col)


ppc_plot <- ppc_plot + theme_pubr() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.2, "cm"),
        axis.title.y = element_text(size = 16, face = "bold"))

ppc_plot
ggsave(paste0(path_to_figures, "ppc_plot_with_weighting_each_strategy.pdf"), height = 5, width = 10,
       dpi = 600, plot = ppc_plot)


# SMC effectiveness plot with rain mirror
ggsave(filename = paste0(path_to_figures, "p_model_fit_with_rain_smc.pdf"),
       plot = p_model_fit_with_rain_smc,
       height = 7,
       width = 11,
       dpi = 500)
