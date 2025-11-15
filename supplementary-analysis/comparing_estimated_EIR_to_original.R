##############################
##  Load Required Libraries ##
##############################
library(ggplot2)
library(dplyr)
library(tidyr)
library(malclimsim)

#################################################
##  Climate–entomological function definitions ##
#################################################
path_to_repo <- "C:/Users/putnni/Documents/git-hub-repositories/moissala-smc/"
source(paste0(path_to_repo, "paper/supplementary-analysis/original_climate_equation_definitions.R"))

##################################################################
##  Finding EIR-temperature parameters that maximized posterior ##
##################################################################
mcmc_results_out_dir = "C:\\Users\\putnni\\switchdrive\\Chad\\mcmc-results\\"
mcmc_results <- readRDS(paste0(mcmc_results_out_dir, "mcmc_results_final.rds"))

# Sampling from the posterior distribution
n_samples <- 500
param_samples <- sample_mcmc_steps(mcmc_results$coda_pars, n_samples)
param_samples <- param_samples[,c("sigma_LT", "sigma_RT")]

max_posterior_params <- extract_max_posterior_params(mcmc_results)
best_sigma_LT <- max_posterior_params["sigma_LT"]
best_sigma_RT <- max_posterior_params["sigma_RT"]

# Final simplified EIR function with parameterized saturation on X
combined_eir <- function(temp, rain, T_opt, sigma_LT, sigma_RT,
                         R_opt, k, X = 0.1, c_X, d_X, alpha = 1) {
  # Temperature effect (asymmetric Gaussian)
  temp_effect <- ifelse(
    temp <= T_opt,
    exp(-((temp - T_opt)^2) / (2 * sigma_LT^2)),
    exp(-((temp - T_opt)^2) / (2 * sigma_RT^2))
  )

  # Rainfall effect (logistic)
  rain_effect <- 1 / (1 + exp(-k * (rain - R_opt)))

  # Saturated X effect
  X_effect <- X / (c_X + d_X * X)

  # Final EIR
  eir <- alpha * temp_effect * rain_effect * X_effect
  return(eir)
}

#############
## Settings ##
#############
T_opt      <- 27.818
rain_fixed <- 1
X_fixed    <- 0.1
temp_seq   <- seq(18, 35, by = 0.1)

#############################################
## 1. Compute Mechanistic & Best‐Fit Curves ##
#############################################
# 1a. Mechanistic EIR (normalized to [0,1])
eir_mech <- sapply(temp_seq, function(t) {
  EIR_f(t, p_EA_R = p_EA_R_f(rain_fixed), X = X_fixed)
})
eir_mech <- eir_mech / max(eir_mech)

# 1b. Simplified best‐fit EIR (normalized to [0,1])
eir_inferred <- combined_eir(
  temp     = temp_seq,
  rain     = rain_fixed,
  T_opt    = T_opt,
  sigma_LT = best_sigma_LT,
  sigma_RT = best_sigma_RT,
  R_opt    = rain_fixed,
  k        = mcmc_results$param_inputs$k1,
  X        = X_fixed,
  c_X      = 0.01,
  d_X      = 1.0,
  alpha    = 1
)
eir_inferred <- eir_inferred / max(eir_inferred)

#########################################
## 2. Simulate Posterior‐Draw Curves   ##
#########################################
# generate 500 normalized EIR curves
eir_mat <- apply(param_samples, 1, function(par) {
  tmp <- combined_eir(
    temp     = temp_seq,
    rain     = rain_fixed,
    T_opt    = T_opt,
    sigma_LT = par["sigma_LT"],
    sigma_RT = par["sigma_RT"],
    R_opt    = rain_fixed,
    k        = mcmc_results$param_inputs$k1,
    X        = X_fixed,
    c_X      = 0.01,
    d_X      = 1.0,
    alpha    = 1
  )
  tmp / max(tmp)
})

#########################################################
## 3. 95% Credible Ribbon from Posterior Simulations  ##
#########################################################
eir_ci <- apply(eir_mat, 1, quantile, probs = c(0.025, 0.975))
ribbon_df <- tibble(
  Temperature = temp_seq,
  ymin        = eir_ci[1, ],
  ymax        = eir_ci[2, ]
)

##########################################
## 4. Prepare Data for Line Plotting    ##
##########################################
eir_temp_df <- tibble(
  Temperature = temp_seq,
  `Ukawuba 2022`           = eir_mech,
  `Simplified (estimated)` = eir_inferred
) %>%
  pivot_longer(
    cols      = c(`Ukawuba 2022`, `Simplified (estimated)`),
    names_to  = "Model",
    values_to = "EIR"
  )

##############################
## 5. Plot: Ribbon + Lines ##
##############################
p_inferred_original <- ggplot() +
  # ribbon first, so lines sit on top
  geom_ribbon(
    data  = ribbon_df,
    aes(x = Temperature, ymin = ymin, ymax = ymax),
    fill  = "grey80",
    alpha = 0.5
  ) +
  # mechanistic & best‐fit lines
  geom_line(
    data = eir_temp_df,
    aes(x = Temperature, y = EIR, color = Model),
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      "Ukawuba 2022"           = "firebrick",
      "Simplified (estimated)" = "steelblue"
    )
  ) +
  labs(
    x     = "Temperature (°C)",
    y     = "Daily EIR (scaled)",
    color = NULL
  ) +
  theme_pubr() +
  theme(
    legend.text     = element_text(size = 13),
    legend.title    = element_text(face = "bold")
  )

p_inferred_original
