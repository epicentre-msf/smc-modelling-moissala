#############################################################
##   Climate-Entomological Modeling: Setup and Preparation ##
#############################################################
# This script shows how the EIR-climate relationships were simplified

##############################
##  Load Required Libraries ##
##############################
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

##############################
##  Reading in climate data ##
##############################
# Latitude and longitude where climate data (rainfall and temperature) is to be saved
# Rainfall data is from CHIRPS and temperature data is from ERA5
lat <- 8.34
lon <- 17.77

path_to_chad <- "C:/Users/putnni/switchdrive/Chad/" # Replace this with location of chad switchdrive folder

# Reading in climate data saved by `save_climate_data`
temp_path <- paste0(path_to_chad, "data/climate-data/era5_moissala.nc")
rain_path <- paste0(path_to_chad, "data/climate-data/chirps_moissala.rds")

met_365 <- process_climate_data(lon, lat, years = 2018:2023, temp_path = temp_path,
                                rain_path = rain_path, months_30_days = FALSE)


#################################################
##  Climate–Entomological Function Definitions ##
#################################################
path_to_repo <- "C:/Users/putnni/Documents/git-hub-repositories/moissala-smc/"
source(paste0(path_to_repo, "paper/supplementary-analysis/original_climate_equation_definitions.R"))

temps <- seq(0, 40, by = 0.1)

# Output temperature at max EIR
temps[which.max(sapply(temps, function(t) EIR_f(t, p_EA_R = 0.8)))]

# Reset layout
par(mfrow = c(1, 1))


#################################################
## EIR vs Temperature Using Original Equations ##
#################################################
# Define temperature range (in °C)
temp_range <- seq(16, 38, by = 0.01)

# Set fixed covariate values
fixed_rain <- 0      # Standardized rainfall
fixed_X    <- 0.1    # Proportion of infectious humans

# Compute mechanistic EIR for each temperature
eir_temp <- sapply(temp_range, function(temp) {
  EIR_f(temp, p_EA_R = p_EA_R_f(fixed_rain), X = fixed_X)
})

# Package results in a data frame
temp_data <- data.frame(Temperature = temp_range, EIR = eir_temp)

optimal_temp <- temp_data[which.max(temp_data$EIR),]$Temperature

# Plot: Mechanistic EIR vs Temperature
p_mech_EIR_temp <- ggplot(temp_data, aes(x = Temperature, y = EIR)) +
  geom_line(color = "red", size = 1) +
  labs(
    title = "",
    x = "Temperature (°C)", y = "Daily EIR"
  ) +
  # annotation of sigma values
  annotate(
    "text",
    x     = 25.5,
    y     = 0.3,
    label = as.expression(
      bquote(T[opt] == .(round(optimal_temp, 2)))
    ), size = 4,
    color = "black",
    hjust = 0
  )  + theme_pubr() + xlim(20, 35) +
  scale_x_continuous(breaks = seq(15, 35, 2)) +
  coord_cartesian(xlim = c(19,35), ylim = c(0, 0.3))

p_mech_EIR_temp

max_temp_for_transmission <- temp_data[which(temp_data$EIR == 0),]$Temperature[1]

########################################################################
## Proportion of "days" (using smoothed estimates) excess of max temp ##
########################################################################
mean(met_365$temp > max_temp_for_transmission)

# Bin temps into factor categories
met_binned <- met_365 %>%
  mutate(
    temp_cat = cut(
      temp,
      breaks = c(-Inf, 26, 28, 30, max_temp_for_transmission, Inf),
      labels = c(
        "below 26°C",
        "26–28°C",
        "28–30°C",
        paste0("30–", max_temp_for_transmission, "°C"),
        paste0("above ", max_temp_for_transmission, "°C")
      ),
      right = FALSE
    )
  )

# Count & compute proportions
prop_df <- met_binned %>%
  count(temp_cat) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop)
  )

# Plot
p_percent_hot_days <- ggplot(prop_df, aes(x = temp_cat, y = prop)) +
  geom_col(fill = "#2c7fb8", color = "black", width = 0.7) +
  geom_text(aes(label = label), vjust = -0.2, size = 4) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(c(0, 0.1))) +
  labs(
    x     = "",
    y     = "Percentage of days",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold")
  )

p_percent_hot_days + theme_pubr()

#################################################################
## Simplified Functional Form for Temperature–EIR Relationship ##
#################################################################
# This simplified form approximates the mechanistic EIR-temperature curve
# using an asymmetric Gaussian centered at T_opt, with different left/right widths.

asymmetric_gaussian <- function(x, sigma_LT, sigma_RT, T_opt = 27.818, scale = 1) {
  scale * ifelse(
    x <= T_opt,
    exp(-((x - T_opt)^2) / (2 * sigma_LT^2)),
    exp(-((x - T_opt)^2) / (2 * sigma_RT^2))
  )
}

##############################################################
## EIR vs Temperature: Original vs Simplified Comparison ##
##############################################################
# Extract fitted parameters
sigma_LT     <- fit$par[1]
sigma_RT     <- fit$par[2]
scale_factor <- max(EIR_true$EIR)  # Bring predictions to original scale

#############
## Settings ##
#############
T_opt      <- 27.818
rain_fixed <- 1
X_fixed    <- 0.1
new_temps  <- seq(18, 35, by = 0.1)

###############################################################
## Fit Asymmetric Gaussian to Output From Original Equations ##
###############################################################
temps <- seq(18, 33, length.out = 100)
EIR_mech <- EIR_f(temps)
EIR_scaled <- EIR_mech / max(EIR_mech)
EIR_true <- data.frame(temp = temps, EIR = EIR_mech, EIR_scaled = EIR_scaled)

fit <- optim(
  par    = c(sigma_LT = 1, sigma_RT = 1),
  fn     = function(params) {
    pred <- asymmetric_gaussian(EIR_true$temp, params[1], params[2], T_opt)
    sum((EIR_true$EIR_scaled - pred)^2)
  },
  method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(10, 10)
)

sigma_LT <- fit$par[1]
sigma_RT <- fit$par[2]

##########################################################
## Compute Mechanistic & Simplified Curves (normalized) ##
##########################################################
# 1a. Mechanistic (normalized to max=1)
EIR_mech_2    <- EIR_f(new_temps, p_EA_R = p_EA_R_f(rain_fixed), X = X_fixed)
EIR_scaled_2  <- EIR_mech_2 / max(EIR_mech_2)

# 1b. Simplified (use scale = 1 so output is on same 0–1 scale)
predicted_EIR <- asymmetric_gaussian(
  new_temps, sigma_LT, sigma_RT, T_opt, scale = 1
)

################################
## Build tidy data for ggplot ##
################################
df <- tibble(
  temp  = new_temps,
  mech  = EIR_scaled_2,
  simp  = predicted_EIR
) %>%
  pivot_longer(
    cols      = c(mech, simp),
    names_to  = "model",
    values_to = "EIR"
  ) %>%
  mutate(
    model = recode(model,
                   mech = "Ukawuba 2022",
                   simp = "Simplified")
  )

##################################################
## Plot: Overlaid Normalized Curves          ##
##################################################
p_approx_orig <- ggplot(df, aes(x = temp, y = EIR, color = model)) +
  geom_line(size = 1.2) +
  annotate(
    "text",
    x     = 21.5,
    y     = 0.75,
    label = as.expression(
      bquote(sigma[LT] == .(round(sigma_LT, 2)))
    ),
    size  = 4, color = "black", hjust = 0
  ) +
  annotate(
    "text",
    x     = 30.7,
    y     = 0.75,
    label = as.expression(
      bquote(sigma[RT] == .(round(sigma_RT, 2)))
    ),
    size  = 4, color = "black", hjust = 0
  ) +
  scale_color_manual(
    values = c("Ukawuba 2022" = "firebrick", "Simplified" = "steelblue")
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Daily EIR (scaled)"
  ) +
  theme_pubr() +
  theme(
    legend.text  = element_text(size = 13),
    legend.title = element_blank()
  )

p_approx_orig
############################
## Saving Relevant Plots  ##
############################
# Location for plots to be saved
path_to_figures <- paste0(path_to_repo, "paper/paper/plots/supplemental-plots/")
source(paste0(path_to_repo, "paper/supplementary-analysis/comparing_estimated_EIR_to_original.R"))

p_mech <- p_mech_EIR_temp + theme_pubr() + ggtitle("A") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))
p_mech

p_hot_days <- p_percent_hot_days + theme_pubr() + ggtitle("B") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

p_hot_days

common_xlim <- c(18, 35)
common_ylim <- c(0, 1.1)
p_app_orig <- p_approx_orig + theme_pubr() + ggtitle("C") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold")) +
  coord_cartesian(xlim = common_xlim, ylim = common_ylim)

p_app_orig

p_inf_orig <- p_inferred_original + theme_pubr() + ggtitle("D") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  coord_cartesian(xlim = common_xlim, ylim = common_ylim)
p_inf_orig


combined_eir_simplified <- (p_mech + p_hot_days) / (p_app_orig + p_inf_orig)

ggsave(filename = paste0(path_to_figures, "eir_simplified_figure.pdf"),
       plot = combined_eir_simplified,
       height = 10,
       width = 13,
       dpi = 500)

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #
# @@@@@@@@@@@@@@@ ADDITIONAL ANALYSIS NOT FOUND IN PAPER  @@@@@@@@@@@@@@@@@@@@ #
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ #

###################################################################
## Fit Saturating Term for Proportion Infectious (X) in EIR Fit ##
###################################################################
# Fix temperature and rainfall (same as in temp fitting section)
temp_fixed <- 27.818
rain_fixed <- 1

# Generate sequence of X values
X_seq <- seq(0.01, 1, by = 0.01)

# Mechanistic EIR vs X
EIR_mech_X <- sapply(X_seq, function(X) {
  EIR_f(temp_fixed, p_EA_R = p_EA_R_f(rain_fixed), X = X)
})

# Rescale mechanistic EIR for fitting
EIR_scaled_X <- EIR_mech_X / max(EIR_mech_X)
EIR_X_true <- data.frame(X = X_seq, EIR = EIR_mech_X, EIR_scaled = EIR_scaled_X)

# Define generalized saturating function: X / (c + d·X)
saturating_X <- function(X, c, d, scale = 1) {
  scale * (X / (c + d * X))
}

# Optimization to find best-fitting (c, d)
fit_X <- optim(
  par = c(c = 1, d = 1),
  fn = function(params) {
    c_val <- params[1]
    d_val <- params[2]
    pred <- saturating_X(EIR_X_true$X, c = c_val, d = d_val, scale = 1)
    sum((EIR_X_true$EIR_scaled - pred)^2)
  },
  method = "L-BFGS-B",
  lower = c(1e-10, 1e-10),
  upper = c(200, 200)
)

# Extract fitted parameters
c_opt <- fit_X$par[1]
d_opt <- fit_X$par[2]
scale_X <- max(EIR_X_true$EIR)

# Predicted values using fitted parameters and original scale
predicted_EIR_X <- saturating_X(EIR_X_true$X, c = c_opt, d = d_opt, scale = scale_X)

##################################################
## EIR vs X: Original vs Simplified Comparison ##
##################################################
# Plot: Mechanistic vs Simplified Fit
plot(EIR_X_true$X, EIR_X_true$EIR, pch = 19, col = "red",
     xlab = "Proportion Infectious (X)", ylab = "EIR",
     main = "Original vs Simplified X Fit")

lines(EIR_X_true$X, predicted_EIR_X, col = "blue", lwd = 2)

# Annotate fitted parameter values
text(0.2, scale_X * 0.8, paste0("c = ", round(c_opt, 5)), col = "blue", pos = 3)
text(0.2, scale_X * 0.7, paste0("d = ", round(d_opt, 5)), col = "blue", pos = 3)

legend("topright",
       legend = c("Original EIR", "Fitted Saturated X"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))


################################################################################
## Final Proposed EIR Function: Gaussian(Temp) × Logistic(Rain) × Saturated X ##
################################################################################
# This function formalizes the proposed simplified EIR function:
# EIR(temp, rain, X) = α × AsymGaussian(temp) × Logistic(rain) × (X / (c + d·X))

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

# Parameters from previous optimizations
alpha     <- max(sapply(seq(20, 35, 0.2), function(t) EIR_f(t, p_EA_R = p_EA_R_f(2), X = 0.1)))
T_opt     <- 27.818
sigma_LT  <- sigma_LT
sigma_RT  <- sigma_RT
R_opt     <- 2.0   # aligned with p_EA_R_f
k         <- 0.5   # logistic slope
c_X       <- c_opt # from earlier X-fit
d_X       <- d_opt
X_fixed   <- 0.1

# Use grid of typical values
temp_vals <- seq(22, 32, by = 0.5)
rain_vals <- seq(-1, 3, by = 0.5)
X_vals <- seq(0.05, 0.2, by = 0.05)

# Mechanistic values
eir_mech_grid <- expand.grid(temp = temp_vals, rain = rain_vals, X = X_vals) %>%
  rowwise() %>%
  mutate(EIR = EIR_f(temp, p_EA_R = p_EA_R_f(rain), X = X)) %>%
  ungroup()

# Simplified values (without alpha)
eir_simp_grid <- expand.grid(temp = temp_vals, rain = rain_vals, X = X_vals) %>%
  rowwise() %>%
  mutate(EIR = combined_eir(temp, rain, T_opt, sigma_LT, sigma_RT,
                            R_opt, k, X, c_X, d_X, alpha = 1)) %>%
  ungroup()

# Estimate scale factor alpha
alpha_opt <- max(eir_mech_grid$EIR) / max(eir_simp_grid$EIR)

######################################
## Scale Simplified EIR with alpha ##
######################################
# Add alpha_opt to simplified values
eir_simp_grid <- eir_simp_grid %>%
  mutate(EIR_scaled = EIR * alpha_opt)

# Merge for comparison
comparison_df <- eir_mech_grid %>%
  rename(Mechanistic = EIR) %>%
  left_join(eir_simp_grid %>% select(temp, rain, X, Simplified = EIR_scaled),
            by = c("temp", "rain", "X")) %>%
  pivot_longer(cols = c(Mechanistic, Simplified), names_to = "Model", values_to = "EIR")

####################################
## Plot: EIR vs Temperature       ##
####################################
ggplot(comparison_df %>% filter(rain == 0, X == 0.1),
       aes(x = temp, y = EIR, color = Model)) +
  geom_line(size = 1) +
  labs(title = "EIR vs Temperature (rain = 0, X = 0.1)",
       x = "Temperature (°C)", y = "Daily EIR") +
  theme_minimal()

####################################
## Plot: EIR vs Rainfall          ##
####################################
ggplot(comparison_df %>% filter(temp == 27, X == 0.1),
       aes(x = rain, y = EIR, color = Model)) +
  geom_line(size = 1) +
  labs(title = "EIR vs Rainfall (temp = 27°C, X = 0.1)",
       x = "Standardized Rainfall", y = "Daily EIR") +
  theme_minimal()

####################################
## Plot: EIR vs Proportion Infectious ##
####################################
ggplot(comparison_df %>% filter(temp == 27, rain == 0),
       aes(x = X, y = EIR, color = Model)) +
  geom_line(size = 1) +
  labs(title = "EIR vs Proportion Infectious (temp = 27°C, rain = 0)",
       x = "Proportion Infectious (X)", y = "Daily EIR") +
  theme_minimal()

######################################
## Final Simplified EIR Parameters  ##
######################################
cat(sprintf("alpha = %.4e # Scaling constant\n", alpha_opt))
cat(sprintf("T_opt = %.2f # Optimal temperature (°C)\n", T_opt))
cat(sprintf("sigma_LT = %.2f # Left-side std dev (temp < T_opt)\n", sigma_LT))
cat(sprintf("sigma_RT = %.2f # Right-side std dev (temp > T_opt)\n", sigma_RT))
cat(sprintf("R_opt = %.2f # Rainfall inflection point (standardized)\n", R_opt))
cat(sprintf("k = %.2f # Logistic slope for rainfall effect\n", k))
cat(sprintf("c_X = %.4f # X saturation intercept\n", c_X))
cat(sprintf("d_X = %.4f # X saturation slope\n", d_X))
