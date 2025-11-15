##########################################################
##  Climate–Entomological Function Definitions           ##
##  --------------------------------------------------  ##
##  These functions define mechanistic relationships     ##
##  between climate variables (temperature, rainfall)    ##
##  and key entomological quantities that affect         ##
##  malaria transmission dynamics.                       ##
##########################################################

# Biting rate helper used internally by several functions
biting_rate_f <- function(temp) pmax(0.017 * temp - 0.165, 0)

# Water temperature
temp_w_f <- function(temp, k = 1, delta_temp = 2) pmax(k * temp + delta_temp, 0)

# Egg-to-adult development time
tau_EA_f <- function(temp_w) ifelse(
  temp_w <= 14.7 | temp_w >= 34, Inf,
  1 / (0.000111 * temp_w * (temp_w - 14.7) * sqrt(34 - temp_w))
)

# Egg-to-adult survivorship (temperature-dependent)
p_EA_T_f <- function(temp_w) pmax(-0.00924 * temp_w^2 + 0.453 * temp_w - 4.77, 0)

# Sporogony (parasite development)
n_f <- function(temp) ifelse(
  temp <= 15.384 | temp >= 35, Inf,
  1 / (0.000112 * temp * (temp - 15.384) * sqrt(35 - temp))
)

# Biting rate
a_f <- biting_rate_f

# Gonotropic cycle (days between bites)
GP_f <- function(temp) {
  a <- biting_rate_f(temp)
  ifelse(a <= 0, Inf, 1 / a)
}

# Lifetime egg production
B_f <- function(temp, e = 50, g = -log(0.93)) {
  a <- biting_rate_f(temp)
  ifelse(a <= 0, 0, e / (exp((1 / a) * g) - 1))
}

# Mosquito density
m_f <- function(temp, e = 50, g = -log(0.93), p_EA_R = 0.50) {
  temp_w <- temp_w_f(temp)
  B_egg <- B_f(temp, e = e, g = g)
  p_EA <- p_EA_T_f(temp_w) * p_EA_R
  tau_EA <- tau_EA_f(temp_w)
  ifelse(tau_EA <= 0 | g <= 0, 0, (B_egg * p_EA) / (tau_EA * g))
}


# Entomological Inoculation Rate (EIR)
EIR_f <- function(temp, p_HM = 0.125, g = -log(0.93), X = 0.1, p_EA_R = 0.30) {
  a <- a_f(temp)
  m <- m_f(temp, p_EA_R = p_EA_R)
  n <- n_f(temp)
  #B <- B_f(temp)
  ifelse(g + (a * p_HM) <= 0 | is.infinite(n), 0,
         (m * a^2 * exp(-g * n) * p_HM * X) / (g + (a * p_HM * X))) |> pmax(0)
}

# Rainfall-based survival
p_EA_R_f <- function(rain, a_R = 0.5, b_R = 2) {
  pmax(1 / (1 + exp(-a_R * (rain - b_R))), 0)
}


######################################################################
##  Visualizing Entomological Parameters as Functions of Temperature##
######################################################################

temps <- seq(0, 40, by = 0.1)

# Helper to clean invalid values and plot
plot_temp_response <- function(var_name, values, color = "black", ylab = var_name) {
  values[!is.finite(values) | values < 0] <- NA
  plot(temps, values, type = "l", col = color, lwd = 2,
       xlab = "Temperature (°C)", ylab = ylab, main = var_name,
       ylim = get_adjusted_ylim(values))
}

get_adjusted_ylim <- function(values, buffer = 0.1) {
  finite_values <- values[is.finite(values)]
  lower <- min(finite_values, na.rm = TRUE)
  upper <- quantile(finite_values, 0.95, na.rm = TRUE)
  range <- upper - lower
  c(lower - buffer * range, upper + buffer * range)
}

# Set 3x3 grid layout
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), cex.axis = 0.8, cex.lab = 0.8)

# Generate and plot each variable
plot_temp_response("Water Temperature", temp_w_f(temps), "blue")
plot_temp_response("Egg-Adult Development Time", tau_EA_f(temp_w_f(temps)), "green", "Development Time (days)")
plot_temp_response("Egg-to-Adult Survivorship", p_EA_T_f(temp_w_f(temps)), "orange", "Survivorship")
plot_temp_response("Sporogony", n_f(temps), "purple", "Sporogony (days)")
plot_temp_response("Mosquito Biting Rate", a_f(temps), "red", "Biting Rate")
plot_temp_response("Gonotropic Cycle", GP_f(temps), "brown", "Cycle Length (days)")
plot_temp_response("Lifetime Number of Eggs", B_f(temps), "cyan", "Lifetime Eggs")
plot_temp_response("Mosquito Density", sapply(temps, m_f), "darkgreen", "Density")
plot_temp_response("Entomological Inoculation Rate", sapply(temps, function(t) EIR_f(t, p_EA_R = 0.8)), "magenta", "EIR")

