####################################
##      Load Required Libraries   ##
####################################
library(deSolve)  # For solving differential equations
library(tidyr)    # For data manipulation (pivot_longer, etc.)
library(ggplot2)  # For plotting
library(dplyr)    # For data manipulation (mutate, etc.)

####################################
##     Demographic ODE Model      ##
####################################
# This function defines the system of differential equations modeling
# population dynamics split between children (C) and adults (A)
demography <- function(t, x, parms){
  C = x[1]  # Number of children
  A = x[2]  # Number of adults

  b = parms["b"]     # Birth (and death) rate
  age = parms["age"] # Aging rate (1 / average time as a child)

  # dC: new children born, aging into adulthood, and dying
  # dA: children aging into adulthood, and adults dying
  dC = (A + C) * b - age * C - b * C
  dA = age * C - b * A

  res = c(dC, dA)
  list(res)
}

####################################
##     Simulation Function        ##
####################################
simulate_demo <- function(parameters){
  parms = c(parameters["b"], parameters["age"])  # Extract parameters
  init <- c(parameters["initC"], parameters["initA"])  # Initial conditions

  # Simulate the system over 10,000 days
  temps <- seq(0, 10000)
  solveSIR <- lsoda(y = init, times = temps, func = demography, parms = parms)
  solution = as.data.frame(solveSIR)
  names(solution) = c("time", "C", "A")  # Rename columns for clarity

  return(solution)
}

##############################################
##     Estimate Child Proportion at Equilibrium   ##
##############################################
# This function simulates the system for a grid of birth rates,
# then calculates the stabilized proportion of children (C / (C + A))
calc_prop_C <- function(grid){
  prop_c_vec <- array(dim = length(b_vec))  # Vector to store child proportions

  for(i in 1:length(b_vec)){
    theta_init = c("b" = b_vec[i],
                   "age" = 1 / (5 * 365),  # Average child age: 5 years
                   "initC" = 30,
                   "initA" = 130)

    simul = simulate_demo(parameters = theta_init)
    prop_c <- (simul$C / (simul$C + simul$A))[10000]  # Get final time point
    prop_c_vec[i] <- prop_c
  }

  return(prop_c_vec)
}

####################################
##     Run Estimation Procedure   ##
####################################
b_vec <- (10:80) / (365 * 1000)  # Grid of birth/death rates to test
prop_C_vec <- calc_prop_C(b_vec)  # Corresponding child proportions

# Find the index of the birth rate that gives ~19% children
best_ind <- which.min(abs(prop_C_vec - 0.19))
best_val <- b_vec[best_ind]  # This is approximately 1 / (365 * 21)
(10:80)[best_ind]  # Return raw index value from original vector (for inspection)

####################################
##     Final Simulation & Plot    ##
####################################
# Use the estimated birth rate and simulate final model
theta_init = c("b" = best_val,
               "age" = 1 / (5 * 365),  # Average child age = 5 years
               "initC" = 10,
               "initA" = 140)
simul = simulate_demo(parameters = theta_init)

# Plot the proportion of children and adults over time
simul %>%
  mutate(N = C + A,
         prop_C = C / N,
         prop_A = A / N) %>%
  pivot_longer(cols = c(prop_C, prop_A)) %>%
  ggplot() +
  geom_line(aes(x = time, y = value, color = name), lwd = 1)

