#############
## Overview #
#############
# Here the goal is to calculate a population growth rate, by age group
# We assume the growth rate in the district of Moissala is the same as in Chad overall
# This is because we do not have reliable population data for Moissala

##########
## Setup #
##########
library(readxl)
library(dplyr)
path_to_chad <- "C:/Users/putnni/switchdrive/Chad/"

####################
## Reading in data #
####################
# Reading in Chad population data
# Data comes from UN, World Population Prospects 2024 (link to data: https://ourworldindata.org/explorers/population-and-demography?country=~TCD&indicator=Population+by+broad+age+group&Sex=Both+sexes&Age=Total&Projection+scenario=High)
chad_pop_data <- read.csv(paste0(path_to_chad, "data/demographic-data/population-by-broad-age-group.csv"))

head(chad_pop_data)

####################
## Data processing #
####################
chad_pop_data <- chad_pop_data %>% filter(Year %in% 2018:2024)
chad_pop_data <- chad_pop_data %>% reframe(year = Year,
                                  `<5 yrs` = X0.4.years,
                                  `>=5 yrs` = X5.14.years + X15.24.years + X25.64.years + X65..years,
                                  total = `<5 yrs` + `>=5 yrs`)

#####################################
## Calculating growth rate for Chad #
#####################################
yearly_growth_rate <- function(p0, pn, n){
  return(log(pn / p0) / n)
}

# Yearly growth rates

# total
p0_tot = chad_pop_data$total[1]
pn_tot = chad_pop_data$total[nrow(chad_pop_data)]
n_tot <- max(chad_pop_data$year) - min(chad_pop_data$year)

r_tot <- yearly_growth_rate(p0_tot, pn_tot, n_tot)

# <5
p0_u5 = chad_pop_data$`<5 yrs`[1]
pn_u5 = chad_pop_data$`<5 yrs`[nrow(chad_pop_data)]
n_u5 <- max(chad_pop_data$year) - min(chad_pop_data$year)

r_u5 <- yearly_growth_rate(p0_u5, pn_u5, n_u5)

# >=5
p0_o5 = chad_pop_data$`>=5 yrs`[1]
pn_o5 = chad_pop_data$`>=5 yrs`[nrow(chad_pop_data)]
n_o5 <- max(chad_pop_data$year) - min(chad_pop_data$year)

r_o5 <- yearly_growth_rate(p0_o5, pn_o5, n_o5)

cat("r_tot =", round(r_tot, 4))
cat("r_u5 =", round(r_u5, 4))
cat("r_o5 =", round(r_o5, 4))

# Daily growth rates
cat("daily multiplier total:", (1 + r_tot)^(1/365))
cat("daily multiplier u5:", (1 + r_u5)^(1/365))
cat("daily multiplier o5:", (1 + r_o5)^(1/365))

r_daily_df <- data.frame(r_daily_tot = (1 + r_tot)^(1/365),
           r_daily_u5 = (1 + r_u5)^(1/365),
           r_daily_o5 = (1 + r_o5)^(1/365))

saveRDS(r_daily_df, paste0(path_to_chad, "data/model-inputs/daily_growth_rates_chad.rds"))


1 + yearly_growth_rate(p0_o5, pn_o5, n_o5 * 365.25)
1 + yearly_growth_rate(p0_u5, pn_u5, n_u5 * 365.25)
1 + yearly_growth_rate(p0_tot, pn_tot, n_u5 * 365.25)
