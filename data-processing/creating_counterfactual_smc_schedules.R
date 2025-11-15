library(readxl)
library(lubridate)
path_to_chad <- "C:/Users/putnni/switchdrive/Chad/"
SMC_raw <- read_excel(file.path(path_to_chad, "Data/data-raw/CPS/CPS_coverage_imput_2018_with_2023_no_duplicates.xlsx"))

# Calculating average coverage per round
cov_by_round <- SMC_raw %>% group_by(smc_round) %>%
  filter(year(date_start) != 2019, smc_couv_tot > 0) %>%
  summarize(smc_couv_tot = mean(smc_couv_tot))

# SMC in 2019
SMC_in_2019 <- SMC_raw
SMC_in_2019[year(SMC_in_2019$date_start) == 2019,]$smc_couv_tot <- cov_by_round[1:4,]$smc_couv_tot
saveRDS(SMC_in_2019, file.path(path_to_chad, "Data/smc-inputs/SMC_in_2019.rds"))

# Four rounds of SMC in 2021 and 2022
SMC_2021_2022_cf <- SMC_raw
SMC_2021_2022_cf <- SMC_2021_2022_cf %>%
  mutate(smc_couv_tot = if_else(
    (year(date_start) %in% 2021:2022) & (month(date_start) == 11),
    0,
    smc_couv_tot
  ))

saveRDS(SMC_2021_2022_cf, file.path(path_to_chad, "Data/smc-inputs/SMC_2021_2022_4rounds.rds"))

# SMC starting in July in 2023
SMC_2023_July <- SMC_raw
smc_couv_2023 <- SMC_2023_July %>% filter(year(date_start) == 2023, month(date_start) %in% 6:10) %>% pull(smc_couv_tot)

SMC_2023_July <- SMC_2023_July %>%
  mutate(
    smc_couv_tot = replace(smc_couv_tot,
                           which(year(date_start) == 2023 & month(date_start) == 6),
                           0),
    smc_couv_tot = replace(smc_couv_tot,
                           which(year(date_start) == 2023 & month(date_start) %in% 7:11),
                           smc_couv_2023)
  )

saveRDS(SMC_2023_July, file.path(path_to_chad, "Data/smc-inputs/SMC_2023_July.rds"))
