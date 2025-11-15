##########################
## Wrapper function     ##
##########################
process_msf_cases <- function(data, service_type, save_path, cds_name_filter = NULL, monthly = FALSE, start_date = NULL) {
  # Optional filtering by cds_name
  if (!is.null(cds_name_filter)) {
    data <- data %>% filter(cds_name %in% cds_name_filter)
  }

  # Aggregate n_cases by date_sunday, service, and age_group
  data <- data %>%
    group_by(date_sunday, service, age_group) %>%
    summarise(n_cases = sum(n_cases, na.rm = TRUE), .groups = "drop")

  # Filter by age group and service
  data_u5 <- data %>% filter(age_group == "0to4", service == service_type) %>% select(date_sunday, n_cases)

  data_combined <- data_u5

  # Rename and create additional columns
  colnames(data_combined) <- c("date_ymd", "inc_C")

  if (monthly) {
    data_combined <- data_combined %>%
      mutate(month_start = floor_date(date_ymd, unit = "month")) %>%
      group_by(month_start) %>%
      summarise(
        inc_C = sum(inc_C),
        .groups = "drop"
      ) %>%
      rename(date_ymd = month_start)
  } else {
    data_combined <- data_combined %>%
      mutate(week_no = 0:(n() - 1)) %>%
      select(date_ymd, week_no, inc_C)
  }

  if(!is.null(start_date)){
    data_combined <- data_combined %>% filter(date_ymd >= start_date)
  }

  # Save to file
  saveRDS(data_combined, save_path)

  return(data_combined)
}

###############
## Setup     ##
###############
library(dplyr)
library(lubridate)

path_to_chad <- "C:/Users/putne/switchdrive/Chad/" # Replace with location of Chad repo
MSF_district <- readRDS(paste0(path_to_chad, "Data/data-raw/MSF/data_district_msf_excel_dhis2_opd_ipd_2012_2023_age.rds"))
case_data_path <- paste0(path_to_chad, "Data/cases-data/MSF-data/")

# Outpatient department data - weekly
opd_week_path <- paste0(case_data_path, "opd_district_MOISSALA_2018_2023.rds")
process_msf_cases(MSF_district, "opd",
                  save_path = opd_week_path,
                  cds_name_filter = NULL, monthly = FALSE,
                  start_date = "2018-01-01")
