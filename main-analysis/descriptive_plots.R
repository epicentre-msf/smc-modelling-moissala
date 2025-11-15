#######################
## LOADING LIBRARIES ##
#######################
library(readxl)
#library(malclimsim)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(gtable)
library(gridExtra)
library(ggtext)

path_to_chad <- "C:/Users/putne/switchdrive/Chad/"
path_to_repo <- "C:/Users/putne/OneDrive/Documents/git-repos/moissala-smc/"

devtools::load_all(path = "C:/Users/putne/OneDrive/Documents/git-repos/malclimsim")

source(paste0(path_to_repo, "paper-code/aux-functions/aux_functions.R"))
source(paste0(path_to_repo, "paper-code/aux-functions/plotting_functions.R"))

###########
## SETUP ##
###########
out_dir <- file.path(path_to_repo, "paper-plots/descriptive-plots/")

plot_width <- 8
plot_height <- 3
plot_dpi <- 1000

start_date <- as.Date("2018-01-01")
end_date   <- as.Date("2023-12-01")
years_analysis <- 2018:2023

##########################################
## DATA LOADING
## -- Load cases, SMC dates, and rainfall
##########################################
SMC_df <- read_excel(file.path(path_to_chad, "Data/data-raw/CPS/CPS_coverage_imput_2018_with_2023_no_duplicates.xlsx")) %>%
  mutate(date_start = as.Date(date_start))

obs_cases <- readRDS(file.path(path_to_chad, "Data/cases-data/MSF-data/opd_district_MOISSALA_2012_2023.rds")) %>%
  filter(date_ymd %in% as.Date((start_date : end_date)))


# Reading in climate data saved by `save_climate_data`
temp_path <- paste0(path_to_chad, "data/climate-data/era5_moissala.nc")
rain_path <- paste0(path_to_chad, "data/climate-data/chirps_moissala.rds")

monthly_temp_df <- extract_era5(lat = 8.34, lon = 17.77, path_to_file = temp_path) %>% filter(year(Date) %in% years_analysis)
monthly_temp_df$Date <- as.Date(monthly_temp_df$Date)
rain_df <- readRDS(rain_path) %>% filter(year(date) %in% years_analysis)

colnames(monthly_temp_df) <- c("dates", "temp")
colnames(rain_df) <- c("dates", "rainfall")

monthly_rain_df <- rain_df %>%
  mutate(year = year(dates),
         month = month(dates)) %>%
  group_by(year, month) %>%
  summarise(rainfall = sum(rainfall, na.rm = TRUE), .groups = "drop") %>%
  mutate(dates = as.Date(sprintf("%04d-%02d-01", year, month))) %>%
  select(dates, rainfall) %>%
  arrange(dates)

monthly_rain_df <- monthly_rain_df %>% filter(dates >= min(obs_cases$date_ymd))
monthly_temp_df <- monthly_temp_df %>% filter(dates >= min(obs_cases$date_ymd))

plot(monthly_rain_df$rainfall, type = "b")
plot(monthly_temp_df$temp, type = "b")

###############################
## SHARED DEFINITIONS
###############################
date_breaks = "2 months"
y_axis_title_size = 14
x_axis_label_size = 16

###############################
## PLOT: Weekly inc_C with SMC Shading (Line + Bars)
###############################
plot_width <- 7
plot_height <- 4
base_size = 12

# Prepare SMC season rectangles
smc_rects <- SMC_df %>%
  filter(year(date_start) %in% years_analysis) %>%
  mutate(year = year(date_start)) %>%
  group_by(year) %>%
  summarize(
    xmin = min(date_start),
    xmax = max(date_start) + days(30),
    .groups = "drop"
  )

# Prepare incidence data
inc_df <- obs_cases %>%
  select(date_ymd, inc_C)

p_inc_line <- plot_incidence(
  obs_df             = inc_df,
  smc_df             = SMC_df,
  years_analysis     = years_analysis,
  type               = "line",
  fill_rect          = "gray30",
  rect_alpha         = 0.7,
  line_color         = "black",
  line_size          = 1,
  base_size          = 12,
  axis_text_x_size   = 15,
  axis_title_y_size  = 13,
  ylab               = "Weekly Cases (<5 yrs)"
)


p_inc_point <- plot_incidence(
  obs_df             = inc_df,
  smc_df             = SMC_df,
  years_analysis     = years_analysis,
  type               = "point",
  fill_rect          = "gray30",
  rect_alpha         = 0.7,
  point_color        = "black",
  point_size         = 2,
  base_size          = 12,
  axis_text_x_size   = 0,
  axis_text_y_size   = 13,
  axis_title_y_size  = 13,
  ylab               = "Weekly Cases (<5 yrs)"
)

######################################
## Incidence plot with SMC timeline ##
######################################
p_inc_tl <- (p_inc_line + theme_pubr() + theme(plot.title = element_text(face = "bold", size = 19)))

# Prepare your label data with Date x’s
smc_tl_df <- tibble(
  change = c(
    "SMC halted",
    "Fifth round of<br>SMC added",
    "First SMC round<br>shifted to June"
  ),
  x = as.Date(c("2019-09-01", "2021-09-15", "2023-09-01")),
  y = c(3700, 3200, 2400)
)

# Overlay with geom_richtext, nudging above the top edge
p_inc_tl <- p_inc_tl +
  geom_richtext(
    data        = smc_tl_df,
    aes(x = x, y = y, label = change),
    inherit.aes = FALSE,
    # BOX STYLE:
    #fill        = "#FFFFFFCC",   # semi-transparent white
    fill        = "#FFFFFFFF",   # semi-transparent white
    label.color = "#000000",     # dark gray border
    label.padding = grid::unit(c(6, 8, 6, 8), "pt"),
    label.r     = grid::unit(3, "pt"),
    # TEXT STYLE:
    size        = 7,             # bigger than default
    family      = "roboto",      # or "serif", "Arial", etc.
    color       = "#000000",     # dark text
    vjust       = 1.2
  )

p_inc_tl
#############
## Climate ##
#############``
p_monthly_clim <- plot_climate_monthly(rain_df = monthly_rain_df,
  temp_df = monthly_temp_df,
  plot_type = c("both"),
  separate = FALSE,
  rain_color = "steelblue",
  temp_color = "tomato",
  title_rain = "Monthly Rainfall",
  title_temp = "Monthly Temperature",
  ylab_rain = "Monthly Cumulative Rainfall (mm)",
  ylab_temp = "Average Monthly Temperature (°C)",
  title_size = 13,
  axis_text_x_size = 15,
  axis_text_y_size = 13,
  line_size = 1,
  bar_alpha = 0.8)

##################
## SMC Coverage ##
##################
SMC_df <- SMC_df %>% filter(!date_start == "2023-11-07")

smc_cov_bar <- plot_smc_coverage(df = SMC_df,
                  palette = c("steelblue", "darkorange"),
                  title = "C",
                  ylab = "Coverage",
                  xlab = "Year-Round",
                  dodge_width = 0.7,
                  errorbar_width = 0.5,
                  title_size = 14,
                  axis_text_size = 10,
                  legend_title = "Coverage Type")

##################
## Saving Plots ##
##################
# Single inc plot with smc strategy timeline
p_inc_tl <- p_inc_tl + theme(axis.title.y = element_text(size = 16, face = "bold"))

ggsave(file.path(out_dir, "p_inc_tl.png"), p_inc_tl,
       width = 14, height = 5, dpi = 1000)

# Combined plots
p_inc <- (p_inc_point + theme_pubr() + theme(axis.text.x = element_blank())) + ggtitle("A") + theme(plot.title = element_text(face = "bold", size = 19))

p_clim <- (p_monthly_clim + xlab("Year") + theme_pubr()) + ggtitle("B") + theme(plot.title = element_text(face = "bold", size = 19))

p_inc_anom_temp <- p_inc / p_clim

# SMC coverage bar plot
smc_cov_bar <- smc_cov_bar + theme_pubr() + theme(
  #plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
  axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(face = "bold", size = 19))



p_inc_clim_smc <- (p_inc_anom_temp / smc_cov_bar) + plot_layout(heights = c(2,1.25,1))
p_inc_clim_smc

ggsave(file.path(out_dir, "inc_clim_smc.pdf"), p_inc_clim_smc,
       width = 14, height = 13, dpi = 1000)


ggsave(file.path(out_dir, "inc_anom_temp.pdf"), p_inc_anom_temp,
       width = 10, height = 11, dpi = 1000)

smc_cov_bar <- smc_cov_bar + ggtitle("")

ggsave(file.path(out_dir, "smc_coverage_year_round.pdf"), smc_cov_bar,
       width = 11, height =5 , dpi = 700)
