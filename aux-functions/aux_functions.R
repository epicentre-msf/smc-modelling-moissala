# JS_aux_functions.R

# Load necessary libraries
library(dplyr)
library(tibble)

# Calculate cases averted given baseline and intervention scenario draws
calculate_cases_averted <- function(baseline_counts, intervention_counts) {
  cases_averted <- baseline_counts - intervention_counts
  return(cases_averted)
}

# Compute effectiveness (relative reduction) draw‐wise
effectiveness_fn <- function(baseline_counts, intervention_counts) {
  effectiveness <- (baseline_counts - intervention_counts) / baseline_counts
  return(effectiveness)
}

# Prepare a tidy data.frame of raw simulation draws
prep_plot_df <- function(sim_results) {
  df_list <- lapply(names(sim_results), function(name) {
    sim <- sim_results[[name]]
    tibble(
      scenario = name,
      time     = sim$time,
      draw     = sim$counts
    )
  })
  df <- bind_rows(df_list)
  return(df)
}

# Summarize draws into mean + CI
prep_ci_df <- function(plot_df, prob = 0.95) {
  ci_df <- plot_df %>%
    group_by(scenario, time) %>%
    summarize(
      mean  = mean(draw),
      lower = quantile(draw, (1 - prob) / 2),
      upper = quantile(draw, 1 - (1 - prob) / 2),
      .groups = "drop"
    )
  return(ci_df)
}

# Wrapper returning both raw draws and CIs
prepare_plot_inputs <- function(sim_results, prob = 0.95) {
  raw <- prep_plot_df(sim_results)
  ci  <- prep_ci_df(raw, prob)
  list(raw = raw, ci = ci)
}


# Prepare raw + CI data frames (unchanged)
prepare_ppc_data <- function(summary_df, legend_levels) {
  raw <- summary_df %>%
    select(date_ymd, median, label) %>%
    rename(value = median) %>%
    mutate(label = factor(label, levels = legend_levels))
  ci  <- summary_df %>%
    select(date_ymd, lower, upper, label) %>%
    mutate(label = factor(label, levels = legend_levels))
  list(raw = raw, ci = ci)
}

# Simulation ribbon + line layer
add_simulation_layer <- function(p, ppc_data) {
  p +
    geom_ribbon(data = ppc_data$ci,
                aes(x = date_ymd, ymin = lower, ymax = upper, fill = label),
                alpha = 0.3) +
    geom_line(data = ppc_data$raw,
              aes(x = date_ymd, y = value, color = label),
              size = 0.8)
}

# Observations layer (points + dashed line)
add_observed_layer <- function(p, obs_df, obs_col = "inc_C") {
  p +
    geom_point(data = obs_df,
               aes(x = date_ymd, y = .data[[obs_col]]),
               color = "black", size = 1) +
    geom_line(data = obs_df,
              aes(x = date_ymd, y = .data[[obs_col]]),
              color = "black", linetype = "dashed")
}

# Common theme + scales, *no* legend positioning here
apply_ppc_theme <- function(p) {
  p +
    scale_color_brewer(type = "qual", palette = "Dark2", name = "Scenario") +
    scale_fill_brewer(type = "qual", palette = "Dark2", name = "Scenario") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
}


# Extract the legend grob
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  leg <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[leg]]
}

create_ppc_plot <- function(summary_df,
                            obs_df         = NULL,
                            smc_patterns   = NULL,    # your expanded_patterns list
                            smc_day        = 1,       # day of month for SMC rounds
                            xlim           = NULL,
                            ylim           = NULL,
                            xlab           = NULL,
                            ylab           = "Weekly Cases",
                            scenario_levels,
                            scenario_cols) {         # e.g. SMC_COLORS[scenario_levels]

  library(dplyr)
  library(ggplot2)
  library(purrr)

  # 1) Build exactly one box **per scenario per year**, from first SMC date to last SMC date + 30d
  rects <- NULL
  if (!is.null(smc_patterns)) {
    valid <- intersect(scenario_levels, names(smc_patterns))
    rects <- purrr::map_dfr(valid, function(lbl) {
      months_by_year <- smc_patterns[[lbl]]
      years_vec      <- as.integer(names(months_by_year))

      sched <- gen_smc_schedule(
        start_date       = min(summary_df$date_ymd),
        end_date         = max(summary_df$date_ymd),
        years            = years_vec,
        months_active    = months_by_year,
        coverage         = rep(1,12),
        smc_day_of_month = smc_day
      )

      # keep only the dates where a round occurs
      round_dates <- sched$dates[sched$SMC == 1]
      if (length(round_dates) == 0) return(NULL)

      data.frame(
        label = lbl,
        # for each year in those round_dates, one xmin/xmax row
        purrr::map_dfr(split(round_dates, year(round_dates)), function(dts) {
          data.frame(
            xmin = min(dts),
            xmax = max(dts) + 30,       # extend exactly 30 days after last round
            fill = scenario_cols[lbl],
            year = unique(year(dts))    # keep for debugging if you like
          )
        })
      )
    })
  }

  # 2) Prepare your simulation summary + CI
  if(!is.null(summary_df)){
    raw <- summary_df %>%
      select(date_ymd, median, label) %>%
      rename(value = median) %>%
      mutate(label = factor(label, levels = scenario_levels))

    ci <- summary_df %>%
      select(date_ymd, lower, upper, label) %>%
      mutate(label = factor(label, levels = scenario_levels))

    sim_raw <- filter(raw, label != "Observed")
    sim_ci  <- filter(ci,  label != "Observed")

  }

  # 3) Start the ggplot
  p <- ggplot()

  # 3a. Draw SMC boxes (if any)
  if (!is.null(rects) && nrow(rects)>0) {
    p <- p + geom_rect(
      data   = rects,
      aes(xmin = xmin, xmax = xmax),
      ymin   = -Inf, ymax   = Inf,
      fill   = rects$fill,
      alpha  = 0.35
    )
  }

  # 3b. Ribbon + line
  if (!is.null(summary_df)) {
    if(nrow(sim_ci)>0){
      p <- p + geom_ribbon(
        data  = sim_ci,
        aes(x = date_ymd, ymin = lower, ymax = upper, fill = label, group = label),
        alpha = 0.3)
    }
  }
  if (!is.null(summary_df)) {
    if(nrow(sim_raw)>0){
      p <- p + geom_line(
        data = sim_raw,
        aes(x = date_ymd, y = value, color = label, group = label),
        linewidth = 1.3)

    }
  }

  # 3c. Optional Observed
  if (!is.null(obs_df)) {
    p <- p + geom_point(
      data = obs_df,
      aes(x = date_ymd, y = inc_C, color = "Observed"),
      size = 1.2
    )
  }

  # 4) Axes
  p <- p + scale_x_date(
    breaks      = "6 months",
    date_labels = "%b\n%Y",
    limits      = xlim,
    expand      = expansion(mult = c(0.01, 0.01))
  )
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

  # 5) Colors + theme
  p <- p +
    scale_color_manual(name   = "Scenario", values = scenario_cols) +
    scale_fill_manual (values = scenario_cols, guide = "none") +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      axis.text.x     = element_text(hjust = 0.5)
    )
}

# 3. Legend‐only plot, re‐using the same vector
create_scenario_legend_plot <- function(line_size = 2, scenario_levels = scenario_levels,
                                        scenario_cols = scenario_cols) {

  # dummy data so geom_line() shows up in legend
  dummy <- expand.grid(
    label = factor(scenario_levels, levels = scenario_levels),
    x     = c(1,2)
  )
  dummy$y <- 1

  ggplot(dummy, aes(x = x, y = y, color = label)) +
    geom_line(size = line_size) +
    # **exact same** manual scale
    scale_color_manual(name   = "",
                       values = scenario_cols,
                       guide  = guide_legend(override.aes = list(size = line_size))) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 12)
    )
}



create_shared_legend <- function(labels, values, title = NULL, line_size = 2, legend_position = "bottom") {
  # Dummy data with one row per label
  dummy <- tibble(
    label = factor(labels, levels = labels),
    x     = 1,
    y     = 1  # all on same line — doesn't matter, legend only
  )

  # Build consistent override styles
  guide_overrides <- lapply(labels, function(lab) {
    if (lab == "Observed") {
      list(shape = 16, size = 3, linetype = 0)  # point only
    } else {
      list(shape = NA, size = line_size, linetype = 1)  # line only
    }
  })
  names(guide_overrides) <- labels

  guide_aes <- bind_rows(guide_overrides, .id = "label")
  rownames(guide_aes) <- guide_aes$label
  guide_aes$label <- NULL

  ggplot(dummy, aes(x = x, y = y, color = label)) +
    # Draw only one point for "Observed"
    geom_point(data = filter(dummy, label == "Observed"), size = 3) +
    geom_line(data = filter(dummy, label != "Observed"), linewidth = line_size) +
    scale_color_manual(
      name   = title,
      values = values,
      guide  = guide_legend(override.aes = guide_aes)
    ) +
    theme_void() +
    theme(
      legend.position = legend_position,
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 11)
    )
}



# Define basic monthly SMC patterns used across scenarios
pattern_4_rounds_june <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
pattern_4_rounds_july <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0)
pattern_5_rounds_june <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0)
pattern_5_rounds_july <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0)
pattern_none          <- rep(0, 12)

# Define basic scenario labels with corresponding monthly patterns
patterns <- list(
  "4 rounds (June start)" = pattern_4_rounds_june,
  "4 rounds (July start)" = pattern_4_rounds_july,
  "5 rounds (June start)" = pattern_5_rounds_june,
  "5 rounds (July start)" = pattern_5_rounds_july
)

# Define target years of simulation
target_years <- 2018:2023

# Expand each pattern to apply uniformly across all target years
expanded_patterns <- lapply(patterns, function(pat) {
  setNames(replicate(length(target_years), pat, simplify = FALSE), as.character(target_years))
})

# Define the actual (historical) SMC timing schedule by year
true_schedule <- list(
  "2018" = pattern_4_rounds_july,
  "2019" = pattern_none,
  "2020" = pattern_4_rounds_july,
  "2021" = pattern_5_rounds_july,
  "2022" = pattern_5_rounds_july,
  "2023" = pattern_5_rounds_june
)

# Add the true SMC schedule to the list of evaluated patterns
expanded_patterns[["True SMC Schedule"]] <- true_schedule

# Counterfactual 1: use 4 rounds in July for 2021 and 2022 (rest as observed)
cf_4rounds_july_21_22 <- true_schedule
cf_4rounds_july_21_22[["2021"]] <- pattern_4_rounds_july
cf_4rounds_july_21_22[["2022"]] <- pattern_4_rounds_july
expanded_patterns[["4 rounds July in 2021–2022"]] <- cf_4rounds_july_21_22

# Counterfactual 2: use 5 rounds in July for 2023 (instead of June)
cf_5rounds_july_2023 <- true_schedule
cf_5rounds_july_2023[["2023"]] <- pattern_5_rounds_july
expanded_patterns[["July Start (2023)"]] <- cf_5rounds_july_2023

smc_in_2019 <- list(
  "2018" = pattern_4_rounds_july,
  "2019" = pattern_4_rounds_july,
  "2020" = pattern_4_rounds_july,
  "2021" = pattern_5_rounds_july,
  "2022" = pattern_5_rounds_july,
  "2023" = pattern_5_rounds_june
)


expanded_patterns[["SMC in 2019"]] <- smc_in_2019

create_ppc_plot_animated <- function(summary_df,
                                     obs_df         = NULL,
                                     smc_patterns   = NULL,
                                     smc_day        = 1,
                                     xlim           = NULL,
                                     ylim           = NULL,
                                     xlab           = NULL,
                                     ylab           = "Weekly Cases",
                                     scenario_levels,
                                     scenario_cols) {
  library(dplyr)
  library(ggplot2)
  library(purrr)

  # Handle factor levels
  if(!is.null(summary_df)){
    summary_df <- summary_df %>%
    mutate(
      label = factor(label, levels = scenario_levels),
      frame = factor(frame, levels = unique(summary_df$frame))
    )
    }

  # SMC boxes — skipped in animation for simplicity (or add conditionally if needed)
  # You can reintroduce this if you animate scenarios with SMC patterns

  # Separate observed + simulated
  raw <- summary_df %>%
    select(date_ymd, value = median, label, frame)

  ci <- summary_df %>%
    select(date_ymd, lower, upper, label, frame)

  sim_raw <- filter(raw, label != "Observed")
  sim_ci  <- filter(ci,  label != "Observed")
  obs_raw <- filter(raw, label == "Observed")

  # Start ggplot
  p <- ggplot() +
    # Ribbon for CIs
    geom_ribbon(
      data = sim_ci,
      aes(x = date_ymd, ymin = lower, ymax = upper, fill = label, group = interaction(label, frame)),
      alpha = 0.3
    ) +
    # Model fit line
    geom_line(
      data = sim_raw,
      aes(x = date_ymd, y = value, color = label, group = interaction(label, frame)),
      linewidth = 1.2
    ) +
    # Observed points + dashed line
    geom_point(
      data = obs_raw,
      aes(x = date_ymd, y = value, color = label),
      size = 1.2
    ) +
    geom_line(
      data = obs_raw,
      aes(x = date_ymd, y = value, color = label, group = interaction(label, frame)),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # Axes + scale
    scale_x_date(
      breaks = "6 months",
      date_labels = "%b\n%Y",
      limits = xlim,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(ylim = ylim) +
    scale_color_manual(name = "Scenario", values = scenario_cols) +
    scale_fill_manual(values = scenario_cols, guide = "none") +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      axis.text.x     = element_text(hjust = 0.5)
    ) +
    gganimate::transition_states(
      states = frame,
      transition_length = 1,
      state_length = 2
    ) +
    labs(title = "{closest_state}")

  return(p)
}


save_time_series_plot <- function(summary_df, scenario_levels, filename, xlim, ylim,
                          smc_patterns = expanded_patterns, smc_day = 15,
                          width, height, obs_df = NULL,
                          x_axis_label_size = 12.5, y_axis_title_size = 14.5,
                          path_to_figures, save = FALSE,
                          add_x_lab = FALSE) {

  p <- create_ppc_plot(
    summary_df      = summary_df,
    obs_df          = obs_df,
    smc_patterns    = smc_patterns,
    smc_day         = smc_day,
    xlim            = xlim,
    ylim            = ylim,
    xlab            = "",
    ylab            = "Weekly Cases (<5 yrs)",
    scenario_levels = scenario_levels,
    scenario_cols   = SMC_COLORS[scenario_levels]
  ) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = x_axis_label_size),
      axis.title.y    = element_text(size = y_axis_title_size)
    )

  if (add_x_lab) {
    p <- p + scale_x_date(
      date_labels = "%b\n%Y",
      date_breaks = "2 month",
      expand      = expansion(mult = c(0.01, 0.01))
    )
  }

  ggsave(
    file.path(path_to_figures, filename),
    p,
    width  = width,
    height = height,
    dpi    = 1000
  )
}


save_case_rainfall_plot <- function(summary_df, rain_df, scenario_levels, filename,
                                    xlim, ylim, width, height, smc_patterns) {
  require(ggplot2)
  require(cowplot)
  require(lubridate)
  require(zoo)
  require(dplyr)

  # Filter rainfall and compute rolling average
  rain_roll30_weekly <- rain_df %>%
    filter(date >= xlim[1], date <= xlim[2]) %>%
    arrange(date) %>%
    mutate(
      roll30     = zoo::rollmean(rainfall, 30, fill = NA, align = "right"),
      week_start = floor_date(date, "week", week_start = 1)
    ) %>%
    group_by(week_start) %>%
    summarize(avg30 = last(na.omit(roll30)), .groups = "drop")

  # Common scale x for both panels
  common_scale_x <- scale_x_date(
    date_breaks = "3 month",
    date_labels = "%b\n%Y",
    limits      = xlim,
    expand      = c(0, 0)
  )

  # Top panel: cases
  p_cases <- create_ppc_plot(
    summary_df      = summary_df,
    obs_df          = NULL,
    smc_patterns    = smc_patterns,
    smc_day         = 15,
    xlim            = xlim,
    ylim            = ylim,
    xlab            = NULL,
    ylab            = "Weekly Cases (<5 yrs)",
    scenario_levels = scenario_levels,
    scenario_cols   = SMC_COLORS[scenario_levels]
  ) +
    common_scale_x +
    theme(
      legend.position = "none",
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.margin     = margin(0, 0, 0, 0),
      axis.title.y    = element_text(size = 11.5)
    )

  # Bottom panel: rainfall mirrored
  p_rain <- ggplot(rain_roll30_weekly, aes(week_start, -avg30)) +
    geom_col(fill = "#56B4E9", width = 7) +
    common_scale_x +
    scale_y_continuous(labels = function(x) abs(x), expand = c(0, 0)) +
    labs(x = NULL, y = "30-day Rainfall (mm)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin   = margin(3, 3, 3, 3),
      axis.text.x   = element_text(size = 12.5),
      axis.title.y  = element_text(size = 11.5)
    )

  # # Legend
  # shared_legend <- create_shared_legend(
  #   labels = scenario_levels,
  #   values = SMC_COLORS[scenario_levels]
  # )

  # Stack
  combined <- plot_grid(
    p_cases,
    p_rain,
    ncol        = 1,
    rel_heights = c(1, 0.6),
    align       = "v",
    axis        = "lr"
  )

  ggsave(
    filename,
    combined,
    width  = width,
    height = height,
    dpi    = 1000
  )
}
############################
# Load necessary libraries #
############################

#' Plot Publication-Ready MCMC Trace Panel
#'
#' Combines all trace plots into one multi-panel figure, shared x-axis label and legend.
#'
#' @param results Object with `coda_pars` or `pars` and `n_chains`.
#' @param output_dir Output directory for plots.
#' @param file_prefix Filename prefix.
#' @param save Logical; whether to save output.
#' @param save_format File format ("pdf", "png", etc.).
#' @param param_exclude Character vector of parameter names to exclude.
#' @param trace_colors Optional vector of chain colors.
#' @param plot_width,plot_height Size of final plot in inches.
#' @param base_size Base text size.
#' @param line_size Line width for trace.
#' @param line_alpha Transparency.
#' @param ncol_traces Number of columns in grid layout.
#' @param ggtheme A ggplot2 theme.
#'
#' @return Invisibly returns the patchwork plot object.
plot_mcmc_trace <- function(results,
                                  output_dir = "mcmc_plots",
                                  file_prefix = "trace",
                                  save = TRUE,
                                  save_format = "pdf",
                                  param_exclude = c("log_likelihood", "log_prior", "log_posterior"),
                                  trace_colors = NULL,
                                  plot_width = 12,
                                  plot_height = 3.5,
                                  base_size = 12,
                                  line_size = 0.4,
                                  line_alpha = 0.8,
                                  ncol_traces = 3,
                                  ggtheme = ggplot2::theme_classic(base_size = 12)) {

  if (save && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if ("coda_pars" %in% names(results)) {
    mcmc_chains <- results$coda_pars
    if (inherits(mcmc_chains, "mcmc")) {
      mcmc_chains <- coda::as.mcmc.list(list(mcmc_chains))
    }
  }

  if (!is.null(results[[1]]$pars)) {
    pars <- results[[1]]$pars
    n_chains <- results$n_chains
    n_samples <- nrow(pars) / n_chains
    param_names <- colnames(as.matrix(mcmc_chains[[1]]))

    # Remove log-quantities or internal parameters
    if (!is.null(param_exclude)) {
      param_names <- setdiff(param_names, param_exclude)
    }

    chain_list <- split(pars, rep(1:n_chains, each = n_samples))
    chains <- lapply(chain_list, function(x) {
      coda::as.mcmc(matrix(x, nrow = n_samples, ncol = ncol(pars),
                           dimnames = list(NULL, param_names)))
    })
    mcmc_chains <- coda::as.mcmc.list(chains)
  } else {
    stop("Could not identify MCMC chains in results object.")
  }

  plot_list <- list()

  for (p in param_names) {
    df <- do.call(rbind, lapply(seq_along(mcmc_chains), function(i) {
      data.frame(
        iteration = 1:n_samples,
        var1 = mcmc_chains[[i]][, p],
        chain = factor(i),
        param = p
      )
    }))

    p_plot <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = var1, color = chain)) +
      ggplot2::geom_line(size = line_size, alpha = line_alpha) +
      ggplot2::labs(x = "Iteration", y = NULL, title = p) +
      ggtheme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = element_text(face = "bold", size = 12)
      )

    if (!is.null(trace_colors)) {
      p_plot <- p_plot + ggplot2::scale_color_manual(values = trace_colors)
    }

    plot_list[[p]] <- p_plot
  }

  # Combine plots
  combined_body <- patchwork::wrap_plots(plot_list, ncol = ncol_traces) &
    ggplot2::theme(legend.position = "right")

  final_combined <- combined_body +
    patchwork::plot_layout(guides = "collect",
                           axis_titles = "collect")
  # Save or print
  if (save) {
    ggsave_fn <- paste0(file_prefix, "_trace_all.", save_format)
    ggplot2::ggsave(
      filename = file.path(output_dir, ggsave_fn),
      plot = final_combined,
      width = plot_width,
      height = plot_height * ceiling(length(plot_list) / ncol_traces)
    )
  } else {
    print(final_combined)
  }

  invisible(final_combined)
}


# 3. Legend‐only plot, re‐using the same vector
create_scenario_legend_plot <- function(line_size = 2, scenario_levels = scenario_levels,
                                        scenario_cols = scenario_cols) {

  # dummy data so geom_line() shows up in legend
  dummy <- expand.grid(
    label = factor(scenario_levels, levels = scenario_levels),
    x     = c(1,2)
  )
  dummy$y <- 1

  ggplot(dummy, aes(x = x, y = y, color = label)) +
    geom_line(size = line_size) +
    # **exact same** manual scale
    scale_color_manual(name   = "",
                       values = scenario_cols,
                       guide  = guide_legend(override.aes = list(size = line_size))) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 12)
    )
}



create_shared_legend <- function(labels, values, title = NULL, line_size = 2, legend_position = "bottom") {
  # Dummy data with one row per label
  dummy <- tibble(
    label = factor(labels, levels = labels),
    x     = 1,
    y     = 1  # all on same line — doesn't matter, legend only
  )

  # Build consistent override styles
  guide_overrides <- lapply(labels, function(lab) {
    if (lab == "Observed") {
      list(shape = 16, size = 3, linetype = 0)  # point only
    } else {
      list(shape = NA, size = line_size, linetype = 1)  # line only
    }
  })
  names(guide_overrides) <- labels

  guide_aes <- bind_rows(guide_overrides, .id = "label")
  rownames(guide_aes) <- guide_aes$label
  guide_aes$label <- NULL

  ggplot(dummy, aes(x = x, y = y, color = label)) +
    # Draw only one point for "Observed"
    geom_point(data = filter(dummy, label == "Observed"), size = 3) +
    geom_line(data = filter(dummy, label != "Observed"), linewidth = line_size) +
    scale_color_manual(
      name   = title,
      values = values,
      guide  = guide_legend(override.aes = guide_aes)
    ) +
    theme_void() +
    theme(
      legend.position = legend_position,
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 11)
    )
}



#' Build smc_day_of_month matrix from SMC start dates
#'
#' @param smc_data   A data.frame or tibble with a `date_start` column of Dates (or POSIXt).
#'                   Each row is one SMC round; months with no entries imply no SMC.
#' @param years      Integer vector of years you care about (in the order you will pass to gen_smc_schedule).
#'
#' @return An integer matrix with nrow = length(years), ncol = 12, rownames = years.
#'         Entry [i, j] is the day of month of the SMC round for year = years[i], month = j,
#'         or NA if none occurred that month.
#'
#' @examples
#' smc_data <- tibble::tibble(
#'   date_start = as.Date(c("2019-07-15","2019-08-10","2020-07-20"))
#' )
#' years <- 2019:2020
#' day_mat <- build_smc_day_matrix(smc_data, years)
#' # month matrix of active months (1 = SMC, 0 = none):
#' active_mat <- ( !is.na(day_mat) ) * 1
#'
build_smc_day_matrix <- function(smc_data, years) {
  if (!"date_start" %in% names(smc_data)) {
    stop("`smc_data` must have a `date_start` column.")
  }
  # Extract date parts
  dates <- as.Date(smc_data$date_start)
  yrs   <- as.integer(format(dates, "%Y"))
  mons  <- as.integer(format(dates, "%m"))
  dys   <- as.integer(format(dates, "%d"))

  # Initialize matrix of NAs
  mat <- matrix(
    NA_integer_,
    nrow = length(years),
    ncol = 12,
    dimnames = list(as.character(years), sprintf("%02d", 1:12))
  )

  # Fill in days where we have a date
  for (i in seq_along(dates)) {
    y <- yrs[i]; m <- mons[i]; d <- dys[i]
    row_idx <- which(years == y)
    if (length(row_idx)==1 && m >=1 && m <=12) {
      mat[row_idx, m] <- d
    }
  }

  return(mat)
}
