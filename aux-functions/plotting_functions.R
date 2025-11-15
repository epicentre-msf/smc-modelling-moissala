###########################
## MCMC Diagnostic Plots ##
###########################
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
save_mcmc_trace <- function(results,
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
                            param_order = NULL,
                            ggtheme = ggplot2::theme_classic(base_size = 12)) {

  if (save && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  latex_param_labels <- function(param) {
    label_map <- c(
      "lag_R" = "$j_2$", "lag_T" = "$j_1$", "qR" = "$q_{R}$"
      , "sigma_LT" = "$\\sigma_{LT}$", "sigma_RT" = "$\\sigma_{RT}$",
      "k1" = "$\\tau$", "size_1" = "$\\phi$", "s" = "s",
      "eff_SMC" = "$eff_{SMC}$", "b" = "b"
    )
    latex2exp::TeX(ifelse(param %in% names(label_map), label_map[[param]], param), bold = TRUE)
  }

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
      mat <- matrix(x, nrow = n_samples, ncol = ncol(pars),
             dimnames = list(NULL, param_names))
      if (!is.null(param_order)) {
        mat <- mat[,param_order]
      }
      coda::as.mcmc(mat)
    })

    mcmc_chains <- coda::as.mcmc.list(chains)
  } else {
    stop("Could not identify MCMC chains in results object.")
  }

  plot_list <- list()

  param_names <-  colnames(mcmc_chains[[1]])


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
      #ggplot2::labs(x = "Iteration", y = NULL, title = p) +
      ggplot2::labs(x = "Iteration", y = NULL, title = latex_param_labels(p)) +
      #ggtheme +
      theme_pubr() +
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

#' Publication-Ready Posterior Distribution Plots
#'
#' Plots posterior distributions for estimated parameters with optional prior overlays and true values.
#'
#' @param results_list List of MCMC outputs, each with a `$pars` matrix.
#' @param params Character vector of parameter names to include.
#' @param nrow,ncol Number of rows and columns in facet layout.
#' @param show_true Logical. If TRUE, overlay true values.
#' @param true_vals Named vector of true parameter values.
#' @param show_prior Logical. If TRUE, overlay prior distributions.
#' @param prior_n Number of samples to draw for priors.
#' @param run_labels Labels for each MCMC run.
#' @param plot_type "density" or "histogram".
#' @param title Optional plot title.
#' @param fill_alpha Fill transparency.
#' @param line_size Line thickness.
#' @param base_size Base text size for plot theme.
#' @param palette Color palette for runs.
#'
#' @return A patchwork object.
#' @export
save_posterior_prior_plot <- function(
    results_list,
    params,
    nrow = 3,
    ncol = length(params),
    show_true = FALSE,
    true_vals = NULL,
    show_prior = FALSE,
    prior_n = 1000,
    run_labels = NULL,
    plot_type = "density",
    title = NULL,
    fill_alpha = 0.5,
    line_size = 1,
    base_size = 12,
    show_legend = TRUE,
    prior_fill = "grey50",
    prior_alpha = 0.2,
    palette = RColorBrewer::brewer.pal(8, "Set2"),
    save_path = NULL,      # e.g., "figures/posteriors.pdf"
    save_width = 10,       # in inches
    save_height = 6,       # in inches
    save_dpi = 1000         # only used for raster (e.g., PNG)
) {
  # Validate input
  if (!is.list(results_list)) stop("'results_list' must be a list")
  if (is.null(run_labels)) run_labels <- paste("Run", seq_along(results_list))
  if (length(run_labels) != length(results_list)) stop("Length of 'run_labels' must match 'results_list'")
  if (!plot_type %in% c("density", "histogram")) stop("'plot_type' must be 'density' or 'histogram'")

  # Helper: LaTeX-friendly parameter labels
  latex_param_labels <- function(param) {
    label_map <- c(
      "lag_R" = "$j_2$", "lag_T" = "$j_1$", "qR" = "$q_{R}$",
      "z" = "$z$", "sigma_LT" = "$\\sigma_{LT}$", "sigma_RT" = "$\\sigma_{RT}$",
      "k1" = "$\\tau$", "fT_C" = "$\\nu_1$", "size_1" = "$\\phi$", "size_2" = "$size_{2}$",
      "eff_SMC" = "$eff_{SMC}$"
    )
    latex2exp::TeX(ifelse(param %in% names(label_map), label_map[[param]], param), bold = TRUE)
  }

  # Combine posterior draws
  df_list <- lapply(seq_along(results_list), function(i) {
    df <- as.data.frame(results_list[[i]][[1]]$pars)
    df$Run <- run_labels[i]
    df
  })
  post_df <- do.call(rbind, df_list)
  post_long <- tidyr::pivot_longer(post_df, -Run, names_to = "param", values_to = "value")
  post_long <- dplyr::filter(post_long, param %in% params)

  # Handle priors
  prior_long <- NULL
  if (show_prior) {
    priors <- return_default_priors()[params]
    prior_long <- lapply(names(priors), function(p) {
      post_vals <- post_long$value[post_long$param == p]
      x_vals <- seq(quantile(post_vals, 0.0001), quantile(post_vals, 0.9999), length.out = prior_n)
      y_vals <- exp(priors[[p]]$prior(x_vals))
      y_vals <- pmax(0, y_vals)
      y_vals <- y_vals / max(y_vals, na.rm = TRUE) * max(density(post_vals)$y, na.rm = TRUE)
      data.frame(param = p, value = x_vals, density = y_vals)
    }) |> dplyr::bind_rows()
  }

  # Plot
  color_map <- setNames(palette[seq_along(run_labels)], run_labels)
  plot_list <- lapply(params, function(p) {
    p_post <- dplyr::filter(post_long, param == p)
    g <- ggplot(p_post, aes(x = value, fill = Run, color = Run))

    g <- if (plot_type == "histogram") {
      g + geom_histogram(aes(y = ..density..), bins = 100, alpha = fill_alpha, position = "identity")
    } else {
      g + geom_density(alpha = fill_alpha, size = line_size)
    }

    # Add prior as a filled density ribbon
    if (show_prior && !is.null(prior_long)) {
      p_prior <- dplyr::filter(prior_long, param == p)
      g <- g +
        geom_ribbon(
          data = p_prior,
          aes(x = value, ymin = 0, ymax = density),
          inherit.aes = FALSE,
          fill = prior_fill,
          alpha = prior_alpha
        )
    }

    if (show_true && !is.null(true_vals) && p %in% names(true_vals)) {
      g <- g + geom_vline(xintercept = true_vals[p], linetype = "dashed", color = "black")
    }

    g + theme_pubr(base_size = base_size) +
      labs(title = latex_param_labels(p), y = "Density", x = NULL) +
      scale_fill_manual(values = color_map) +
      scale_color_manual(values = color_map) +
      theme(legend.position = (if(show_legend) "right" else "none"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            axis.title.y = element_text(face = "bold", size = 12))
  })

  final_plot <- patchwork::wrap_plots(plot_list, nrow = nrow, ncol = ncol, guides = "collect",
                        axis_titles = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, size = base_size + 2, face = "bold"),
                    axis.title.y = element_text(face = "bold", base_size))
    )

  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot = final_plot,
      width = save_width,
      height = save_height,
      dpi = save_dpi,
      units = "in"
    )
  }

  return(invisible(final_plot))
}

#' Save Correlation Plot of MCMC Samples with LaTeX Labels
#'
#' Generates and optionally saves a publication-ready correlation matrix of MCMC parameters
#' with fully customizable sampling, smoothing, ordering, themes, and LaTeX strip labels.
#'
#' @param results list: PMCMC result object (e.g., `mcmc_results_final`), containing either
#'  - `results$mcmc_run$pars`: raw matrix with metadata in first three columns, or
#'  - `results$coda_pars`: coda::mcmc object.
#' @param params character: Parameter names to include. Defaults to all columns after metadata.
#' @param param_order character: Optional ordering of `params`; unmatched names are dropped.
#' @param sample_size integer: Number of rows to sample for speed (default: 1000).
#' @param corr_method character: Correlation method for upper panels ("pearson", "spearman").
#' @param smooth_method character: Smoothing method for lower panels ("loess", "gam").
#' @param corr_size numeric: Font size for correlation coefficients (default: 4).
#' @param smooth_alpha numeric: Transparency for smoothing lines (default: 0.3).
#' @param smooth_size numeric: Line size for smoothing lines (default: 0.5).
#' @param palette character: Color vector for points (default: RColorBrewer Set2 palette).
#' @param ggtheme ggplot2 theme: Applied to the entire ggpairs object (default: theme_bw(12)).
#' @param strip_size numeric: Font size of facet-strip text (default: 7).
#' @param axis_text_angle numeric: Angle for x-axis strip labels (default: 45).
#' @param axis_breaks integer: Number of axis breaks via scales::pretty_breaks (default: 2).
#' @param title character: Optional overall plot title.
#' @param save logical: If TRUE, saves the plot to disk (default: FALSE).
#' @param output_dir character: Directory to save plots (default: "mcmc_plots").
#' @param file_prefix character: Prefix for saved filename (appends "_corr").
#' @param save_format character: File format for saving ("pdf", "png", etc., default: "pdf").
#' @param width numeric: Width in inches for saved plot (default: 7).
#' @param height numeric: Height in inches for saved plot (default: 7).
#' @param dpi integer: Resolution for raster formats (default: 300).
#' @return Invisibly returns the `GGally::ggpairs` object.
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 theme_bw theme ggtitle element_text scale_x_continuous scale_y_continuous ggsave labs
#' @importFrom scales pretty_breaks
#' @importFrom latex2exp TeX
#' @export
save_corr_plot <- function(
    results,
    params        = NULL,
    sample_size   = 1000,
    corr_method   = "pearson",
    smooth_method = "loess",
    hist_fill     = "grey80",
    point_alpha   = 0.3,
    smooth_size   = 0.5,
    cor_size      = 6,
    strip_size    = 12,
    ggtheme       = theme_bw(base_size = 12),
    save          = FALSE,
    output_dir    = "mcmc_plots",
    file_prefix   = "corr_manual",
    save_format   = "png",
    width         = 8,
    height        = 8,
    dpi           = 300
) {
  # 1) Extract & sample
  df_raw <- if (!is.null(results$mcmc_run$pars)) {
    as.data.frame(results$mcmc_run$pars)
  } else if (!is.null(results$coda_pars)) {
    as.data.frame(as.matrix(results$coda_pars))
  } else stop("No parameter data found.")
  if (is.null(params)) params <- colnames(df_raw)
  df <- df_raw[, params, drop = FALSE]
  set.seed(123)
  df <- df %>% slice_sample(n = min(nrow(df), sample_size))

  # 2) LaTeX label map & helper

  label_map <- c(
    lag_R    = "$j_2$",
    lag_T    = "$j_1$",
    qR       = "$q_R$",
    b = "$b$",
    s = "$s$",
    z        = "$z$",
    sigma_LT = "$\\sigma_{LT}$",
    sigma_RT = "$\\sigma_{RT}$",
    k1       = "$\\tau$",
    fT_C     = "$\\nu_1$",
    size_1   = "$\\phi$",
    size_2   = "$size_2$",
    eff_SMC  = "$eff_{SMC}$"
  )

  # 0) helper: safe mk_label
  mk_label <- function(x) {
    if (! x %in% names(label_map)) {
      # fall back to raw name
      return(latex2exp::TeX(x, bold = TRUE))
    }
    latex2exp::TeX(label_map[[x]], bold = TRUE)
  }

  blank <- ggplot() + theme_void()
  N     <- length(params)
  grid  <- vector("list", (N+1)^2)

  for(i in 0:N) {
    for(j in 0:N) {
      idx <- i * (N + 1) + j + 1

      # A) only top‐right is blank
      if (i == 0 && j == N) {
        grid[[idx]] <- blank

        # B) top row (columns 1..N–1): column labels
      } else if (i == 0 && j >= 1 && j < N) {
        grid[[idx]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.1,
                   label = mk_label(params[j]), size = strip_size,
                   hjust = 0.5, vjust = 0) +
          theme_void()

        # C) right column (rows 1..N): row labels
      } else if (j == N && i >= 1) {
        grid[[idx]] <- ggplot() +
          annotate("text", x = 0.9, y = 0.5,
                   label = mk_label(params[i]), angle = -90,
                   size = strip_size,
                   hjust = 1, vjust = 0.5) +
          theme_void()

        # D) interior panels
      } else if (i >= 1 && j >= 1) {
        if (j > i) {
          # upper triangle: correlation text
          cors <- cor(df[[params[j]]], df[[params[i]]], method = corr_method)
          lab  <- sprintf("%.2f", cors)
          grid[[idx]] <- ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = lab, size = cor_size) +
            theme_void()
        } else if (j < i) {
          # lower triangle: scatter + smooth
          grid[[idx]] <- ggplot(df, aes_string(x = params[j], y = params[i])) +
            geom_point(alpha = point_alpha) +
            geom_smooth(method = smooth_method, size = smooth_size, se = FALSE) +
            theme_classic() +
            #ggtheme +
            theme(
              axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank()
            )
        } else {
          # diagonal: histogram
          grid[[idx]] <- ggplot(df, aes_string(x = params[i])) +
            geom_histogram(fill = hist_fill, color = NA) +
            theme_classic() +
            #ggtheme +
            theme(
              axis.title = element_blank(),
              axis.text  = element_blank(),
              axis.ticks = element_blank()
            )
        }
      } else {
        # all other corner cells (including top-left) get blank
        grid[[idx]] <- blank
      }
    }
  }

  final <- wrap_plots(grid, ncol = N+1)

  # 5) Save if requested
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    fn <- file.path(output_dir, paste0(file_prefix, ".", save_format))
    ggsave(fn, final, width = width, height = height, dpi = dpi)
  }

  invisible(final)
}


######################
## SMC Impact Plots ##
######################
create_time_series_plot <- function(
    summary_df,
    obs_df                 = NULL,
    smc_day_of_month_list,       # named list of matrices, or a single matrix
    xlim                   = NULL,
    ylim                   = NULL,
    xlab                   = NULL,
    ylab                   = "Weekly Cases",
    scenario_levels,
    scenario_cols,
    scenario_linetypes = NULL,
    legend_order = NULL,
    x_axis_label_size      = 12.5,
    y_axis_title_size      = 14.5,
    axis_text_size         = 10,
    legend_position        = "top",
    legend_text_size       = 12,
    add_x_lab              = FALSE,
    plot_title             = "",
    title_size             = 12,
    # NEW:
    smc_rect_shape         = c("rectangle","bar"),
    smc_bar_height         = 0.07,  # height as fraction of y-range (for "bar")
    smc_bar_gap            = 0.01,  # gap as fraction of y-range (for "bar")
    smc_alpha              = 0.35   # transparency for SMC shading
) {
  library(dplyr); library(ggplot2); library(purrr); library(lubridate); library(rlang)

  smc_rect_shape <- match.arg(smc_rect_shape)

  if (!is.null(scenario_linetypes)) {
    if (!all(scenario_levels %in% names(scenario_linetypes))) {
      stop("scenario_linetypes must be a named vector including all scenario_levels.")
    }
  } else {
    scenario_linetypes <- setNames(rep("solid", length(scenario_levels)), scenario_levels)
  }

  # 1) Normalize input: ensure we have a named list of day-matrices
  rects <- NULL
  if (!is.null(smc_day_of_month_list)) {
    if (is.matrix(smc_day_of_month_list)) {
      # same matrix for all scenarios
      smc_day_of_month_list <- set_names(
        replicate(length(scenario_levels), smc_day_of_month_list, simplify = FALSE),
        scenario_levels
      )
    }
    # check names
    missing <- setdiff(scenario_levels, names(smc_day_of_month_list))
    if (length(missing) > 0) {
      stop("Missing smc_day_of_month matrices for scenarios: ", paste(missing, collapse = ", "))
    }

    # 2) Build one rect per SMC stretch (grouped within a year), per scenario
    rects <- map_dfr(scenario_levels, function(lbl) {
      day_mat <- smc_day_of_month_list[[lbl]]
      years   <- as.integer(rownames(day_mat))
      months_active <- (!is.na(day_mat)) * 1L

      sched <- gen_smc_schedule(
        start_date        = min(summary_df$date_ymd),
        end_date          = max(summary_df$date_ymd),
        years             = years,
        months_active     = months_active,
        coverage          = rep(1, 12),
        smc_day_of_month  = day_mat
      )

      rounds <- sched %>%
        filter(SMC == 1) %>%
        mutate(year = year(dates)) %>%
        group_by(year) %>%
        summarise(
          label = lbl,
          xmin  = min(dates),
          xmax  = max(dates) + days(30),
          .groups = "drop"
        ) %>%
        mutate(fill = scenario_cols[lbl])

      if (nrow(rounds) == 0) return(NULL)
      rounds
    })

    # ensure rects$label is a factor in scenario order for stable stacking
    if (!is.null(rects) && nrow(rects) > 0) {
      rects$label <- factor(rects$label, levels = scenario_levels)
    }
  }

  # 3) Prepare simulation lines + ribbons
  sim_raw <- summary_df %>%
    filter(label != "Observed") %>%
    transmute(date_ymd = date_ymd, value = median,
              label = factor(label, levels = scenario_levels))

  sim_ci <- summary_df %>%
    filter(label != "Observed") %>%
    transmute(date_ymd = date_ymd, lower = lower, upper = upper,
              label = factor(label, levels = scenario_levels))

  # 4) Build the plot
  p <- ggplot()

  # --- Determine y-range early so we can place top bars properly ---
  # If user supplied ylim, use that; else compute from data
  if (is.null(ylim)) {
    y_candidates <- c(
      if (nrow(sim_ci) > 0) c(sim_ci$lower, sim_ci$upper) else numeric(0),
      if (!is.null(obs_df)) obs_df$inc_C else numeric(0),
      if (nrow(sim_raw) > 0) sim_raw$value else numeric(0)
    )
    y_candidates <- y_candidates[is.finite(y_candidates)]
    if (length(y_candidates) == 0) {
      y_min <- 0; y_max <- 1
    } else {
      y_min <- min(y_candidates, na.rm = TRUE)
      y_max <- max(y_candidates, na.rm = TRUE)
      if (y_min == y_max) { y_min <- 0; y_max <- y_max + 1 }
    }
  } else {
    y_min <- ylim[1]; y_max <- ylim[2]
  }
  y_range <- y_max - y_min
  if (y_range <= 0) y_range <- 1

  # 4a) SMC bands/bars for each scenario
  if (!is.null(rects) && nrow(rects) > 0) {

    if (smc_rect_shape == "rectangle") {
      # Full-height background rectangles (original behavior)
      p <- p + geom_rect(
        data   = rects,
        aes(xmin = xmin, xmax = xmax),
        ymin   = -Inf, ymax = Inf,
        fill   = rects$fill,
        alpha  = smc_alpha,
        inherit.aes = FALSE
      )

    } else if (smc_rect_shape == "bar") {
      # Short bars at top, stacked by scenario order to avoid overlap
      # Compute per-scenario vertical slots
      uniq_labels <- levels(rects$label)
      uniq_labels <- uniq_labels[uniq_labels %in% unique(as.character(rects$label))]

      n_slots <- length(uniq_labels)
      h <- smc_bar_height * y_range
      g <- smc_bar_gap    * y_range

      slot_map <- setNames(seq_along(uniq_labels), uniq_labels)

      rects_bar <- rects %>%
        mutate(slot = slot_map[as.character(label)]) %>%
        mutate(
          ymax = y_max - (slot - 1) * (h + g),
          ymin = pmax(y_max - slot * (h + g) + g, y_min)  # keep inside range
        )

      p <- p + geom_rect(
        data   = rects_bar,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill   = rects_bar$fill,
        alpha  = smc_alpha,
        inherit.aes = FALSE
      )
    }
  }

  # 4b) CI ribbon + median line
  if (nrow(sim_ci) > 0) {
    p <- p + geom_ribbon(
      data = sim_ci,
      aes(x = date_ymd, ymin = lower, ymax = upper,
          fill = label, group = label),
      alpha = 0.3
    )
  }
  if (nrow(sim_raw) > 0) {
    p <- p + geom_line(
      data = sim_raw,
      aes(x = date_ymd, y = value,
          color = label, linetype = label, group = label),
      linewidth = 1.3
    )
  }

  # 4c) Observed points
  if (!is.null(obs_df)) {
    p <- p + geom_point(
      data = obs_df,
      aes(x = date_ymd, y = inc_C, color = "Observed"),
      size = 1.2
    )
  }

  # 5) Axes & limits
  p <- p + scale_x_date(
    breaks      = "6 months",
    date_labels = "%b\n%Y",
    limits      = xlim,
    expand      = expansion(mult = c(0.01, 0.01))
  )
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  # 6) Scales & theme
  p <- p +
    guides(color = guide_legend(override.aes = list(linetype = scenario_linetypes))) +
    scale_color_manual(name = "Scenario", values = scenario_cols, breaks = legend_order) +
    scale_linetype_manual(name = "Scenario", values = scenario_linetypes, breaks = legend_order) +
    scale_fill_manual(values = scenario_cols, guide = "none") +
    labs(x = xlab, y = ylab, title = plot_title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = title_size, hjust = 0.5),
      legend.position  = legend_position,
      legend.title     = element_blank(),
      legend.text      = element_text(size = legend_text_size, face = "bold"),
      axis.text        = element_text(size = axis_text_size),
      axis.title.x     = element_text(size = x_axis_label_size),
      axis.title.y     = element_text(size = y_axis_title_size),
      axis.text.x      = element_text(hjust = 0.5)
    )

  if (add_x_lab) {
    p <- p + scale_x_date(
      date_labels = "%b\n%Y",
      date_breaks = "2 months",
      expand      = expansion(mult = c(0.01, 0.01))
    )
  }

  return(p)
}


#######################
## Descriptive Plots ##
#######################
plot_anom <- function(met_df, ylab = "Rainfall Anomaly (Z-score)",
                      fill_color = "steelblue",
                      line_color = "black",
                      title = "Standardized Rainfall Anomaly",
                      title_size = 13, axis_text_size = 10,
                      line_size = 1) {

  ggplot(met_df, aes(x = dates, y = anom)) +
    geom_col(fill = fill_color) +
    geom_hline(yintercept = 0, linetype = "dashed", color = line_color) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(x = NULL, y = ylab, title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = title_size - 1),
      axis.text = element_text(size = axis_text_size)
    )
}


plot_temp <- function(met_df, ylab = "Temperature (°C)",
                      line_color = "tomato",
                      title = "Temperature Over Time",
                      title_size = 13, axis_text_size = 10,
                      line_size = 1) {
  ggplot(met_df, aes(x = dates, y = temp)) +
    geom_line(color = line_color, linewidth = line_size) +
    labs(x = NULL, y = ylab, title = title) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = title_size - 1),
      axis.text = element_text(size = axis_text_size)
    )
}

plot_incidence <- function(obs_df,
                           smc_df           = NULL,
                           years_analysis   = NULL,
                           type             = "line",  # "point", "line", "both", "bar"
                           fill_rect        = "gray30",
                           rect_alpha       = 0.7,
                           line_color       = "black",
                           line_size        = 1,
                           point_color      = "black",
                           point_size       = 2,
                           bar_fill         = "black",
                           bar_alpha        = 0.8,
                           ylab             = "Weekly Cases (<5 yrs)",
                           x_breaks         = "1 year",
                           x_labels         = "%Y",
                           base_size        = 12,
                           axis_text_x_size = 12,
                           axis_text_y_size = 12,
                           axis_title_y_size = 14,
                           plot_title       = NULL) {
  # Prepare SMC rectangles
  smc_rects <- NULL
  if (!is.null(smc_df) && !is.null(years_analysis)) {
    smc_rects <- smc_df %>%
      filter(year(date_start) %in% years_analysis) %>%
      mutate(year = year(date_start)) %>%
      group_by(year) %>%
      dplyr::summarize(
        xmin = min(date_start),
        xmax = max(date_start) + days(30),
        .groups = "drop"
      )
  }

  # Begin ggplot
  p <- ggplot()

  # Add rectangles if present
  if (!is.null(smc_rects)) {
    p <- p + geom_rect(data = smc_rects,
                       aes(xmin = xmin, xmax = xmax),
                       ymin = -Inf, ymax = Inf,
                       fill = fill_rect, alpha = rect_alpha)
  }

  # Add incidence layer(s)
  if (type %in% c("line", "both")) {
    p <- p + geom_line(data = obs_df,
                       aes(x = date_ymd, y = inc_C),
                       color = line_color,
                       linewidth = line_size)
  }
  if (type %in% c("point", "both")) {
    p <- p + geom_point(data = obs_df,
                        aes(x = date_ymd, y = inc_C),
                        color = point_color,
                        size = point_size)
  }
  if (type == "bar") {
    p <- p + geom_col(data = obs_df,
                      aes(x = date_ymd, y = inc_C),
                      fill = bar_fill,
                      alpha = bar_alpha)
  }

  # Final theming
  p <- p +
    scale_x_date(date_breaks = x_breaks, date_labels = x_labels) +
    labs(x = NULL, y = ylab, title = plot_title) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      axis.text.x     = element_text(size = axis_text_x_size),
      axis.text.y     = element_text(size = axis_text_y_size),
      axis.title.y    = element_text(size = axis_title_y_size),
      legend.title    = element_blank()
    )

  return(p)
}
plot_climate_monthly <- function(rain_df,
                                 temp_df,
                                 plot_type = c("both", "rain", "temp"),
                                 separate = FALSE,
                                 rain_color = "steelblue",
                                 temp_color = "tomato",
                                 title_rain = "Monthly Rainfall",
                                 title_temp = "Monthly Temperature",
                                 title_combined = "",
                                 ylab_rain = "Monthly Cumulative Rainfall (mm)",
                                 ylab_temp = "Average Monthly Temperature (°C)",
                                 title_size = 13,
                                 axis_text_x_size = 15,
                                 axis_text_y_size = 15,
                                 line_size = 1,
                                 bar_alpha = 0.8) {

  plot_type <- match.arg(plot_type)

  if (plot_type %in% c("rain", "both")) {
    p_rain <- ggplot(rain_df, aes(x = dates, y = rainfall)) +
      geom_col(fill = rain_color, alpha = bar_alpha) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(x = NULL, y = ylab_rain, title = title_rain) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = title_size - 1),
        axis.text.x = element_text(size = axis_text_x_size),
        axis.text.y = element_text(size = axis_text_y_size)
      )
  }

  if (plot_type %in% c("temp", "both")) {
    p_temp <- ggplot(temp_df, aes(x = dates, y = temp)) +
      geom_point(color = temp_color, size = 2, shape = 15) +
      geom_line(color = temp_color, linewidth = line_size) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      labs(x = NULL, y = ylab_temp, title = title_temp) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = title_size - 1),
        axis.text.x = element_text(size = axis_text_x_size),
        axis.text.y = element_text(size = axis_text_y_size)
      )
  }

  if (plot_type == "both" && !separate) {
    # Join rainfall and temperature data

    df_combined <- full_join(rain_df, temp_df, by = "dates")

    p_dual <- ggplot(df_combined, aes(x = dates)) +
      geom_col(aes(y = rainfall), fill = rain_color, alpha = bar_alpha) +
      geom_line(aes(y = scales::rescale(temp, to = range(rain_df$rainfall, na.rm = TRUE))),
                color = temp_color, linewidth = line_size) +
      geom_point(aes(y = scales::rescale(temp, to = range(rain_df$rainfall, na.rm = TRUE))),
                color = temp_color, size = 2, shape = 15) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      scale_y_continuous(
        name = ylab_rain,
        sec.axis = sec_axis(~scales::rescale(., from = range(rain_df$rainfall, na.rm = TRUE), to = range(temp_df$temp, na.rm = TRUE)),
                            name = ylab_temp)
      ) +
      labs(x = NULL, title = title_combined) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        axis.title.y.left  = element_text(size = title_size - 1, color = rain_color),
        axis.title.y.right = element_text(size = title_size - 1, color = temp_color),
        axis.text.y = element_text(size = axis_text_y_size),
        axis.text.x = element_text(size = axis_text_x_size),
        axis.text.y.right = element_text(color = temp_color),
        axis.text.y.left = element_text(color = rain_color)
      )

    return(p_dual)
  }

  if (plot_type == "rain") return(p_rain)
  if (plot_type == "temp") return(p_temp)
  if (plot_type == "both" && separate) return(p_rain / p_temp)
}

plot_smc_coverage <- function(df,
                              palette = c("steelblue", "darkorange"),
                              title = "SMC Coverage by Round and Year",
                              ylab = "Coverage",
                              xlab = "Year-Round",
                              dodge_width = 0.7,
                              errorbar_width = 0.25,
                              title_size = 14,
                              axis_text_size = 10,
                              legend_title = "Coverage Type") {

  # Ensure year-round identifier
  df <- df %>%
    mutate(year = year(date_start),
           round = smc_round,
           year_round = paste0(year, "-R", round)) %>%
    select(year_round, smc_couv_card, smc_couv_card_lower, smc_couv_card_upper,
           smc_couv_tot, smc_couv_tot_lower, smc_couv_tot_upper)

  # Pivot to long format
  df_long <- df %>%
    pivot_longer(
      cols = -year_round,
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(type = case_when(
      grepl("^smc_couv_card(_.*)?$", metric) ~ "Card only",
      grepl("^smc_couv_tot(_.*)?$", metric) ~ "Card + verbal"
    ),
    bound = case_when(
      grepl("_lower$", metric) ~ "lower",
      grepl("_upper$", metric) ~ "upper",
      TRUE ~ "estimate"
    )) %>%
    select(-metric) %>%
    pivot_wider(names_from = bound, values_from = value)

  # Plot
  ggbarplot(df_long,
            x = "year_round", y = "estimate",
            add = "none",
            fill = "type",
            palette = palette,
            position = position_dodge(dodge_width),
            xlab = xlab, ylab = ylab,
            title = title,
            legend.title = legend_title) +
    geom_errorbar(aes(ymin = lower, ymax = upper, group = type),
                  width = errorbar_width,
                  position = position_dodge(dodge_width)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.text = element_text(size = axis_text_size),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]

  plots.grobs <- lapply(list(...), ggplotGrob)

  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })

  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })

  plots.grobs.eq.widths.aligned
}

############################
## Optimal strategy plots ##
############################
#' Plot SMC cases averted barplot
#'
#' @param df A tibble containing at least:
#'   - pattern (chr): strategy label
#'   - cases_averted (num)
#'   - smc_month_rounds (chr)
#'   - smc_start_day (num, e.g. 1 vs 15)
#' @param orientation One of "vertical" or "horizontal"
#' @return A ggplot object
plot_smc_barplot <- function(df,
                             fill_colors,
                             orientation = c("vertical","horizontal"),
                             plot_title = "Cases Averted by SMC Timing Strategy",
                             x_axis_lab = NULL,
                             y_axis_lab = "Cases Averted") {
  orientation <- match.arg(orientation)

  df2 <- df %>%
    mutate(
      strategy = fct_reorder(strategy, cases_averted),
      start_day = factor(smc_start_day, levels = c(15,1),
                         labels = c("15th","1st"))
    )

  p <- ggplot(df2, aes(
    x = strategy ,
    y = cases_averted,
    fill = smc_month_rounds,
    pattern = start_day
  )) +
    #geom_col_pattern(aes(fill = smc_month_rounds, pattern_density = start_day)) +
    geom_col_pattern(aes(fill = smc_month_rounds, pattern_density = start_day),
                     colour = "black",
                     #pattern = "stripe",
                     pattern_fill    = "black",
                     #pattern_angle   = 45,
                     pattern_density = 0.3,
                     #pattern_spacing = 0.05,
                     pattern_key_scale_factor = 0.2,
    ) +
    # optional error bars
    { if(all(c("lower","upper") %in% names(df2)))
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)
      else NULL } +
    scale_pattern_manual(
      name   = "SMC Start Day",
      values = c("15th" = "none",   # solid
                 "1st"  = "stripe",
                 "4R July" = "none")) +
    scale_fill_manual(name = "Strategy",
                      values = fill_colors,
    ) +
    guides(
      # Force the fill legend keys to be solid
      fill = guide_legend(
        override.aes = list(
          pattern = "none"
        )
      ),
      # And for the pattern legend, keep your hatch override
      pattern = guide_legend(
        override.aes = list(
          fill    = "white",
          colour  = "black",
          pattern = c("none", "stripe"),
          pattern_fill    = "black",
          pattern_angle   = 45,
          pattern_density = 0.3,
          pattern_spacing = 0.05
        )
      )
    ) +
    labs(x = x_axis_lab, y = y_axis_lab, title = plot_title) +
    theme_pubr(base_size = 14) +
    theme(
      legend.position = "right",
      #legend.key.size = unit(1, "cm"),
      legend.text     = element_text(size = 12),
      legend.title = element_text(face = "bold"),
      axis.text.x     = element_text(
        #angle = if(orientation=="vertical") 0 else 45,
        hjust = 1
      )
    )

  if (orientation == "horizontal") {
    p <- p + coord_flip()
  }

  return(p)
}
