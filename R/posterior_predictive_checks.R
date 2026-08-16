#' @export
plot_posterior_pc <- function(model, n_samples = 1000, palette = "Blues",
                              .width = c(0.5,0.8,0.99), bw_kde = 0.5, bw_rate = 5, scale_factor = 1, gridsize = 500) {
  ppc_events <- ppc_draw_events(model,n_samples = n_samples)

  if (!("gp_model") %in% class(model)) stop("Please use a fitted renewr model")

  # Global blink rate
  p1 <- ppc_get_eventrate(model = model,scale_factor=scale_factor)
  # Interevent distribution
  p2 <- ppc_get_interevent_dist(ppc_events = ppc_events,
                                    model = model)
  # Instantaneous blink rate
  p3 <- ppc_get_inst_eventrate(ppc_events = ppc_events,
                                   model = model)

  .ggplot_successive(p1,p2,p3)
}

#' @export
ppc_draw_events <- function(model, n_samples) {

  post_data <- ppc_draw_traces(model, n_samples = n_samples)
  samples <- post_data$sampled_traces
  params <- post_data$survival_params

  family <- if (!is.null(params$k) & !is.null(params$shape)) {
    "gengamma"
  } else if (!is.null(params$k)) {
    "gamma"
  } else if (!is.null(params$shape)) {
    "weibull"
  } else if (!is.null(params$sigma)) {
    "lognormal"
  } else {
    "exponential"
  }

  ppc_events <- samples |> group_split(sample) |>
    imap(\(samp,i) {
      ppc <- do.call(simulate_events,
                         c(list(trace_data = samp, family = family, lerp = NULL),
                                lapply(params, function(mat) mat[i])))
      ppc$events}
      ) |>
    list_rbind(names_to = "sample")

  return(ppc_events)
}

#' @export
ppc_get_eventrate <- function(model,ppc_events=NULL,scale_factor=1,.width=c(0.5,0.8,0.99),palette = "Purples",return_plot=TRUE, n_samples=1000) {
  if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples)
  }

  br_global <- model$events |>
    dplyr::group_by(ind) |>
    dplyr::reframe(ind = unique(ind),
                   group = unique(group),
                   br_global = n() / model$settings$duration * scale_factor)

  br_global_ppc <- ppc_events |>
    dplyr::group_by(ind,sample) |>
    dplyr::reframe(sample = unique(sample),
                   ind = unique(ind),
                   group = unique(group),
                   br_global = n() / model$settings$duration * scale_factor)

  if (return_plot) {
    p <- ggplot2::ggplot() +
      ggdist::stat_interval(data = br_global_ppc, ggplot2::aes(x = factor(ind), y = br_global),
                            .width = .width, size = 6) +
      ggplot2::geom_segment(data = br_global,
                            ggplot2::aes(x = as.numeric(factor(ind)) - 0.2,
                                         xend = as.numeric(factor(ind)) + 0.2,
                                         y = br_global,yend = br_global),
                            color = "#014636",size = 1.2) +
      ggplot2::scale_color_brewer(palette = palette, name = "CI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Individual",
                    y = "Event rate (blinks/min)",
                    color = "Interval")
    return(p)
  } else {
  return(list(br_global=br_global,
              br_global_ppc = br_global_ppc))
  }
}

#' @export
ppc_get_interevent_dist <- function(model, ppc_events = NULL, n_samples = 1000, .width=c(0.5,0.8,0.99), bw_kde = 1.2, return_plot = TRUE) {
  if (is.null(ppc_events) & is.null(model)) {
    stop("Please provide either posterior generated events or a fitted Renewr model")
  } else if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples)
  }

  dt_high <- stats::quantile(ppc_events$dt, 0.99)
  iei_dist <- ppc_events |>
    dplyr::group_by(ind, sample) |>
    dplyr::reframe(
      x = density(dt, n = 512, from = 0, to = quantile(dt,0.9), bw = bw_kde)$x,
      y = density(dt, n = 512, from = 0, to = quantile(dt,0.9), bw = bw_kde)$y) |>
    dplyr::group_by(ind, x) |>
    ggdist::mean_qi(y, .width = .width) |>
    dplyr::ungroup()

  if (return_plot) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = model$events, ggplot2::aes(x = dt, y = ggplot2::after_stat(density),group=ind),
        fill = "darkgrey", color = "white", alpha = 1,
        binwidth = .fd_binwidth(model$events$dt)) +
      ggdist::geom_lineribbon(data = iei_dist,
                              ggplot2::aes(x=x,y=y,
                                           ymin = .lower, ymax = .upper,
                                           fill = forcats::fct_rev(ordered(.width)), group = interaction(ind, .width)),
                              alpha = 0.6) +
      ggplot2::facet_wrap(~ ind, scale="free_x") +
      ggplot2::scale_fill_brewer(palette = palette, name = "CI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Inter-event Interval (IBI  (s))", y = "Density")
    return(p)
  }
  else {
    return(iei_dist)
  }
}

#' @export
ppc_get_inst_eventrate <- function(ppc_events = NULL,model,
                                       bw_rate = 5,gridsize = 500, n_samples = 1000, .width = c(0.5,0.8,0.99), return_plot = TRUE, palette = "Purples",show_events=TRUE) {

  if (is.null(ppc_events) & is.null(model)) {
    stop("Please provide either posterior generated events or a fitted Renewr model")
  } else if (is.null(ppc_events)) {
    ppc_events <- ppc_draw_events(model,n_samples = n_samples)
  }

  time_grid <- seq(0,model$settings$duration,
                   length.out = gridsize)

  inst_rate_ppc <- ppc_events |>
    dplyr::group_by(ind, sample) |>
    dplyr::summarise(
      rate = list(sapply(time_grid, function(t)
        sum(dnorm(t - event_times, sd = bw_rate)))), .groups = "drop") |>
    dplyr::mutate(x = list(time_grid)) |>
    tidyr::unnest(c(x, rate)) |>
    dplyr::group_by(ind, x) |>
    ggdist::mean_qi(rate, .width = .width) |>
    dplyr::ungroup()

  inst_rate <- model$events |>
    dplyr::group_by(ind) |>
    dplyr::summarise(rate = list(sapply(time_grid, function(t)
      sum(dnorm(t - event_times, sd = bw_rate)))), .groups = "drop") |>
    dplyr::mutate(x = list(time_grid)) |>
    tidyr::unnest(c(x,rate))

  if (return_plot) {
    p_inst_event_rate <- ggplot2::ggplot(inst_rate_ppc) +
      ggdist::geom_lineribbon(aes(x = x, y = rate, ymin = .lower, ymax = .upper,
                                  fill = forcats::fct_rev(ordered(.width)), group = interaction(ind, .width)),
                              alpha = 0.8, linewidth = 0.6) +

      ggplot2::scale_fill_brewer(palette = palette, name = "CI Width") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time", y = "Instantaneous event rate (s)",
                    fill = "Credible Interval") +
      ggplot2::geom_line(data=inst_rate,
                         ggplot2::aes(x = x, y = rate, group=ind),
                         linewidth=1,color = "black", alpha = 0.5) +
      ggplot2::facet_wrap(~ind, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::coord_cartesian(xlim = c(0, model$settings$duration))
    if (!show_events) return(p_inst_event_rate)

      p_ob_events <- ggplot(model$events, ggplot2::aes(x = event_times, y = factor(1), group = ind)) +
      ggplot2::geom_tile(width = 0.3, height = 0.5) +
      ggplot2::labs(x = "Time", y = "Ind") +
      ggplot2::facet_wrap(~ind) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank()) +
      ggplot2::coord_cartesian(xlim = c(0, model$settings$duration))

    p_ppc_events <- ggplot2::ggplot(ppc_events, ggplot2::aes(x = event_times, y = factor(sample),
                                                                     group = interaction(ind,sample))) +
      ggplot2::geom_tile(width = 0.5, height = 1) +
      ggplot2::labs(x = "Time", y = "Ind") +
      ggplot2::facet_wrap(~ind) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank()) +
      ggplot2::coord_cartesian(xlim = c(0, model$settings$duration))

    p <- cowplot::plot_grid(p_inst_event_rate, p_ob_events,
                             p_ppc_events , align = "v",
                             ncol = 1, rel_heights = c(1.5,0.4,0.8),
                             rel_widths = c(1,0.5,1))
    return(p)
  } else {
  return(list(inst_rate_ppc = inst_rate_ppc, inst_rate = inst_rate))
  }
}

#' @noRd
.fd_binwidth <- function(x) {
  2 * IQR(x, na.rm = TRUE) / (length(x)^(1/3))
}

#' @noRd
.ggplot_successive <- function(...) {
  input_list <- list(...)

  if (length(input_list) == 1 && is.list(input_list[[1]])) {
    plots <- input_list[[1]]
  } else {
    plots <- input_list
  }

  for (i in seq_along(plots)) {
    print(plots[[i]])

    if (i < length(plots)) {
      res <- readline(prompt = "Hit <Return> to see next plot (or 'q' to quit): ")
      if (tolower(trimws(res)) == "q") break
    }
  }
}


