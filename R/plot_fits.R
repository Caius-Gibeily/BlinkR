
#' @export
#' @method plot gp_model
plot.gp_model <- function(model, level = c("ind", "group", "global"), .width= c(0.50,0.8,0.99),
                          show_ci = TRUE, facet = show_ci, show_events = TRUE, dev_only = FALSE,collapse_level = FALSE,
                          show_traces = TRUE, palette = "Purples", width = 0.1, height = 0.2, size = 10, ...) {


  level <- match.arg(level)
  model_classes <- class(model)

  if (any(c("one_ind", "one_group") %in% model_classes) && level == "global") {
    stop("Global trace not present for models with one individual or group.")
  }

  tidy_model <- reconstruct_traces(model,level=level,.width=.width,dev_only=dev_only,...)

  y_elements <- c(tidy_model$y)
  if (show_ci) y_elements <- c(y_elements, tidy_model$ymin, tidy_model$ymax)
  if (show_traces && !is.null(model$traces)) y_elements <- c(y_elements, model$traces$y)
  global_y <- range(y_elements)
  global_x <- c(0,model$settings$duration)

  facet_var <- level
  group_aes <- level

  build_base_plot <- function(fit_sub, trace_sub = NULL, ev_sub = NULL) {

    p <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = global_x,
                               ylim = global_y) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time (t)", y = "Value")

    if (show_ci && !is.null(fit_sub) && nrow(fit_sub) > 0) {
      p <- p + ggdist::geom_lineribbon(
        data = fit_sub,
        ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax,
                     fill = factor(.width), group = .data[[group_aes]]),
        alpha = 0.6, linewidth = 0.5
      ) +
        ggplot2::scale_fill_brewer(palette = palette, direction = -1, name = "CI Width")
    } else if (!is.null(fit_sub) && nrow(fit_sub) > 0) {
      p <- p + ggplot2::geom_line(
        data = fit_sub,
        ggplot2::aes(x = x, y = y, group = .data[[group_aes]]),
        linewidth = 0.8
      )
    }

    if (show_traces && !is.null(trace_sub) && nrow(trace_sub) > 0) {
      p <- p + ggplot2::geom_line(
        data = trace_sub,
        ggplot2::aes(x = x, y = y_offset, group = .data[[facet_var]]),
        inherit.aes = FALSE, linewidth = 0.7, color = "black", alpha = 0.5
      )
    }

    if (show_events && !is.null(ev_sub) && nrow(ev_sub) > 0) {
      if (facet_var == "global") sub_facet_var <- "group"
      else if (facet_var == "group") sub_facet_var <- "ind"
      else sub_facet_var = "ind"
      y_aes_events <- "ind"


      event_aes <- ggplot2::aes(x = event_times,
                                y = factor(.data[[y_aes_events]]),
                                group = .data[[sub_facet_var]],
                                color = factor(.data[[y_aes_events]]),
                                fill = factor(.data[[y_aes_events]]))
      y_label <- "Events"
      y_label <- tools::toTitleCase(y_aes_events)


      p_below <- ggplot2::ggplot(data = ev_sub, mapping = event_aes) +
        ggplot2::coord_cartesian(xlim = global_x) +
        ggplot2::geom_tile(width = width, height = height) +
        ggplot2::scale_color_viridis_d(guide = "none") +
        ggplot2::scale_fill_viridis_d(guide = "none") +
        ggplot2::labs(x = "Time (t)", y = y_label) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = size))

      if (facet && facet_var == "ind") {
        p_below <- p_below + ggplot2::theme(axis.text.y = ggplot2::element_blank())
      }

      return((p / p_below) + patchwork::plot_layout(heights = c(3, 1)))
    }

    return(p)
  }

  if (!is.null(model$traces) & level=="ind") trace_data <- model$traces
  else if (!is.null(model$traces) & level=="group") trace_data <- model$group_traces
  else if (!is.null(model$traces) & level=="global") trace_data <- model$global_trace
  else trace_data <- NULL
  event_data <- if (!is.null(model$events)) model$events else NULL

  if (!facet || level == "global") {
    combined_plot <- build_base_plot(tidy_model, trace_data, event_data)
    return(combined_plot)

  } else {
    plots_list <- list()
    unique_facets <- unique(tidy_model[[facet_var]])

    for (f_id in unique_facets) {
      fit_sub <- tidy_model |> dply_filter_equal(facet_var, f_id)
      trace_sub <- trace_data |> dply_filter_equal(facet_var, f_id)
      ev_sub  <- event_data |> dply_filter_equal(facet_var, f_id)

      p_f <- build_base_plot(fit_sub, trace_sub, ev_sub)
      p_f <- p_f + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]),scales = "fixed")

      plots_list[[paste0("Facet_", f_id)]] <- p_f
    }

    combined_plot <- patchwork::wrap_plots(plots_list) +
      patchwork::plot_layout(guides = "collect")

    return(combined_plot)
  }
}

#' @export
#' @method plot gp_prior_pc
plot.gp_prior_pc <- function(model, level = c("ind", "group", "global"), .width= c(0.50,0.8,0.99),
                          show_ci = TRUE, facet = show_ci, collapse_level = FALSE, show_generated_traces = TRUE,
                          palette = "Purples", n_samples = 100, ...) {
  level <- match.arg(level)
  model_classes <- class(model)

  if (any(c("one_ind", "one_group") %in% model_classes) && level == "global") {
    stop("Global trace not present for models with one individual or group.")
  }
  tidy_model <- reconstruct_traces(model,level=level,.width = .width,...)
  prior_pc_traces <- ppc_draw_traces(model, level = level, n_samples = n_samples)

  y_elements <- c(tidy_model$y)
  if (show_ci) y_elements <- c(y_elements, tidy_model$ymin, tidy_model$ymax)

  global_y <- range(y_elements)
  global_x <- c(0,model$settings$duration)

  facet_var <- group_aes <- level

  build_base_plot <- function(fit_sub, trace_sub = NULL) {

    p <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = global_x,
                               ylim = global_y) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time (t)", y = "Value")

    if (show_ci && !is.null(fit_sub) && nrow(fit_sub) > 0) {
      p <- p + ggdist::geom_lineribbon(
        data = fit_sub,
        ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax,
                     fill = factor(.width), group = .data[[group_aes]]),
        alpha = 0.6, linewidth = 0.5
      ) +
        ggplot2::scale_fill_brewer(palette = palette, direction = -1, name = "CI Width")
    } else if (!is.null(fit_sub) && nrow(fit_sub) > 0) {
      p <- p + ggplot2::geom_line(
        data = fit_sub,
        ggplot2::aes(x = x, y = y_offset, group = .data[[group_aes]]),
        linewidth = 0.8
      )
    }

    if (show_generated_traces && !is.null(trace_sub) && nrow(trace_sub) > 0) {
      p <- p + ggplot2::geom_line(data = trace_sub,
        ggplot2::aes(x = x, y = y, group = interaction(sample,.data[[facet_var]])),
        inherit.aes = FALSE, linewidth = 0.7, color = "darkgrey", alpha = 0.9)
    }
    return(p)
  }

  trace_data <- prior_pc_traces$sampled_traces

  if (!facet || level == "global") {
    combined_plot <- build_base_plot(tidy_model, trace_data)
    return(combined_plot)

  } else {
    plots_list <- list()
    unique_facets <- unique(tidy_model[[facet_var]])

    for (f_id in unique_facets) {
      fit_sub <- tidy_model |> dply_filter_equal(facet_var, f_id)
      trace_sub <- trace_data |> dply_filter_equal(facet_var, f_id)

      p_f <- build_base_plot(fit_sub, trace_sub)
      p_f <- p_f + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]),scales = "fixed")

      plots_list[[paste0("Facet_", f_id)]] <- p_f
    }

    combined_plot <- patchwork::wrap_plots(plots_list) +
      patchwork::plot_layout(guides = "collect")

    return(combined_plot)
  }

}


#' @noRd
dply_filter_equal <- function(df, col, value) {
  if (is.null(df) || nrow(df) == 0 || !(col %in% names(df))) return(NULL)
  df[df[[col]] == value, , drop = FALSE]
}


