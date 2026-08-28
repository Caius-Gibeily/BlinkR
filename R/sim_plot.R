#' Plot Simulation Traces
#'
#' @param x An object of class \code{sim_traces}.
#' @param level Character vector specifying which levels to plot.
#'   Options include \code{"individual"}, \code{"group"}, or \code{"global"}. Defaults to all three.
#' @param facet Logical. If \code{TRUE}, splits plots up by group or individual.
#' @param separate_pages Logical. If \code{TRUE} and \code{facet = TRUE}, outputs a layout
#'   across multiple pages via \code{gridExtra}.
#' @param ... Additional arguments passed to methods (ignored).
#'
#' @export
#' @method plot sim_traces
plot.sim_traces <- function(sim_data, level = c("individual", "group", "global"),
                            facet = TRUE, separate_pages = FALSE, width = 0.1, height = 0.2, size = 10,...) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for stacking subplots.")
  }

  plots_list <- list()
  level <- match.arg(level, several.ok = TRUE)

  ind_data <- sim_data$traces
  grp_data <- sim_data$group_traces

  if (!is.null(sim_data$global_trace)) {
    global_data <- sim_data$global_trace
  } else if (!is.null(grp_data)) {
    x_col <- if ("x" %in% names(grp_data)) "x" else "time"
    y_col <- if ("y_offset" %in% names(grp_data)) "y_offset" else "y"

    global_data <- grp_data |>
      dplyr::group_by(.data[[x_col]]) |>
      dplyr::summarise(y = mean(.data[[y_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::rename(x = 1)
  } else {
    global_data <- NULL
  }

  event_data <- NULL
  if ("sim_events" %in% class(sim_data) && !is.null(sim_data$events)) {
    event_data <- sim_data$events
  }

  n_groups <- sim_data$sim_parameters$n_groups
  purple_palette <- colorRampPalette(brewer.pal(9, "Purples"))(n_groups)

  build_base_plot <- function(ind_sub = NULL, grp_sub = NULL, ev_sub = NULL, title_suffix = "") {
    p <- ggplot() + theme_minimal() + labs(x = "Time (t)", y = "Eta(t)")

    if ("individual" %in% level && !is.null(ind_sub) && nrow(ind_sub) > 0) {
      i_x <- if ("x" %in% names(ind_sub)) "x" else "time"
      i_y <- if ("y_offset" %in% names(ind_sub)) "y_offset" else "y"

      p <- p + geom_line(data = ind_sub,
                         aes(x = .data[[i_x]], y = .data[[i_y]], group = interaction(group, ind),
                             color = as.factor(ind)), alpha = 0.6, linewidth = 0.5) +
        scale_color_viridis_d(guide = "none")
    }

    if ("group" %in% level && !is.null(grp_sub) && nrow(grp_sub) > 0) {
      g_x <- "x"
      g_y <- if ("y_offset" %in% names(grp_sub)) "y_offset" else "y"

      p <- p + geom_line(data = grp_sub,
                         aes(x = .data[[g_x]], y = .data[[g_y]], group = group),
                         color = "black",linewidth = 1) +
        labs(color = "Group")
    }

    if ("global" %in% level && !is.null(global_data) && nrow(global_data) > 0) {
      gl_x <- "x"
      gl_y <- if ("y_offset" %in% names(global_data)) "y_offset" else "y"

      p <- p + geom_line(data = global_data,
                         aes(x = .data[[gl_x]], y = .data[[gl_y]]),
                         color = "black", linewidth = 1.4, linetype = "dashed")
    }

    if (title_suffix != "") {
      p <- p + ggtitle(title_suffix)
    }

    if (!is.null(ev_sub) && nrow(ev_sub) > 0) {
      p_below <- ggplot(data = ev_sub,
                        aes(x = event_times,
                            y = factor(ind),
                            group = interaction(ind,group),
                            color = interaction(ind,group),
                            fill = interaction(ind,group))) +
        geom_tile(width = width, height = height) +
        scale_color_viridis_d(guide = "none") +
        scale_fill_viridis_d(guide = "none") +
        labs(x = "Time (t)", y = "Individual",) +
        theme(axis.text = element_text(size = size)) +
        theme_minimal()

      # Use patchwork to stack them cleanly, allocating more space to the upper trace
      return(p / p_below + plot_layout(heights = c(2, 1)))
    }

    return(p)
  }


  unique_groups <- unique(ind_data$group)

  if (!facet) {
    p <- build_base_plot(ind_data, grp_data, event_data, title_suffix = "Simulation Traces (All Layered)")
    return(p)

  } else {

    for (g in unique_groups) {
      i_sub <- ind_data |> dplyr::filter(group == g)
      g_sub <- grp_data |> dplyr::filter(group == g)

      ev_sub <- NULL
      if (!is.null(event_data)) {
        ev_sub <- event_data |> dplyr::filter(group == g)
      }

      p_g <- build_base_plot(i_sub, g_sub, ev_sub, title_suffix = paste("Group:", g))
      plots_list[[paste0("Group_", g)]] <- p_g
    }

    if (separate_pages) {
      if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for multi-page rendering.")
      }
      return(gridExtra::marrangeGrob(plots_list, nrow = 1, ncol = 1))

    } else {

      combined_plot <- patchwork::wrap_plots(plots_list, ) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(title = "Trace Simulation")

      return(combined_plot)
    }
  }
}
