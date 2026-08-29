
#' @export
get_nlpd <- function(model, level = "ind", resolution = 0.1,
                     prior = FALSE, dev_only=FALSE, rescale = TRUE) {

  if (!rescale) warning("NLPD should be computed between traces with the same scale.")

  summary_traces <- summarise_traces(model,level = level, prior = prior,
                                 dev_only = dev_only, resolution = resolution, rescale = rescale)


  traces <- switch(level,
                   "ind" = model$traces,
                   "group" = model$group_traces,
                   "global" = model$global_trace)
  if (rescale) {
  scale_factor <- switch(model$sim_parameters$family,
                         "exponential" = 1,
                         "gamma" = model$sim_parameters$`k[1]`,
                         "weibull" = gamma(1 + (1/model$sim_parameters$`shape[1]`)),
                         "gengamma" =  gamma((model$sim_parameters$`shape[1]`+1)/model$sim_parameters$`k[1]`)/
                           gamma(model$sim_parameters$`shape[1]`/model$sim_parameters$`k[1]`))
  } else scale_factor <- 1

  if (level == "global") {
    summary_traces$y_ground <- approx(x = traces$x,
      y = traces$y_offset,
      xout = summary_traces$x)$y
  } else {
    summary_traces <- summary_traces |>
      dplyr::group_by(.data[[level]]) |>
      dplyr::group_modify(~ {
        current_id <- .y[[level]]
        traces_sub <- dplyr::filter(traces, .data[[level]] == current_id)

        .x$y_ground <- approx(x = traces_sub$x,
          y = traces_sub$y_offset,
          xout = .x$x)$y
        .x}) |>
      dplyr::ungroup()
  }

  summary_traces <- summary_traces |>
    tidyr::drop_na(y_ground, means, sds) |>
    dplyr::mutate(dnorms = dnorm(y_ground + log(scale_factor),
                                 means, sds, log = TRUE)) |>
    dplyr::group_by(.data[[level]]) |>
    dplyr::summarise(mean_NLPD = -mean(dnorms),
                     sum_NLPD = -sum(dnorms),
                     .groups="drop")
  return(summary_traces)
}






