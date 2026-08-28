#' Simulate inhomogeneous renewal processes
#'
#' @param traces Tibble or Data.Frame object. A dataframe with columns containing group and individual id and their simulated latent trajectory.
#' @param group Integer. Number of individuals.
#' @param ind Integer. Number of groups.
#' @param family Character vector. Choice of survival family from which to draw renewal events. Exponential, gamma, log-normal and generalised gamma are supported.
#' @param scale Equal to the inverse of rate. The scale is the time-varying
#' @param mu Location parameter of the log-normal and generalised gamma families. It is an alternative parameterisation of the generalised gamma family relative to the Stacey (1962) shape, k and scale parameters.
#' @param shape Weibull shape parameter. When shape = 1 and k != 1, the distribution becomes a gamma distribution.
#' @param k Numeric. Group-level coefficient standard deviation.
#' @param sigma Numeric. Individual-level coefficient standard deviation.
#' @param Q Numeric. Mean coefficients for the spline basis.
#' @param baseline Numeric. Group-level coefficient standard deviation.
#' @param baseline_sd Numeric. Individual-level coefficient standard deviation.
#' @param lerp Numeric. Mean coefficients for the spline basis.
#'
#' @return A tibble containing simulated event times for each group.
#' @export
simulate_events <- function(trace_data, group = group, ind = ind, family = c("exponential","gamma","weibull",
                                                 "log-normal","gengamma"),shape = 1, k = 2, sigma = 1,
                            Q = 0, baseline = 0, lerp = 100) {

  # extract trace data from list if the object is already a list
  if (is.list(trace_data) & !is.data.frame(trace_data)) traces <- trace_data$traces

  # else, extract trace from dataframe and convert to a list container
  else if (is.data.frame(trace_data)) {
    traces <- trace_data
    trace_data <- list(traces = traces)
  }
  sim_struct <- trace_data$traces <- traces |>
    group_by({{group}},{{ind}}) |> group_keys()
  n_ind <- nrow(sim_struct)
  n_group <- unique(sim_struct$group) |> length()
  if (length(baseline) == 1) {
    # if a global mu offset is specified, translate all traces by this fixed value

    traces <- trace_data$traces <- traces |>
      mutate(y_offset = y + baseline,
             scale = 1/exp(-(y_offset)))
    if (!is.null(trace_data$group_traces)) {
      # if the group trace is also present (multi-ind or multi-group models), add baseline to group traces
      trace_data$group_traces <- trace_data$group_traces |>
        mutate(y_offset = y + baseline)
    }
    if (!is.null(trace_data$global_traces)) {
      # if the global trace is also present (multi-group model), add baseline to group traces
      trace_data$global_trace <- trace_data$global_trace |>
        mutate(y_offset = y + baseline)

      trace_data$sim_parameters$mu_global <- baseline
    }

    mu_inds <- as.list(setNames(rep(baseline,n_ind), paste0("mu_ind[",sim_struct$ind, "]")))
    mu_groups <- as.list(setNames(rep(baseline, n_group),
                                  paste0("mu_group[",
                                         unique(sim_struct$group), "]")))

  } else if (length(baseline) == n_ind) {
    # if the number of offsets supplied equals the number of individuals, add each
    # separately to each individual
    traces <- trace_data$traces <- traces |> group_by(ind) |>
      mutate(y_offset = y + baseline[ind],
             scale = exp(y_offset)) |>
      dplyr::ungroup()

    # group baselines
    baseline_groups <- sim_struct |>
      dplyr::bind_cols(baseline=baseline) |>
      group_by(group) |>
      summarise(baseline = mean(baseline)) |>
      dplyr::pull(baseline)

    trace_data$group_traces <- trace_data$group_traces |> group_by(group) |>
      mutate(y_offset = y + baseline_groups) |>
      dplyr::ungroup()

    mu_inds <- as.list(setNames(baseline, paste0("mu_ind[", sim_struct$ind, "]")))
    mu_groups <- as.list(setNames(baseline_groups,
                                  paste0("mu_group[",
                                         unique(sim_struct$group), "]")))

    if (!is.null(trace_data$global_trace)) {
      baseline_global <- sim_struct |>
        dplyr::bind_cols(baseline=baseline) |>
        group_by(group) |>
        summarise(baseline = mean(baseline)) |>
        dplyr::pull(baseline)

      trace_data$global_trace <- trace_data$global_trace |>
        mutate(y_offset = y + baseline_global)

      trace_data$sim_parameters$mu_global <- baseline_global
    }

  } else stop("Please ensure a single baseline offset or a number of baseline offsets equal to the number of individuals")

  trace_data$sim_parameters <- trace_data$sim_parameters |>
    append(c(mu_inds,mu_groups))


  if (family != "log-normal") {
    if (family == "exponential") shape <- k <- 1
    else if (family == "gamma") shape <- 1
    else if (family == "weibull") k <- 1

    event_times <- traces |>
      group_split(ind) |>
      map(\(x) .simulate_renewal(x, modulant = scale,
                                  shape = shape, k = k, lerp = lerp)) |>
      list_rbind(names_to = "ind") |> mutate(ind = as.factor(ind)) |>
      left_join(sim_struct, by = "ind")|> relocate(group,.after=ind)

    trace_data$sim_parameters <- trace_data$sim_parameters |>
      append(list(`shape[1]` = shape, `k[1]`= k))

  } else if (family == "log-normal") {
    if (is.null(mu)) {
    traces <- traces |>
      mutate(mu = log(scale) + log(k)/sqrt(shape))
    }

    event_times <- traces |>
      group_split(ind) |>
      map(\(x) .simulate_renewal(x, modulant = mu,
                                  sigma = sigma, Q = 0, lerp = lerp)) |>
      list_rbind(names_to = "ind") |> mutate(ind = as.factor(ind)) |>
      left_join(sim_struct, by = "ind") |> relocate(group,.after=ind)
  } else stop("Please choose a survival function from ",
              "exponential, gamma, weibull, log-normal, gengamma")
  if (is.list(trace_data)) {
    sim_data <- trace_data %>% append(
      list(events = event_times))
    class(sim_data) <- c("sim_traces",class(sim_data))
  } else {
    sim_data <- list(traces = trace_data,
                     events = event_times)
  }
  class(sim_data) <- c("sim_events",class(sim_data))
  return(sim_data)

}

#' @noRd
.simulate_renewal <- function(trace,modulant,shape,k,sigma,Q,lerp = 20) {

  if (!is.null(lerp)) {
    trace_parts <- trace |> pull({{modulant}})
    trace_parts <- trace_parts |> approx(n = max(trace$x) * lerp)

    modulant <- trace_parts$y
    time <- trace_parts$x

  } else {
    modulant <- trace |> pull({{modulant}})
    time <- trace$x
  }

  if (!missing(shape) && !missing(k)) {
    events <- simulate_renewal_orig(time,modulant,shape,k)
  } else if (!missing(sigma) && !missing(Q)) {
    events <- simulate_renewal(time,modulant,sigma,Q)
  }
  renewal_events <- tibble(event_times = events,
                        dt = diff(c(0,events)))
  return(renewal_events)
}
#
#   event_times <- tibble(t_event = 1,
#                         dt = mean({{modulant}}),
#                         censored = TRUE)
#   t_diff <- time[2] - time[1]
#   idx <- 1
#
#   event_times <- numeric(1000)
#   event_times[1] <- last_blink <- 1
#
#   for (ti in time[-1]) {
#
#     dt = ti - last_blink
#
#     if (!missing(shape) && !missing(k)) {
#       h_t <-  t_diff * flexsurv::dgengamma.orig(dt,shape = shape, scale = modulant[which.min(abs(time - ti))],
#                                        k = k) /
#         (1 - flexsurv::pgengamma.orig(dt,shape = shape, scale = modulant[which.min(abs(time - ti))], k = k))
#       } else if (!missing(sigma) && !missing(Q)) {
#       h_t <-  t_diff * flexsurv::dgengamma(dt,mu = modulant[which.min(abs(time - ti))],
#                                            sigma = sigma, Q = Q) /
#         (1 - flexsurv::pgengamma(dt,mu = modulant[which.min(abs(time - ti))],
#                                       sigma = sigma, Q = Q))
#
#     }
#
#       if (runif(1) <= h_t) {
#
#         idx <- idx + 1
#         if (idx > length(event_times)) length(event_times) <- length(event_times) * 1.5
#
#         event_times[idx] <- ti
#         last_blink <- ti
#       }
#
#   }
#   event_times = event_times[event_times != 0]
#   renewal_events <- tibble(event_times = event_times,
#                            dt = diff(c(0,event_times)))
# }
#

