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
simulate_events <- function(traces, group = group, ind = row_id, family = c("exponential","gamma","weibull",
                                                 "log-normal","gengamma"),
                            scale = NULL, mu = NULL, shape = 1, k = 2, sigma = 1,
                            Q = 0, baseline = 0, baseline_sd = 0, lerp = 10) {
  sim_struct <- traces |>
    group_by({{group}},{{ind}}) |> group_keys()

  if (is.null(scale)) {
    if (length(baseline) == 1) {
      traces <- traces |>
        mutate(scale = 1/exp(-(y + baseline)))
    } else {
      traces <- traces |> group_by(row_id) |>
        mutate(scale = exp(y + rnorm(n(),0,baseline_sd)))
    }

  }

  if (family != "log-normal") {
    if (family == "exponential") shape <- k <- 1
    else if (family == "gamma") shape <- 1
    else if (family == "weibull") k <- 1

    event_times <- traces |>
      group_split(row_id) |>
      map(\(.x) .simulate_renewal(.x, modulant = scale,
                                  shape = shape, k = k, lerp = lerp)) |>
      list_rbind(names_to = "row_id") |>
      left_join(sim_struct, by = "row_id")|> relocate(group,.after=row_id)


  } else if (family == "log-normal") {
    if (is.null(mu)) {
    traces <- traces |>
      mutate(mu = log(scale) + log(k)/sqrt(shape))
    }

    event_times <- traces |>
      group_split(row_id) |>
      map(\(.x) .simulate_renewal(.x, modulant = mu,
                                  sigma = sigma, Q = 0, lerp = lerp)) |>
      list_rbind(names_to = "row_id") |>
      left_join(sim_struct, by = "row_id") |> relocate(group,.after=row_id)
  } else stop("Please choose a survival function from ",
              "exponential, gamma, weibull, log-normal, gengamma")
  return(event_times)

}

#' @noRd
.simulate_renewal <- function(trace,modulant,shape,k,sigma,Q,lerp = 10) {


  trace_parts <- trace |>
    dplyr::select({{modulant}}) |> dplyr::pull()
  trace_parts <- trace_parts |> approx(n = length(trace_parts) * lerp)

  modulant <- trace_parts$y
  time <- trace_parts$x

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

