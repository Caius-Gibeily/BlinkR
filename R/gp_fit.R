#' Simulate inhomogeneous renewal processes
#'
#' @param sigma Numeric. Individual-level coefficient standard deviation.
#' @param Q Numeric. Mean coefficients for the spline basis.
#' @param baseline Numeric. Group-level coefficient standard deviation.
#' @param baseline_sd Numeric. Individual-level coefficient standard deviation.
#' @param lerp Numeric. Mean coefficients for the spline basis.
#'
#' @return A tibble containing simulated event times for each group.
#' @export

gp_fit <- function(events, ...) {
  UseMethod("gp_fit")
}

#' @export
#' @method gp_fit data.frame
gp_fit.data.frame <- function(events,...) {
  event_data <- list(events = events)
  gp_fit(event_data,...)
}

#' @export
#' @method gp_fit list
gp_fit.list <- function(event_data, duration, group_id = "group", ind_id = "ind_id",
                       event_times = "event_times", dt = "dt",

                       M = 40, family = c("exponential","gamma","weibull","log-normal","gengamma"),
                       kernel = c("squared_exp", "matern12","matern32","matern52","periodic"),
                       priors = NULL, L_factor = 1.2, w0 = 0.5,
                       chains = 4, iter = 2000, warmup = 1000, adapt_delta = 0.95, max_treedepth = 12,...) {
  kernel <- match.arg(kernel)
  family <- match.arg(family)

  if (kernel == "periodic") {
    warning("Periodic kernel uses half the number of basis functions (M). Doubling M")
    M = M * 2
  }

  events <- event_data$events

  ind_str <- tryCatch(rlang::as_name(rlang::ensym(ind_id)),
                      error = function(e) NULL)
  group_str <- tryCatch(rlang::as_name(rlang::ensym(group_id)),
                        error = function(e) NULL)

  event_times_str <- tryCatch(rlang::as_name(rlang::ensym(event_times)),
                     error = function(e) NULL)

  dt_str <- tryCatch(rlang::as_name(rlang::ensym(dt)),
                    error = function(e) NULL)

  if (missing(duration)) {
    duration <- events %>% reframe(max_t = max(.data[[event_times_str]])) %>% pull(max_t)
    warning("Setting duration to the latest event time. If this is incorrect, please set duration.")
  }
  if (is.null(ind_str) || !ind_str %in% names(events)) {
    warning("Individual ID not provided or not found in 'event_data'. Defaulting to 'one_ind' model")
    subclass <- "one_ind"
    I <- 1
    ind_id <- rep(1, nrow(events))
    g_id <- rep(1, nrow(events))
    g_membership <- 1
  } else {
    I <- events %>%
      dplyr::distinct(.data[[ind_str]]) %>%
      dplyr::pull() %>%
      length()

    if (I == 1) {
      subclass <- "one_ind"
      ind_id <- rep(1, nrow(events))
      g_id <- rep(1, nrow(events))
      g_membership <- 1
    } else {
      ind_id <- events %>% dplyr::select(.data[[ind_str]]) %>%
        dplyr::pull()
      if (is.null(group_str) || !group_str %in% names(events)) {
        warning("Group ID not provided or not found in 'events'. Defaulting to 'one_group' model.")
        subclass <- "one_group"
        g_id <- rep(1, nrow(events))
        g_membership <- rep(1, I)

      } else {

        G <- events %>%
          dplyr::distinct(.data[[group_str]]) %>%
          dplyr::pull() %>%
          length()
        g_id <- events %>% dplyr::select(.data[[group_str]]) %>%
          dplyr::pull()

        if (G == 1) {
          subclass <- "one_group"
          g_membership <- rep(1, I)
        } else {
          subclass <- "multi_group"
          g_membership <- events %>%
            distinct(.data[[group_str]],.data[[ind_str]]) %>%
            pull(group)

        }
      }
    }
  }

  if (is.null(event_times_str) || !event_times_str %in% names(events)) {
    stop("Please make sure that event times are present and that the column ID you specified has the same name.")
  }

  if (is.null(dt_str) || !dt_str %in% names(events)) {
    events <- events %>%
      mutate(dt = c(0,diff(.data[[event_times_str]])))
  }

  model_data <- list(
    event_times = events %>%
      pull(all_of(event_times_str)),
    traces = event_data$traces, # for simulations
    dt = events %>% pull(all_of(dt_str)),
    ind_id = as.numeric(as.character(ind_id)),
    g_id = as.numeric(as.character(g_id)),
    g_membership = g_membership ,
    I = I,
    settings = list(duration = duration,
                    M = M,
                    L_factor = L_factor,
                    w0 = w0,
                    kernel = match.arg(kernel),
                    family = match.arg(family),
                    priors = priors),
    stan_runtime = list(chains = chains,
                        iter = iter, warmup = warmup,
                        adapt_delta = adapt_delta,
                        max_treedepth = max_treedepth)
  )

  class(model_data) <- c(subclass, "gp_model")
  gp_fit(model_data)

}

#' @export
#' @method gp_fit gp_model
gp_fit.gp_model <- function(model_data, ...) {

  if ("one_ind" %in% class(model_data)) {
    message("Fitting a single-individual Stan model")
    prior_frame <- parse_priors.one_ind(model_data)
    file <- "inst/stan/hsgp_one_ind.stan"
  } else if ("one_group" %in% class(model_data)) {
    message("Fitting a single-group hierarchical Stan model")
    prior_frame <- parse_priors.one_group(model_data)
    file <- "inst/stan/hsgp_one_group.stan"
  } else if ("multi_group" %in% class(model_data)) {
    message("Fitting a multi-group hierarchical Stan model")
    prior_frame <- parse_priors.multi_group(model_data)
    file <- "inst/stan/hsgp_multi_group.stan"
  }

  model_data$settings <- append(model_data$settings,
                                 list(variables=prior_frame %>%
                                        pull(prior_variable)))
  flags <- .get_flags(model_data)

  stan_dat <- list(
    N_total = length(model_data$event_times),
    I = model_data$I,
    ind_id = model_data$ind_id,
    G = model_data$G,
    g_id = model_data$g_id,
    g_membership = model_data$g_membership,
    I_per_group = tabulate(model_data$g_membership),

    M = model_data$settings$M,
    L_factor = model_data$settings$L_factor,
    w0 = model_data$settings$w0,
    t_ev = model_data$event_times,
    dt = model_data$dt,

    distributions = as.double(prior_frame$distribution_id),
    N_params = nrow(prior_frame),
    params = as.matrix(prior_frame[c("param_1","param_2")]),
    kernel = match(model_data$settings$kernel,c("squared_exp",
                                                "matern12","matern32","matern52","periodic")),
    family = match(model_data$settings$family, c("exponential","gamma","weibull",
                                                 "log-normal","gengamma"))
  )
  stan_dat <- append(stan_dat, flags)

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  fit <- rstan::stan(
    file = file,
    data = stan_dat,
    chains = model_data$stan_runtime$chains,
    iter = model_data$stan_runtime$iter,
    warmup = model_data$stan_runtime$warmup,
    control = list(adapt_delta = model_data$stan_runtime$adapt_delta,
                   max_treedepth = model_data$stan_runtime$max_treedepth))

  model_data$fit <- fit
  return(model_data)
}

#' @noRd
.get_flags <- function(model_data) {
  flags <- list()
  if (model_data$settings$family != "lognormal") {
    flags[["include_mu_lognormal"]] = 0
    flags[["include_sigma_lognormal"]] = 0
  }

  if (model_data$settings$family == "exponential") {
    flags[["include_k"]] = 0
    flags[["include_shape"]] = 0
  } else if (model_data$settings$family == "gamma") {
    flags[["include_k"]] = match("k",model_data$settings$variables)
    flags[["include_shape"]] = 0
  } else if (model_data$settings$family == "weibull") {
    flags[["include_k"]] = 0
    flags[["include_shape"]] = match("shape",model_data$settings$variables)
  } else if (model_data$settings$family == "gengamma") {
    flags[["include_k"]] = match("k",model_data$settings$variables)
    flags[["include_shape"]] = match("shape",model_data$settings$variables)
  }

  return(flags)

}


