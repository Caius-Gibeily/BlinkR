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

gp_fit <- function(traces, ...) {
  UseMethod("gp_fit")
}

#' @export
#' @method gp_fit tbl
gp_fit.tbl <- function(traces, group_id = "group", ind_id = "ind_id",
                       event_times = "event_times", dt = "dt",

                       M = 40, family = c("exponential","gamma","weibull","log-normal","gengamma"),
                       kernel = c("squared_exp", "matern12","matern32","matern52","periodic"),
                   priors = NULL, L_factor = 1.2) {

  ind_str <- tryCatch(rlang::as_name(rlang::ensym(ind_id)),
                      error = function(e) NULL)
  group_str <- tryCatch(rlang::as_name(rlang::ensym(group_id)),
                        error = function(e) NULL)

  event_times_str <- tryCatch(rlang::as_name(rlang::ensym(event_times)),
                     error = function(e) NULL)

  dt_str <- tryCatch(rlang::as_name(rlang::ensym(dt)),
                    error = function(e) NULL)

  if (is.null(ind_str) || !ind_str %in% names(traces)) {
    warning("Individual ID not provided or not found in 'traces'. Defaulting to 'one_ind' model")
    subclass <- "one_ind"
    ind_id <- rep(1, nrow(traces))
    g_id <- rep(1, nrow(traces))
    I <- 1
  } else {
    I <- traces %>%
      dplyr::distinct(.data[[ind_str]]) %>%
      dplyr::pull() %>%
      length()

    if (I == 1) {
      subclass <- "one_ind"
      ind_id <- rep(1, nrow(traces))
      g_id <- rep(1, nrow(traces))
    } else {
      ind_id <- traces %>% dplyr::select(.data[[ind_str]]) %>%
        dplyr::pull()
      if (is.null(group_str) || !group_str %in% names(traces)) {
        warning("Group ID not provided or not found in 'traces'. Defaulting model to 'one_group'.")
        subclass <- "one_group"
        g_id <- rep(1, nrow(traces))

      } else {

        G <- traces %>%
          dplyr::distinct(.data[[group_str]]) %>%
          dplyr::pull() %>%
          length()
        g_id <- traces %>% dplyr::select(.data[[group_str]]) %>%
          dplyr::pull()

        if (G == 1) {
          subclass <- "one_group"

        } else {
          subclass <- "multi_group"
        }
      }
    }
  }

  if (is.null(event_times_str) || !event_times_str %in% names(traces)) {
    stop("Please make sure that event times are present and that your column ID is the same.")
  }

  if (is.null(dt_str) || !dt_str %in% names(traces)) {
    traces <- traces %>%
      mutate(dt = c(0,diff(.data[[event_times_str]])))
  }
  g_membership <- traces %>%
    distinct(.data[[group_str]],.data[[ind_str]])

  model_setup <- list(
    event_times = traces %>%
      pull(.data[[event_times_str]]),

    dt = traces %>% pull(.data[[dt_str]]),
    ind_id = ind_id,
    g_id = g_id,
    g_membership = g_membership %>% pull(group),
    I = I,
    settings = list(M = M,
                    L_factor = L_factor,
                    kernel = match.arg(kernel),
                    family = match.arg(family),
                    priors = priors)
  )

  class(model_setup) <- c(subclass, "gp_fit_model")
  gp_fit(model_setup)

}

#' @export
#' @method gp_fit one_ind
gp_fit.one_ind <- function(model_setup, ...) {
  message("Fitting a single-individual Stan model")
  prior_frame <- parse_priors.one_ind(model_setup)


  model_setup$settings <- append(model_setup$settings,
                                 list(variables=prior_frame %>% pull(prior_variable)))
  flags <- .get_flags(model_setup)
  print(flags)

  stan_dat <- list(
    N_total = length(model_setup$event_times),
    M = model_setup$settings$M,
    L_factor = model_setup$settings$L_factor,
    t_ev = model_setup$event_times,
    dt = model_setup$dt,

    distributions = as.double(prior_frame$distribution_id),
    N_params = nrow(prior_frame),
    params = as.matrix(prior_frame[c("param_1","param_2")]),
    kernel = match(model_setup$settings$kernel,c("squared_exp",
                                        "matern12","matern32","matern52","periodic")),
    family = match(model_setup$settings$family, c("exponential","gamma","weibull",
                                                    "log-normal","gengamma"))
  )
  stan_dat <- append(stan_dat, flags)

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- stan(
    file = "inst/stan/hsgp_one_ind.stan",
    data = stan_dat,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)

  )
  return(fit)

}

#' @export
#' @method gp_fit one_group
gp_fit.one_group <- function(model_setup,flags, ...) {
  message("Fitting a single-group hierarchical Stan model")
  prior_frame <- parse_priors.one_group(model_setup)

  model_setup$settings <- append(model_setup$settings,
                                 list(variables=prior_frame %>%
                                        pull(prior_variable)))
  flags <- .get_flags(model_setup)
  print(flags)

  stan_dat <- list(
    N_total = length(model_setup$event_times),
    I = model_setup$I,
    ind_id = model_setup$ind_id,
    M = model_setup$settings$M,
    L_factor = model_setup$settings$L_factor,
    t_ev = model_setup$event_times,
    dt = model_setup$dt,

    distributions = as.double(prior_frame$distribution_id),
    N_params = nrow(prior_frame),
    params = as.matrix(prior_frame[c("param_1","param_2")]),
    kernel = match(model_setup$settings$kernel,c("squared_exp",
                                                 "matern12","matern32","matern52","periodic")),
    family = match(model_setup$settings$family, c("exponential","gamma","weibull",
                                                  "log-normal","gengamma"))
  )
  stan_dat <- append(stan_dat, flags)

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- stan(
    file = "inst/stan/hsgp_one_group.stan",
    data = stan_dat,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)

  )
  return(fit)
}

#' @export
#' @method gp_fit multi_group
gp_fit.multi_group <- function(model_setup,...) {
  message("Fitting a multi-group hierarchical Stan model")
  prior_frame <- parse_priors.multi_group(model_setup)

  model_setup$settings <- append(model_setup$settings,
                                 list(variables=prior_frame %>%
                                        pull(prior_variable)))
  flags <- .get_flags(model_setup)

  stan_dat <- list(
    N_total = length(model_setup$event_times),
    I = model_setup$I,
    ind_id = model_setup$ind_id,
    G = model_setup$G,
    g_id = model_setup$g_id,
    g_membership = model_setup$g_membership,
    I_per_group = tabulate(model_setup$g_membership),
    M = model_setup$settings$M,
    L_factor = model_setup$settings$L_factor,
    t_ev = model_setup$event_times,
    dt = model_setup$dt,

    distributions = as.double(prior_frame$distribution_id),
    N_params = nrow(prior_frame),
    params = as.matrix(prior_frame[c("param_1","param_2")]),
    kernel = match(model_setup$settings$kernel,c("squared_exp",
                                                 "matern12","matern32","matern52","periodic")),
    family = match(model_setup$settings$family, c("exponential","gamma","weibull",
                                                  "log-normal","gengamma"))
  )
  stan_dat <- append(stan_dat, flags)

  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- stan(
    file = "inst/stan/hsgp_multi_group.stan",
    data = stan_dat,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.95, max_treedepth = 12)

  )
  return(fit)
}


#' @noRd
.get_flags <- function(model_setup) {
  flags <- list()
  if (model_setup$settings$family != "lognormal") {
    flags[["include_mu_lognormal"]] = 0
    flags[["include_sigma_lognormal"]] = 0
  }

  if (model_setup$settings$family == "exponential") {
    flags[["include_k"]] = 0
    flags[["include_shape"]] = 0
  } else if (model_setup$settings$family == "gamma") {
    flags[["include_k"]] = match("k",model_setup$settings$variables)
    flags[["include_shape"]] = 0
  } else if (model_setup$settings$family == "weibull") {
    flags[["include_k"]] = 0
    flags[["include_shape"]] = match("shape",model_setup$settings$variables)
  } else if (model_setup$settings$family == "gengamma") {
    flags[["include_k"]] = match("k",model_setup$settings$variables)
    flags[["include_shape"]] = match("shape",model_setup$settings$variables)
  }

  return(flags)

}


