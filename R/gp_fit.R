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


gp_fit.default <- function(traces, group_id = group_id, ind_id = ind_id, event_times = event_times, dt = dt,
                   M = 40, family = c("exponential","gamma","weibull",
                                 "log-normal","gengamma"), priors = NULL, L_factor = 1.2) {

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

  model_setup <- list(
    event_times = traces %>%
      select(.data[[event_times_str]]),

    dt = traces %>% select(.data[[dt_str]]),
    ind_id = ind_id,
    g_id = g_id,
    I = I,
    settings = list(M = M,
                    L_factor = L_factor,
                    family = match.arg(family),
                    priors = priors)
  )

  class(model_setup) <- c(subclass, "gp_fit_model")

  gp_fit(model_setup,...)
}

gp_fit.one_ind <- function(model_setup, ...) {
  message("Fitting a single-individual Stan model")
  prior_frame <- parse_priors(model_setup)

}

gp_fit.one_group <- function(traces, ...) {
  message("Fitting a single-group hierarchical Stan model")
  prior_frame <- parse_priors(model_setup)
}

gp_fit.multi_group <- function(traces, ...) {
  message("Fitting a multi-group hierarchical Stan model...")
  prior_frame <- parse_priors(model_setup)
}

