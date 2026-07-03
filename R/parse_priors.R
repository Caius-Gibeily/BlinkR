#' @param sigma Numeric. Individual-level coefficient standard deviation.
#' @param Q Numeric. Mean coefficients for the spline basis.
#' @param baseline Numeric. Group-level coefficient standard deviation.
#' @param baseline_sd Numeric. Individual-level coefficient standard deviation.
#' @param lerp Numeric. Mean coefficients for the spline basis.
#'
#' @return A tibble containing simulated event times for each group.
#' @export

parse_priors <- function(model_setup, ...) {
  UseMethod("parse_priors")
}

parse_priors.one_ind <- function(model_setup) {
  valid_variables = c("z","mu","alpha","rho") %>%
    append(.get_family_priors(model_setup$settings$family))

  prior_frame <- .process_prior_frame(path = "inst/extdata/default_priors_one-ind.csv",
                                      valid_variables = valid_variables,
                                      priors = model_setup$settings$priors)
}

parse_priors.one_group <- function(model_setup) {
  valid_variables = c("z_group","z_ind", "mu","alpha_group",
                      "alpha_ind", "rho_group","rho_ind") %>%
    append(.get_family_priors(model_setup$settings$family))

  prior_frame <- .process_prior_frame(path = "inst/extdata/default_priors_one-group.csv",
                                      valid_variables = valid_variables,
                                      priors = model_setup$settings$priors)
}

parse_priors.multi_group <- function(model_setup) {
  valid_variables = c("z_pop","z_group","z_ind","mu_group","mu_ind","alpha_pop",
                      "alpha_group","alpha_ind", "rho_pop","rho_group","rho_ind") %>%
    append(.get_family_priors(model_setup$settings$family))

  prior_frame <- .process_prior_frame(path = "inst/extdata/default_priors_multi-group.csv",
                                      valid_variables = valid_variables,
                                      priors = model_setup$settings$priors)
}

.get_family_priors <- function(family) {

  if (family == "gamma") return("k")
  else if (family == "weibull") return("shape")
  else if (family == "gengamma") return(c("k","shape"))
  else if (family == "lognormal") return(c("mu_lnorm","sigma_lnorm"))
}

.process_prior_frame <- function(path, valid_variables, priors = NULL) {
  # e.g.
  # prior_set <- alist(
  #   z_pop ~ normal(0,1),
  #   z_group ~ normal(0,1),
  #   z_ind ~ normal(0,1),
  #
  #   mu_group ~ normal(0,1),
  #   mu_ind ~ normal(0,1),
  #
  #   alpha_pop ~ normal(0,1),
  #   alpha_group ~ normal(0,1),
  #   alpha_ind ~ normal(0,1),
  #
  #   rho_pop ~ lognormal(0.3,0.2),
  #   rho_group ~ lognormal(0.3,0.2),
  #   rho_ind ~ lognormal(0.3,0.2),
  #
  #   k ~ exp(0.5),
  #   shape ~ exp(0.5)
  # )

  prior_frame <- readr::read_csv(csv_path, show_col_types = FALSE)

  if (missing(priors) || is.null(priors)) {
    return(prior_frame)
  }

  for (f in priors) {
    formula <- as.formula(f)
    prior_variable <- as.character(formula[[2]])

    if (!(prior_variable %in% valid_variables)) {
      stop(paste0("Variable '", prior_variable, "' is not valid for this model type. ",
                  "Allowed parameters: ", paste(valid_variables, collapse = ", ")))
    }

    dist_data <- formula[[3]]
    dist_name <- as.character(dist_data[[1]])

    if (!(dist_name %in% c("normal", "lognormal", "cauchy", "exp"))) {
      stop("Supported distributions: normal, lognormal, cauchy, or exp.")
    }

    dist_args <- as.list(dist_data)[-1]
    if (any(!sapply(dist_args, is.numeric))) {
      stop("Prior parameters must be numeric constants.")
    }

    prior_frame[prior_frame$prior_variable == prior_variable, "distribution"] <- dist_name

    if (length(dist_args) == 1) {
      prior_frame[prior_frame$prior_variable == prior_variable, "param_1"] <- dist_args[[1]]
      prior_frame[prior_frame$prior_variable == prior_variable, "param_2"] <- NA
    } else {
      prior_frame[prior_frame$prior_variable == prior_variable, c("param_1", "param_2")] <- dist_args
    }
  }

  return(prior_frame)
}
