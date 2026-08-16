
#' @export
mcmc_plot <- function(model,...) {
  UseMethod("mcmc_plot")
}

#' @export
#' @method mcmc_plot gp_model
mcmc_plot.gp_model <- function(model,pars=NULL,type = "intervals",scheme="blue",...) {

  all_pars <- names(model$fit)
  pars_filt <- grep(c("z_|beta|eta|log_|diag|mu_raw_|lp_"), all_pars, value = TRUE, invert = TRUE)

  pars_types <- list(mu = c("mu"), alpha = c("alpha"), rho = c("rho"), survival = c("k|shape"))
  if (is.null(pars)) {
    sub_pars <- lapply(pars_types, function(type) {
      grep(type, pars_filt, value = TRUE, invert = FALSE)
    })
  } else {
    matching_pars <- lapply(pars, function(type) {
      grep(type, pars_filt, value = TRUE, invert = FALSE)
    })
    sub_pars <- lapply(pars_types, function(type) {
      grep(type, purrr::flatten(matching_pars), value = TRUE, invert = FALSE)
    })
  }

  plot_fn <- switch(type,
                    "areas" = bayesplot::mcmc_areas,
                    "dens" = bayesplot::mcmc_dens,
                    "dens_chains" = bayesplot::mcmc_dens_chains,
                    "hist" = bayesplot::mcmc_hist,
                    "intervals" = bayesplot::mcmc_intervals,
                    "trace" = bayesplot::mcmc_trace,
                    stop("Unsupported plot type")
  )

  bayesplot::color_scheme_set(scheme = scheme)

  plots <- lapply(sub_pars, function(p) {
    if (!rlang::is_empty(p)) {
      plot_fn(model$fit, pars = p, ...) + ggplot2::theme_minimal()
    }
  }) |> purrr::discard(is.null)

  .ggplot_successive(plots)
}
