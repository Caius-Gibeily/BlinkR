

#' Title
#'
#' @param model
#' @param pars
#' @param prior
#' @param type
#' @param scheme
#' @param show_ground
#' @param height
#' @param linewidth
#' @param color
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
mcmc_plot <- function(model,pars=NULL,prior=FALSE,
                               type = c("areas", "dens", "dens_chains", "hist", "intervals", "trace"),
                      scheme="blue",show_ground = TRUE,height=0.3, linewidth=1,color="black",...) {
  if (length(type) > 1) type = "intervals" # default option

  type <- match.arg(type)

  if (!prior & !is.null(model$fit)) {
    model_dat <- model$fit
  } else if (prior & !is.null(model$prior_pc) | "gp_prior_pc" %in% class(model)) {
    model_dat <- model$prior_pc
  }

  all_pars <- names(model_dat)

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
      p_layer <- plot_fn(model_dat, pars = p, ...) +
        ggplot2::theme_minimal()
      if (show_ground & !is.null(model$sim_parameters) & any(p %in% names(model$sim_parameters))) {

        params <- model$sim_parameters |>
          purrr::keep_at(p) |>
          dplyr::as_tibble() |>
          tidyr::pivot_longer(cols= everything(),
                       names_to="parameter",
                       values_to = "xintercept")
        y_levels <- ggplot2::ggplot_build(p_layer)$layout$panel_params[[1]]$y$get_labels()

        params <- params |>
          dplyr::filter(parameter %in% y_levels) |>
          dplyr::mutate(parameter = factor(parameter, levels = y_levels))

        p_layer +
          ggplot2::geom_linerange(
            data = params,
            aes(x = xintercept, y = parameter, ymin = as.numeric(parameter) - height, ymax = as.numeric(parameter) + height),
            inherit.aes = FALSE,
            color = color, linewidth = linewidth
          )

      } else p_layer
    }

  }) |> purrr::discard(is.null)

  .ggplot_successive(plots)
}
