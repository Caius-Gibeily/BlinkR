
#' @export
#' @method print gp_model
print.gp_model <- function(model, prior = FALSE, pars = NULL,
                           probs = c(0.25, 0.5, 0.75), digits_summary = 2, ...) {

  if (!prior & !is.null(model$fit)) {
    model_dat <- model$fit
  } else if (prior & !is.null(model$prior_pc) | "gp_prior_pc" %in% class(model)) {
    model_dat <- model$prior_pc
  }

  all_pars <- names(model_dat)

  pars_filt <- grep(c("z_|beta|eta|log_|diag|mu_raw_|lp_"), all_pars, value = TRUE, invert = TRUE)

  print(model_dat, pars = pars_filt,
        probs = probs, digits_summary = digits_summary, ...)

}
