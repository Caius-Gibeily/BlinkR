
#' @export
model_compare <- function(...) {
  models <- list(...)
  comp <-loo_compare(lapply(model_loo, models))
  return(comp)
}

#' @export
model_loo <- function(model) {
  if (!("gp_model" %in% class(model))) stop("Please provide a fitted GP model.")
  log_lik_array <- loo::extract_log_lik(model$fit,
                                   parameter_name = "log_lik",
                                   merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(log_lik_array))
  loo <- loo::loo(log_lik_array, r_eff = r_eff)
  return(loo)
}
