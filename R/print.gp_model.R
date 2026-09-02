#' Print table of inferred parameters
#' @inheritParams plot.gp_model
#' @param probs Numeric or numeric vector. Parameter quantiles to print. By default
#' `probs=c(0.25,0.5,0.75)`, printing 25th, 50th and 75th percentile parameter values
#' @param digits_summary Integer. Number of significant digits to round to. By default,
#' this is 2.
#' @param ... Other optional arguments
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

  pars_filt <- grep(c("z_|beta|eta|log_|diag|mu_raw_|lp_"),
    all_pars,
    value = TRUE, invert = TRUE
  )

  print(model_dat,
    pars = pars_filt,
    probs = probs, digits_summary = digits_summary, ...
  )
}

#' @export
#' @method print sim_traces
print.sim_traces <- function(trace_dat) {

  params <- trace_dat$sim_parameters
  param_names <- names(params)
  if (is.null(param_names) || any(param_names == "")) {
    param_names <- paste0("[", seq_along(params), "]")
  }


  cat("# Simulation parameters:\n")

  max_width <- max(nchar(param_names))

  for (i in seq_along(params)) {
    padded_name <- sprintf(paste0("%-", max_width, "s"),
                           param_names[i])

    val <- params[[i]]
    val_str <- if (length(val) > 1) {
      paste0("[", paste(head(val, 5), collapse = ", "),
             if(length(val) > 5) "..." else "", "]")
    } else {
      format(val)
    }

    cat(paste0("  ", padded_name, " : ", val_str, "\n"))
  }

  invisible(params)
}





