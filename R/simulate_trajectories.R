#' Simulate Trajectories
#'
#' @param duration Integer. The length of the time series.
#' @param n_ind Integer. Number of individuals.
#' @param n_groups Integer. Number of groups.
#' @param class Character. Type of simulation: "spline", "boxcar", or "gp".
#' @param df Integer. Degrees of freedom for splines.
#' @param sd_global Numeric. Global standard deviation.
#' @param sd_group Numeric. Group-level standard deviation.
#' @param sd_ind Numeric. Individual-level standard deviation.
#' @param mu_coeffs Numeric vector. Mean coefficients for splines.
#' @param nodes Numeric vector. Node locations for boxcar simulation.
#' @param steps Numeric vector. Step sizes for boxcar simulation.
#' @param scale_factors_group Numeric vector. Scaling factors per group.
#' @param cycles Integer. Number of cycles for boxcar.
#' @param kernel Character. GP kernel choice.
#' @param alpha_global,rho_global,alpha_group,rho_group,alpha_ind,rho_ind,nu,P GP hyper-parameters.
#'
#' @importFrom magrittr |>
#' @importFrom stats rnorm expand.grid findInterval
#' @export
simulate_trajectories <- function(sim, duration = 100, n_ind = 10, n_groups = 1, ...) {
  UseMethod("simulate_trajectories", sim)
}

#' @export
simulate_trajectories.sim_spline <- function(sim, duration = 100, n_ind = 10, n_groups = 1, ...) {
  sim_struct <- tibble::rowid_to_column(tidyr::expand_grid(group = 1:n_groups, ind = 1:n_ind), var = "row_id")

  mu_coeffs <- sim$mu_coeffs
  if (n_ind == 1 || n_groups == 1) mu_coeffs <- 0

  group_traces <- .simulate_n(duration, n = n_groups, class = "spline",
                              df = sim$df, mu_coeffs = mu_coeffs, sd_coeffs = sim$sd_group)

  ind_traces <- .simulate_n(duration, n = n_ind * n_groups, class = "spline",
                            df = sim$df, mu_coeffs = group_traces$coeffs[, sim_struct$group, drop = FALSE],
                            sd_coeffs = sim$sd_ind)

  ind_traces$traces <- ind_traces$traces %>%
    dplyr::left_join(sim_struct, by = "row_id") %>%
    dplyr::relocate(group, ind)

  return(list(group_traces = group_traces, ind_traces = ind_traces))
}

#' @export
simulate_trajectories.sim_boxcar <- function(sim, duration = 100, n_ind = 10, n_groups = 1, ...) {
  sim_struct <- tibble::rowid_to_column(stats::expand.grid(group = 1:n_groups, ind = 1:n_ind), var = "row_id")

  scale_factors_group <- stats::rnorm(n_groups, 1, sim$sd_group %||% 0.5)
  if (n_groups == 1 || n_ind == 1) scale_factors_group <- 1

  scale_factors_ind <- scale_factors_group[sim_struct$group] + stats::rnorm(n_groups * n_ind, 0, sim$sd_ind)

  group_traces <- .simulate_n(duration, n = n_groups, class = "boxcar",
                              nodes = sim$nodes, steps = sim$steps, cycles = sim$cycles,
                              scale_factor = scale_factors_group)

  ind_traces <- .simulate_n(duration, n = n_ind * n_groups, class = "boxcar",
                            nodes = sim$nodes, steps = sim$steps,
                            cycles = sim$cycles, scale_factor = scale_factors_ind)

  ind_traces$traces <- ind_traces$traces %>%
    dplyr::left_join(sim_struct, by = "row_id") %>%
    dplyr::relocate(group, ind)

  return(list(group_traces = group_traces, ind_traces = ind_traces))
}

#' @export
simulate_trajectories.sim_gp <- function(sim, duration = 100, n_ind = 10, n_groups = 1, ...) {
  sim_struct <- tibble::rowid_to_column(stats::expand.grid(group = 1:n_groups, ind = 1:n_ind), var = "row_id")

  alpha_global <- sim$alpha_global
  alpha_group <- sim$alpha_group

  if (n_ind == 1) {
    alpha_global <- alpha_group <- 0
  } else if (n_groups == 1) {
    alpha_global <- 0
  }

  global_trace <- .simulate_n(duration, n = 1, class = "gp",
                              kernel = sim$kernel, alpha = alpha_global,
                              rho = sim$rho_global, nu = sim$nu, P = sim$P, level = "global")

  group_traces <- .simulate_n(duration, n = n_groups, class = "gp",
                              kernel = sim$kernel, alpha = alpha_group,
                              rho = sim$rho_group, nu = sim$nu, P = sim$P, level = "group") %>%
    dplyr::mutate(y = y + global_trace$y)

  ind_traces <- .simulate_n(duration, n = n_ind, class = "gp",
                            kernel = sim$kernel, alpha = sim$alpha_ind,
                            rho = sim$rho_ind, nu = sim$nu, P = sim$P) %>%
    dplyr::left_join(sim_struct, by = "row_id") %>%
    dplyr::left_join(group_traces, by = c("group", "x"), suffix = c("_ind", "_group")) %>%
    dplyr::mutate(y = y_ind + y_group) %>%
    dplyr::relocate(group, ind, .after = "row_id")

  return(list(group_traces = group_traces, ind_traces = ind_traces))
}


.simulate_n <- function(duration, n, class, df = 5, mu_coeffs = 0, sd_coeffs = 0.5,
                        level = "row_id", nodes, steps, cycles, scale_factor = 1, ...) {

  if (class == "spline") {
    if (length(sd_coeffs) != 1 && length(sd_coeffs) != df) {
      stop("Please enter a vector of spline coefficient standard deviations equal to df or a single value.")
    }
    if (length(mu_coeffs) != 1 && .get_first_dim(mu_coeffs) != df) {
      stop("Please enter a vector of mean spline coefficients equal to df or a single value.")
    }

    coeffs <- matrix(stats::rnorm(n * df, mu_coeffs, sd_coeffs), nrow = df, ncol = n)

    traces <- purrr::map(1:n, \(i) .sim_splines(duration, df, coeffs[, i], ...)) |>
      purrr::list_rbind(names_to = level)

    return(list(traces = traces, coeffs = coeffs))

  } else if (class == "boxcar") {
    level_steps <- as.matrix(steps) %*% t(scale_factor)

    traces <- purrr::map(1:n, \(i) .sim_boxcar(duration, nodes, level_steps[, i], cycles)) |>
      purrr::list_rbind(names_to = level)

    return(list(traces = traces, scales = level_steps))

  } else if (class == "gp") {
    traces <- purrr::map(1:n, \(i) .sim_GP(duration, ...)) |>
      purrr::list_rbind(names_to = level)

    return(traces)
  }
}

.sim_splines <- function(duration, df, coeffs) {
  basis_matrix <- splines::bs(seq(1, duration), df = df)
  y <- as.vector(basis_matrix %*% coeffs)
  return(tibble::tibble(x = seq(1, duration), y = y))
}

.sim_GP <- function(duration, kernel = "squared_exp", alpha = 0.5,
                    rho = duration / 5, nu = 3/2, P = duration / 2) {
  cov_mat <- .kernel(duration, kernel, alpha, rho, nu, P)
  y <- MASS::mvrnorm(n = 1, mu = rep(0, duration), Sigma = cov_mat)
  return(tibble::tibble(x = seq(1, duration), y = y))
}

.sim_boxcar <- function(duration, nodes, steps, cycles) {
  if (nodes[length(nodes)] < duration) {
    cycles <- 1
  }
  if (any(diff(nodes) / cycles < 1)) {
    warning("Number of cycles smaller than step grid resolution")
  }

  interval_id <- stats::findInterval(1:(nodes[length(nodes)] / cycles),
                                     nodes / cycles, rightmost.closed = TRUE) |>
    rep_len(duration)

  return(tibble::tibble(x = 1:duration, y = steps[interval_id]))
}

.kernel <- function(duration, kernel, alpha, rho, nu, P) {
  X <- seq(1, duration, by = 1)

  # Replace 'KernelPackage' with the actual package name providing these functions
  if (kernel == "squared_exp") {
    K <- rkriging::Gaussian.Kernel(rho)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "matern") {
    K <- rkriging::Matern.Kernel(rho, nu = nu)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "rational_quad") {
    K <- rkriging::RQ.Kernel(rho, alpha = 1)
    cov_mat <- alpha^2 * rkriging::Evaluate.Kernel(K, X)
  } else if (kernel == "periodic") {
    cov_mat <- alpha^2 * .periodic_kernel(duration, rho, P)
  } else {
    stop("Please choose among squared_exp, matern, rational_quad and periodic")
  }
  return(cov_mat)
}

.periodic_kernel <- function(duration, rho, P) {
  grid <- seq(1, duration, by = 1)
  dist_mat <- abs(outer(grid, grid, "-"))
  cov_mat <- exp(-(2 * sin(pi * dist_mat / P)^2) / rho^2)
  return(cov_mat)
}

.get_first_dim <- function(x) {
  if (is.null(dim(x))) length(x) else dim(x)[1]
}
