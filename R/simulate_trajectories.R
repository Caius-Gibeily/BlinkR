
#' @param name description
#' @example path.R
#' @export

simulate_trajectories <- function(duration = 100, n_ind = 10,
                                  n_groups = 1, class = "spline", df = 5, sd_global = 1,
                                  sd_group = 0.5, sd_ind = 0.25, mu_coeffs = rnorm(df,0,sd_global),
                                  formula = NULL,
                                  # if boxcar selected
                                  nodes = seq(1,duration, length.out = 3), steps = c(0,1), cycles = 2,
                                  # if kernels selected
                                  kernel = "squared_exp", alpha_global = 0.5, rho_global = duration / 2,
                                  alpha_group = 0.5, rho_group = duration / 2,
                                  alpha_ind = 0.25, rho_ind = duration / 2, nu = 3/2, P = 1) {

  sim_struct <- expand.grid(group = 1:n_groups,
                            ind = 1:n_ind) %>% rowid_to_column(var = "row_id")
  if (class == "spline") {

    if (n_ind == 1 | n_groups == 1) mu_coeffs = 0

    group_traces <- .simulate_n(duration, n = n_groups, class = class,
                                df = df, mu_coeffs = mu_coeffs, sd_coeff = sd_group)
    ind_traces <- .simulate_n(duration, n = n_ind * n_groups, class = class,
                              df = df, mu_coeffs = group_traces$coeffs[,sim_struct$group],
                              sd_coeff = sd_ind)
    ind_traces$traces %>% left_join(sim_struct, by = "row_id") %>%
      relocate(group, ind) -> ind_traces$traces
    return(list(group_traces = group_traces,
                ind_traces = ind_traces))

  } else if (class == "boxcar") {

    scale_factors_group <- rnorm(n_groups, 1, sd_group)
    scale_factors_ind <- scale_factors_group[sim_struct$group] + rnorm(n_groups * n_ind, 1, sd_ind)

    group_traces <- .simulate_n(duration, n = n_groups, class = class,
                                nodes = nodes, steps = steps, cycles = cycles,
                                scale_factor = scale_factors_group)
    ind_traces <- .simulate_n(duration, n = n_ind * n_groups, class = class,
                              nodes = nodes, steps = steps,
                              cycles = cycles,scale_factor = scale_factors_ind)
    ind_traces$traces <- ind_traces$traces %>% left_join(sim_struct, by = "row_id") %>%
      relocate(group, ind)

    return(list(group_traces = group_traces,
                ind_traces = ind_traces))

  } else if (class == "gp") {

    if (n_ind == 1) alpha_global = alpha_group = 0
    else if (n_groups == 1) alpha_global = 0

    global_trace <- .simulate_n(duration, n = 1, class = class,
                                kernel = kernel, alpha = alpha_global,
                                rho = rho_global,nu = nu, P = P, level = "global")
    group_traces <- .simulate_n(duration, n = n_groups, class = class,
                                kernel = kernel, alpha = alpha_group,
                                rho = rho_group,nu = nu, P = P, level = "group") %>%
      mutate(y = y + global_trace$y)

    ind_traces <- .simulate_n(duration, n = n_ind, class = class,
                              kernel = kernel, alpha = alpha_ind,
                              rho = rho_ind,nu = nu, P = P) %>%
      left_join(sim_struct, by = "row_id") %>%
      left_join(group_traces,
                by = c("group","x"),
                suffix = c("_ind","_group")) %>%
      mutate(y = y_ind + y_group) %>%
      relocate(group, ind,.after = "row_id") -> ind_traces

    return(list(group_traces = group_traces,
                ind_traces = ind_traces))

  } else if (!(class %in% c("spline","gp","boxcar"))) {

      stop("Please choose a class from spline, gp or boxcar")
    }

  }



.simulate_n <- function(duration, n, class, df = 5, mu_coeffs = 0, sd_coeffs = 0.5,
                        level = "row_id", nodes, steps, cycles, scale_factor = 1, ...) {

  if (class == "spline") {
    if (length(sd_coeffs) != 1 & length(sd_coeffs) != df) {
      stop("Please enter a vector of spline coefficient standard deviations equal to the spline degrees of freedom (df) or a single value.")
    }
    if (length(mu_coeffs) != 1 & .get_first_dim(mu_coeffs) != df) {
      stop("Please enter a vector of the mean spline coefficients equal to the number of spline degrees of freedom (df). Alternatively, provide a single value.")
    }

    coeffs <- matrix(rnorm(n * df, mu_coeffs,
                           sd_coeffs), nrow = df, ncol = n)

    traces <- map(1:n, \(i) .sim_splines(duration,df,
                                             coeffs[,i],...)) %>%
      list_rbind(names_to = level)

    return(list(traces = traces,
                coeffs = coeffs))

  } else if (class == "boxcar") {

    level_steps <- steps %*% t(scale_factor)

    traces <- map(1:n, \(i) .sim_boxcar(duration, nodes, level_steps[,i], cycles)) %>%
      list_rbind(names_to = level)

    return(list(traces = traces,
                scales = level_steps))

  } else if (class == "gp") {

    traces <- map(1:n, \(i) .sim_GP(duration, ...)) %>%
      list_rbind(names_to = level)

    return(traces)

  }

}

.sim_splines <- function(duration, df,
                         coeffs) {

  basis_matrix <- bs(seq(1,duration), df = df)
  y <- as.vector(basis_matrix %*% coeffs)

  return(tibble(x = seq(1,duration),
                y = y))
}

.sim_GP <- function(duration,kernel = "squared_exp", alpha = 0.5,
                    rho = duration / 5, nu = 3/2, P = duration / 2) {
  cov <- .kernel(duration, kernel, alpha, rho, nu, P)
  y <- mvrnorm(n = 1, mu = rep(0,duration), Sigma = cov)

  return(tibble(x = seq(1,duration),
                y = y))
}

.sim_boxcar <- function(duration, nodes, steps, cycles) {
  if (nodes[length(nodes)] < duration) {
    cycles <- 1
  }

  if (any(diff(nodes) / cycles < 1)) warning("Number of cycles smaller than step grid resolution")

  interval_id <- findInterval(1:(nodes[length(nodes)] / cycles),
                              nodes / cycles,rightmost.closed = TRUE) %>%
    rep_len(duration)

  return(tibble(x = 1:duration,
                y = steps[interval_id]))


}

.kernel <- function(duration, kernel, alpha, rho, nu, P) {
  X = seq(1,duration,by=1)

  if (kernel == "squared_exp") {
    K <- Gaussian.Kernel(rho)
    cov <- alpha^2 * Evaluate.Kernel(K,X)
  } else if (kernel == "matern") {
    K <- Matern.Kernel(rho, nu = nu)
    cov <- alpha^2 * Evaluate.Kernel(K,X)
  } else if (kernel == "rational_quad") {
    K <- RQ.Kernel(rho, alpha = 1)
    cov <- alpha^2 * Evaluate.Kernel(K,X)
  } else if (kernel == "periodic") {
    cov <- alpha^2 * .periodic_kernel(duration, rho, P)
  } else {
    stop("Please choose among squared_exp, matern, rational_quad and periodic")
  }
  return(cov)
}


.periodic_kernel <- function(duration, rho, P) {

  compute_periodic <- function(i,j,rho,P) {
    exp(-(2*sin(pi * abs(i-j)/P)^2)/rho^2)
  }
  cov_mat <- matrix(nrow = duration,
                    ncol = duration)
  for (i in 1:duration) {
    for (j in 1:duration) {
      cov_mat[i,j] = cov_mat[j,i] = compute_periodic(i,j,rho, P)
    }
  }
  return(cov_mat)
}

.get_first_dim <- function(x) {
  if (is.null(dim(x))) {
    return(length(x))
  } else {
    return(dim(x)[1])
  }
}




