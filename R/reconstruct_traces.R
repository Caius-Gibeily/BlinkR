

###
# should be modular and cumulative:
## plot features:
# posterior trajectories and bands. If simulation, then also show ground truth trajectories
# prior predictive checks: copies of each stan script but only with generated quantities
# option to select sample trajectories, event rasters/heatmaps, event rates and interevent distributions
#' @export
reconstruct_traces <- function(model,resolution=0.01,.width=c(0.5,0.8,0.99),...) {
  UseMethod("reconstruct_traces")
}

#' @export
#' @method reconstruct_traces one_ind
reconstruct_traces.one_ind <- function(model,.width,...) {

  M <- model$settings$M
  L <- model$settings$L_factor * model$settings$duration
  t_grid <- seq(0,model$settings$duration, length.out = 500)

  if (model$settings$kernel != "periodic") phi_basis <- .phi(t_grid, M, L)
  else phi_basis <- .phi_periodic(t_grid, M, model$settings$w0)

  n_grid <- nrow(phi_basis)
  post_data <- rstan::extract(model$fit)
  n_draws <- length(post_data$rho_group)

  I <- model$I
  G <- 1
  g_membership <- 1

  beta <- post_data$beta

  mu <- post_data$mu

  f <- phi_basis %*% t(beta) |> sweep(MARGIN = 2,
                                       STATS = mu, FUN = "+")


  gp_dat <- f

  tidy_quants <- .tidy_quantiles(t(gp_dat),t_grid,I,g_membership,n_grid,.width,...)
  class(tidy_quants) <- c("one_group","gp_model",class(tidy_quants))
  return(tidy_quants)
}

#' @export
#' @method reconstruct_traces one_group
reconstruct_traces.one_group <- function(model,.width,level="ind",...) {

  M <- model$settings$M
  L <- model$settings$L_factor * model$settings$duration
  t_grid <- seq(0,model$settings$duration, length.out = 500)

  if (model$settings$kernel != "periodic") { phi_basis <- .phi(t_grid, M, L)
  } else phi_basis <- .phi_periodic(t_grid, M, model$settings$w0)

  n_grid <- nrow(phi_basis)
  post_data <- rstan::extract(model$fit)
  n_draws <- length(post_data$rho_group)

  I <- model$I
  G <- 1
  g_membership <- model$g_membership
  beta_group <- post_data$beta_group
  beta_ind <- post_data$beta_ind |> rTensor::as.tensor()

  mu <- post_data$mu

  mu_ind <- post_data$mu_ind
  f_group <- phi_basis %*% t(beta_group)

  f_ind <- ttm(beta_ind,phi_basis,3) |> sweep(MARGIN = c(1,2),
                                               STATS = mu_ind, FUN = "+")
  f_group_ind <- f_ind |> sweep(MARGIN = c(1,3),
                                 STATS = t(f_group), FUN = "+")

  if (level == "ind") {
    gp_dat <- f_group_ind@data
    tidy_quants <- .tidy_quantiles(gp_dat,t_grid,I,g_membership,n_grid,level=level,.width = .width)
  } else if (level == "group") {
    gp_dat <- f_group |> sweep(MARGIN = 2, STATS = mu, FUN = "+")
    tidy_quants <- .tidy_quantiles(t(gp_dat),t_grid,G,g_membership,n_grid,level = level,.width = .width)
  }

  class(tidy_quants) <- c("one_group","gp_model",class(tidy_quants))
  return(tidy_quants)
}

#' @export
#' @method reconstruct_traces multi_group
reconstruct_traces.multi_group <- function(model,.width,level="ind",...) {

  M <- model$settings$M
  L <- model$settings$L_factor * model$settings$duration
  t_grid <- seq(0,model$settings$duration, length.out = 500)

  if (model$settings$kernel != "periodic") phi_basis <- .phi(t_grid, M, L)
  else phi_basis <- .phi_periodic(t_grid, M, model$settings$w0)

  n_grid <- nrow(phi_basis)
  post_data <- rstan::extract(model$fit)
  n_draws <- length(post_data$rho_global)

  I <- model$I
  G <- model$G
  g_membership <- model$g_membership
  gp_draws <- array(NA, dim = c(n_draws, I, n_grid))

  beta_global <- post_data$beta_global
  beta_group <- post_data$beta_group |> as.tensor()
  beta_ind <- post_data$beta_ind |> as.tensor()

  mu_group <- post_data$mu_group
  mu_ind <- post_data$mu_ind

  f_global <- phi_basis %*% t(beta_global)
  f_group <- ttm(beta_group,phi_basis,3) #|> sweep(MARGIN = c(1,2),
                                          #          STATS = mu_group, FUN = "+")
  f_ind <- ttm(beta_ind,phi_basis,3) |> sweep(MARGIN = c(1,2),
                                              STATS = mu_ind, FUN = "+")
  f_group_ind <- (f_group[,g_membership,] + f_ind) |> sweep(MARGIN = c(1,3),
                                                           STATS = t(f_global), FUN = "+")
  if (level == "ind") {
    gp_dat <- f_group_ind@data
    tidy_quants <- .tidy_quantiles(gp_dat,t_grid,I,g_membership,n_grid,level=level,.width = .width)
  } else if (level == "group") {
    f_group <- f_group |> sweep(MARGIN = c(1,2),
                                STATS = mu_group, FUN = "+") |>
      sweep(MARGIN = c(1,3),
            STATS = t(f_global), FUN = "+")
    gp_dat <- f_group@data
    tidy_quants <- .tidy_quantiles(gp_dat,t_grid,G,g_membership,n_grid,level=level,.width = .width)
  } else if (level == "global") {
    gp_dat <- f_global |> sweep(MARGIN = 2, STATS = mu, FUN = "+")
    tidy_quants <- .tidy_quantiles(t(gp_dat),t_grid,1,g_membership,n_grid,level=level,.width = .width)
  }

  class(tidy_quants) <- c("one_group","gp_model",class(tidy_quants))
  return(tidy_quants)
}


.phi <- function(x, M, L) {
  outer(x,seq_len(M), function(x, m)
    sin(pi * m * (x + L) / (2 * L)) / sqrt(L)
  )
}

.phi_periodic <- function(x,M,w0) {
  k <- seq_len(M / 2)
  w0xk <- outer(w0 * x, k)
  return(cbind(cos(w0xk),sin(w0xk)))
}


#' @noRd
.tidy_quantiles <- function(f_data, t_grid, N, g_membership, n_grid, .width = c(0.5, 0.8, 0.99), level = "ind") {

  dims <- dim(f_data)
  gp_flattened <- matrix(f_data, nrow = dims[1], ncol = prod(dims[-1]))

  base_df <- switch(
    level,
    "ind" = tibble(
      ind = as.factor(rep(1:N, times = n_grid)),
      group = as.factor(g_membership[ind]),
      x = rep(t_grid, each = N)
    ), "group" = tibble(
      group = rep(1:N, times = length(t_grid)),
      x = rep(t_grid, each = N)
    ), "global" = tibble(
      global = 1,
      x = t_grid
    ),
    stop("Invalid level: ", level)
  )

  probs <- sort(unique(c(0.5, (1 - .width) / 2, 1 - (1 - .width) / 2)))
  quants <- matrixStats::colQuantiles(gp_flattened, probs = probs)
  quantile_df <- bind_cols(base_df, as_tibble(quants))

  id_cols <- names(base_df)
  medians_df <- quantile_df %>%
    select(all_of(id_cols), y = `50%`)

  long_df <- map_dfr(.width, function(w) {
    low_col  <- paste0(round((1 - w) / 2 * 100, 1), "%")
    high_col <- paste0(round((1 - (1 - w) / 2) * 100, 1), "%")

    quantile_df %>%
      select(all_of(id_cols), ymin = !!sym(low_col), ymax = !!sym(high_col)) %>%
      mutate(.width = w)
  })

  left_join(long_df, medians_df, by = id_cols)
}

#' @noRd
.diagSPD_EQ <- function(alpha, rho, M, L) {
  indices <- seq_len(M)
  factor <- alpha * sqrt(sqrt(2 * pi) * rho)
  exponent <- -0.25 * (rho * pi / (2 * L))^2
  factor * exp(exponent * indices^2)
}

.diagSPD_Matern52 <- function(alpha, rho, M, L) {
  factor <- 16 * (sqrt(5) / rho)^5
  indices <- (pi / (2 * L) * seq_len(M))^2

  denom <- 3 * ((5 / rho^2) + indices)^3
  return(alpha * sqrt(factor / denom))
}

#' @noRd
.diagSPD_Matern12 <- function(alpha, rho, M, L) {
  indices <- 1:M
  factor <- 2.0

  denom <- rho * ((1.0 / rho)^2 + (pi * indices / (2 * L))^2)
  return(alpha * sqrt(factor * (1/denom)))
}

#' @noRd
.diagSPD_Matern32 <- function(alpha, rho, M, L) {

  M_series <- seq_len(M)
  indices <- ((pi / (2 * L)) * M_series)^2

  factor <- 2 * alpha * (sqrt(3) / rho)^1.5
  denom <- 3 / (rho^2) + indices
  return(factor / denom)
}

#' @noRd
.diagSPD_Periodic <- function(alpha, rho, M,...) {
  a <- 1 / (rho^2)
  indices <- 1:M

  log_bessel <- log(besselI(x = a, nu = indices, expon.scaled = TRUE)) + a

  q <- exp(log(alpha) + 0.5 * (log(2) - a + log_bessel))

  # append_row(q, q) is equivalent to concatenating the vector with itself
  return(append_row(q,q))
}
