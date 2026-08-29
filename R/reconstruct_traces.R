

###
# should be modular and cumulative:
## plot features:
# posterior trajectories and bands. If simulation, then also show ground truth trajectories
# prior predictive checks: copies of each stan script but only with generated quantities
# option to select sample trajectories, event rasters/heatmaps, event rates and interevent distributions

#' @export
reconstruct_traces <- function(model, level = c("ind","group","global"), prior = FALSE, dev_only = FALSE, resolution = 0.1) {
  if (dev_only && "one_ind" %in% class(model)) {
    warning("Computing deviations is not possible in a single-individual model with one hierarchical level. Defaulting to individual GP trace")
    dev_only = FALSE
  }

  if (length(level) > 1) level = "ind"

  if (level %in% c("group","global") & "one_ind" %in% class(model)) {
    stop("Only the individual level is available for one-individual models")
  } else if (level == "global" & "one_group" %in% class(model)) {
    stop("only the individual and group levels are available for one-group models")
  }


  # Basis parameters
  M <- model$settings$M
  L <- model$settings$L_factor * model$settings$duration

  t_grid <- seq(0,model$settings$duration,
                by = resolution)
  if (model$settings$kernel != "periodic") {
    phi_basis <- .phi(t_grid, M, L)
  } else {
    phi_basis <- .phi_periodic(t_grid, M,
                               model$settings$w0)
  }

  # extract data - choose whether to draw from prior (if generated) or posterior
  n_grid <- nrow(phi_basis)
  if (!prior & "gp_model" %in% class(model)) {
    post_data <- rstan::extract(model$fit)
  } else if (prior & !is.null(model$prior_pc) | "gp_prior_pc" %in% class(model)) {
    post_data <- rstan::extract(model$prior_pc)
  } else if (prior & is.null(model$prior_pc)) {
    stop("Please ensure your model has prior data. You may do so by setting run = 'prior_pc' or c('prior_pc','fit') in gp_fit()")
  }

  ## Common to all gp_models
  I <- model$I
  G <- model$G
  g_membership <- model$g_membership

  N <- switch(level,
    "ind" = I,
    "group" = G,
    "global" = 1
  )

  recon_data <- list(t_grid = t_grid,
                     N = N,
                     g_membership = g_membership,
                     n_grid = n_grid,
                     level = level)

  survival_params <- switch(model$settings$family,
                            "exponential" = NULL,
                            "gamma" = list(k = post_data$k),
                            "weibull" = list(shape = post_data$shape),
                            "gengamma" = list(k = post_data$k,
                                              shape = post_data$shape),
                            "lognormal" = list(post_data$sigma_lognormal))

  recon_data$survival_params <- survival_params

  beta_ind <- post_data$beta_ind
  if (length(dim(beta_ind)) == 2) {
    beta_ind <- t(beta_ind) |> .add_dim()
  }

  beta_ind <- beta_ind |> as.tensor()

  mu_ind <- post_data$mu_ind

  #> sweep(MARGIN = c(1,2),
  #          STATS = mu_group, FUN = "+")
  f_ind <- ttm(beta_ind,phi_basis,3) |> sweep(MARGIN = c(1,2),
                                              STATS = mu_ind, FUN = "+")

  if (any(c("one_group","multi_group") %in% class(model))) {

    beta_group <- post_data$beta_group

    if (length(dim(beta_group)) == 2) {
      beta_group <- .add_dim(beta_group)

    }

    beta_group <- as.tensor(beta_group)


    mu_group <- post_data$mu_group

    f_group <- ttm(beta_group,phi_basis,3) #|> sweep(MARGIN = c(1,2),
    #          STATS = mu_group, FUN = "+")



    if ("multi_group" %in% class(model)) {
      beta_global <- post_data$beta_global
      mu_global <- post_data$mu_global
      f_global <- phi_basis %*% t(beta_global)

      f_group_ind <- (f_group[,g_membership,] + f_ind) |>
        sweep(MARGIN = c(1,3),
              STATS = t(f_global), FUN = "+")
    } else {
      f_group_ind <- f_ind |> sweep(MARGIN = c(1,3),
                                    STATS = f_group@data, FUN = "+")
    }
  }

  if (level == "ind") {
    if (!dev_only) gp_dat <- f_group_ind@data
    else gp_dat <- f_ind@data
    recon_data$gp_dat <- gp_dat

    return(recon_data)
    #tidy_quants <- .tidy_quantiles(gp_dat,t_grid,I,g_membership,n_grid,level=level,.width = .width)
  } else if (level == "group") {


    if (!dev_only) {
      if ("one_group" %in% class(model)) {
        f_group <- f_group |> sweep(MARGIN = 1,STATS = mu_group, FUN = "+")
      }
      else {
        f_group <- f_group |>
          sweep(MARGIN = c(1,2),STATS = mu_group, FUN = "+") |>
          sweep(MARGIN = c(1,3), STATS = t(f_global), FUN = "+")
      }
      gp_dat <- f_group@data

    } else {
      gp_dat <- f_group@data
    }
    #tidy_quants <- .tidy_quantiles(gp_dat,t_grid,G,g_membership,
    #                               n_grid,level=level,.width = .width)
    recon_data$gp_dat <- gp_dat
    return(recon_data)

  } else if (level == "global") {
    if (dev_only) {
      warning("Computing deviations only can only be performed at lower hierarchy levels (ind, group). Defaulting to global GP trace")
    }

    gp_dat <- f_global |>
      sweep(MARGIN = 2, STATS = mu_global, FUN = "+") |> t()
    recon_data$gp_dat <- recon_data
    return(recon_data)
    #tidy_quants <- .tidy_quantiles(t(gp_dat),t_grid,1,g_membership,
    #                               n_grid,level=level,.width = .width)


  } else stop("Please choose a level from 'ind', 'group' and 'global'")

}


#' @export
tidy_traces <- function(model = NULL, recon_data = NULL, level = c("ind","group","global"),
                        prior = FALSE, dev_only = FALSE, resolution = 0.2, .width=c(0.5,0.8,0.99), rescale = FALSE) {
  if (is.null(model) & is.null(recon_data)) {
    stop("Please ensure you supply either a fitted gp_model or reconstructed traces from a previous gp_model.")
  }

  if (is.null(recon_data)) {
    recon_data <- reconstruct_traces(model, level, prior, dev_only, resolution)
  }
  recon_data$.width <- .width

  if (rescale) {
    recon_data$gp_dat <- .apply_rescaling(model,recon_data)
  }

  recon_data <- recon_data[!names(recon_data) %in% "survival_params"]


  tidy_quants <- do.call(.tidy_quantiles,recon_data)
  return(tidy_quants)
}

#' @export
draw_traces <- function(model = NULL,recon_data = NULL, level="ind",prior=FALSE, dev_only = FALSE,
                            resolution = 0.1, .width=c(0.5,0.8,0.99), n_samples = 1000, rescale = FALSE) {
  if (is.null(model) & is.null(recon_data)) {
    stop("Please ensure you supply either a fitted gp_model or reconstructed traces from a previous gp_model.")
  }
  if (is.null(recon_data)) {
    recon_data <- reconstruct_traces(model, level, prior, dev_only, resolution)
  }

  if (rescale) {
    recon_data$gp_dat <- .apply_rescaling(model,recon_data)
  }

  draw_params <- list(.width = .width,
                      n_samples = n_samples)

  recon_data <- append(recon_data,draw_params)
  recon_data <- recon_data[!(names(recon_data) %in% c("level",".width"))]
  tidy_samples <- do.call(.tidy_samples,recon_data)

  return(tidy_samples)
}

#' @export
summarise_traces <- function(model = NULL,recon_data = NULL, level="ind",prior=FALSE, dev_only = FALSE,
                             resolution = 0.1, rescale = FALSE) {
  if (is.null(model) & is.null(recon_data)) {
    stop("Please ensure you supply either a fitted gp_model or reconstructed traces from a previous gp_model.")
  }
  if (is.null(recon_data)) {
    recon_data <- reconstruct_traces(model, level, prior, dev_only, resolution)
  }


  if (rescale) {
    recon_data$gp_dat <- .apply_rescaling(model,recon_data)
  }

  recon_data <- recon_data[!names(recon_data) %in% "survival_params"]
  summary_traces <- do.call(.tidy_mean_sd, recon_data)

  return(summary_traces)
}

#' @noRd
.add_dim <- function(beta_matrix) {
  dim(beta_matrix) <- c(nrow(beta_matrix),1,ncol(beta_matrix))
  return(beta_matrix)
}


#' @noRd
.tidy_quantiles <- function(gp_dat, t_grid, N, g_membership, n_grid,
                            .width = c(0.5, 0.8, 0.99), level = "ind") {

  dims <- dim(gp_dat)
  gp_flattened <- matrix(gp_dat, nrow = dims[1], ncol = prod(dims[-1]))

  base_df <- switch(
    level,
    "ind" = tibble(
      ind   = as.factor(rep(1:N, times = n_grid)),
      group = as.factor(g_membership[ind]),
      x     = rep(t_grid, each = N)
    ),
    "group" = tibble(
      group = as.factor(rep(1:N, times = length(t_grid))),
      x     = rep(t_grid, each = N)
    ),
    "global" = tibble(
      global = 1,
      x      = t_grid
    ),
    stop("Invalid level: ", level)
  )

  long_df <- map(.width, \(w) {
    probs <- sort(c((1 - w) / 2, 0.5, 1 - (1 - w) / 2))
    quants <- matrixStats::colQuantiles(gp_flattened, probs = probs)

    quants_tibble <- tibble(ymin = quants[, 1],
                            y = quants[, 2],
                            ymax = quants[, 3])

    bind_cols(base_df, quants_tibble) %>%
      mutate(.width = w)
  }) %>%
    list_rbind()

  return(long_df)
}

#' @noRd
.tidy_mean_sd <- function(gp_dat, t_grid, N, g_membership, n_grid,
                            level = "ind") {

  dims <- dim(gp_dat)
  gp_flattened <- matrix(gp_dat, nrow = dims[1], ncol = prod(dims[-1]))

  base_df <- switch(
    level,
    "ind" = tibble(
      ind = as.factor(rep(1:N, times = n_grid)),
      group = as.factor(g_membership[ind]),
      x = rep(t_grid, each = N)),
    "group" = tibble(group = rep(1:N, times = length(t_grid)),
      x = rep(t_grid, each = N)),
    "global" = tibble(global = 1,
      x = t_grid),
    stop("Invalid level: ", level)
  )

  col_means <- tibble(means=matrixStats::colMeans2(gp_flattened))
  col_sds <- tibble(sds = matrixStats::colSds(gp_flattened))

  mean_sd_df <- bind_cols(base_df, col_means) |> bind_cols(col_sds)


  return(mean_sd_df)
}

#' @noRd
.tidy_samples <- function(gp_dat, survival_params, t_grid, N, g_membership, n_grid, n_samples = 100) {
  if (is.null(n_samples)) {
    n_samples <- 100
  }
  dims <- dim(gp_dat)

  if (n_samples > dims[1]) {
    stop("Please choose a number of samples no greater than the number of post-warmup samples.")
  }

  sample_indices <- sample(dims[1], size = n_samples)

  sampled_gp <- gp_dat[sample_indices, 1:N, 1:n_grid]

  if (length(dim(sampled_gp)) == 2) {
    dim(sampled_gp) <- c(nrow(sampled_gp),1,ncol(sampled_gp))
  }
  y_vec <- as.vector(aperm(sampled_gp,c(3,2,1)))

  base_df <- tidyr::expand_grid(
    sample = as.factor(1:n_samples),
    ind = as.factor(1:N),
    grid_idx = seq_along(t_grid)
  ) |>
    mutate(
      x = t_grid[grid_idx],
      group = as.factor(g_membership[ind])
    ) |>
    select(ind, group, sample, x)

  sampled_traces <- bind_cols(base_df, y = y_vec)

  return(list(sampled_traces = sampled_traces,
              survival_params = lapply(survival_params, function(mat) mat[sample_indices])))
}
#' @noRd
.apply_rescaling <- function(model, recon_data) {
  if (model$settings$family == "gamma") {
    rescaled <- recon_data$gp_dat |>
      sweep(MARGIN = 1,
            STATS = log(recon_data$survival_params$k),FUN = "+")
  } else if (model$settings$family == "weibull") {

    rescaled <- recon_data$gp_dat |>
      sweep(MARGIN = 1,
            STATS = log(gamma(1 + (1/recon_data$survival_params$shape))),
            FUN = "+")
  } else if (model$settings$family == "gengamma") {
    rescaled <- recon_data$gp_dat |>
      sweep(MARGIN = 1,
            STATS = log(gamma((recon_data$survival_params$shape+1)/
                                recon_data$survival_params$k)/
                          gamma((recon_data$survival_params$shape)/
                                  recon_data$survival_params$k)),
            FUN = "+")
  }

  return(rescaled)

}

#' @noRd
.phi <- function(x, M, L) {
  outer(x,seq_len(M), function(x, m)
    sin(pi * m * (x + L) / (2 * L)) / sqrt(L)
  )
}

#' @noRd
.phi_periodic <- function(x,M,w0) {
  k <- seq_len(M / 2)
  w0xk <- outer(w0 * x, k)
  return(cbind(cos(w0xk),sin(w0xk)))
}

#' @noRd
.diagSPD_EQ <- function(alpha, rho, M, L) {
  indices <- seq_len(M)
  factor <- alpha * sqrt(sqrt(2 * pi) * rho)
  exponent <- -0.25 * (rho * pi / (2 * L))^2
  factor * exp(exponent * indices^2)
}

#' @noRd
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

