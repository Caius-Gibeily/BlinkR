functions {

  // spectral basis set, phi
  matrix phi(int N, int M, real L, vector x) {
    matrix[N, M] res;

    for (m in 1:M) {
      res[, m] = inv_sqrt(L) * sin(pi() * m * (x + L) / (2 * L));
    }

    return res;
  }

  // Basis indices
  vector basis_indices(int M,real L) {
    vector[M] indices = linspaced_vector(M, 1, M);

    return square(pi() / (2 * L) * indices);
  }

  // Matern 1/2 - sourced from
  // https://epiforecasts.io/EpiNow2/stan/gaussian__process_8stan_source.html#l00090
  vector diagSPD_Matern12(real alpha, real rho, int M, real L) {

    vector[M] denom =
      rho * (square(inv(rho)) + square((pi() / (2 * L)) * basis_indices(M, L)));

    return alpha * sqrt(2.0 ./ denom);
  }

  // Matern 3/2
  vector diagSPD_Matern32(real alpha, real rho, int M, real L) {

    real factor = 2 * alpha * (sqrt(3) / rho)^1.5;
    vector[M] denom = 3 / square(rho) + basis_indices(M, L);

    return factor ./ denom;
  }
  // Matern 5/2
  vector diagSPD_Matern52(real alpha, real rho, int M, real L) {
    real factor = 16 * pow(sqrt(5) / rho, 5);
    vector[M] denom = 3 * pow(5 / square(rho) + basis_indices(M, L), 3);
    return alpha * sqrt(factor ./ denom);
  }

  // Periodic
  vector diagSPD_Periodic(real alpha, real rho, int M) {
    real a = inv_square(rho);
    vector[M] indices = linspaced_vector(M, 1, M);
    vector[M] q = exp(
      log(alpha) + 0.5 *
        (log(2) - a + to_vector(log_modified_bessel_first_kind(indices, a)))
    );
    return append_row(q, q);
  }

  // Matern squared exponential
  vector diagSPD_EQ(real alpha, real rho, int M, real L) {

    real factor = alpha * sqrt(sqrt(2 * pi()) * rho);
    real exponent = -0.25 * (rho * pi() / 2 / L)^2;

    return factor * exp(exponent * square(basis_indices(M, L)));
  }

  real gengamma_lpdf(real x, real shape, real k, real scale) {

    if (x <= 0.0) {
      return negative_infinity();
    }
    real logdens = log(shape) - lgamma(k) + (shape * k - 1.0) * log(x) - (shape * k) * log(scale) - pow(x/scale,shape);
    return logdens;

  }

  void apply_prior_lp(real param, real dist, real arg1, real arg2) {
    if (dist == 1) {
      target += normal_lpdf(param | arg1, arg2);
    } else if (dist == 2) {
      target += lognormal_lpdf(param | arg1, arg2);
    } else if (dist == 3) {
      target += cauchy_lpdf(param | arg1, arg2);
    } else if (dist == 4) {
      target += exponential_lpdf(param | arg1) - exponential_lccdf(1.0 | arg1);
    }
  }

}


data {
  // model settings
  int<lower=1> M;
  int<lower=1> N_params;
  array[N_params,2] real params;
  vector[N_params] distributions;
  int<lower=1, upper=4> kernel;
  int<lower=1, upper=5> family;

  real<lower=0> L_factor;

  int<lower=0, upper=N_params> include_k;
  int<lower=0, upper=N_params> include_shape;
  int<lower=0, upper=N_params> include_sigma_lognormal;
  // data
  int<lower=1> N_total;
  int<lower=1> I;
  array[N_total] int<lower=1,upper=I> ind_id;

  vector[N_total] t_ev;
  vector[N_total] dt;


}

transformed data {

  real L = L_factor * max(t_ev);

  // 3-point Gauss-Legendre quadrature
  vector[3] qx =
    [0.1127016654,
     0.5,
     0.8872983346]';

  vector[3] qw =
    [5.0 / 18.0,
     8.0 / 18.0,
     5.0 / 18.0]';

  vector[3] log_qw = log(qw);

  // precompute basis matrices
  array[3] matrix[N_total, M] PHI_quad;

  for (j in 1:3) {

    vector[N_total] t_quad = t_ev - (1.0 - qx[j]) .* dt;

    PHI_quad[j] = phi(N_total,M,L,t_quad);
  }
}

parameters {

  // population GP
  vector[M] z_group;
  matrix[I-1,M] z_ind_raw;

  real<lower=0> rho_group;
  real<lower=0> alpha_group;

  // Individual-level
  real<lower=0> rho_ind;
  real<lower=0> alpha_ind;


  // indect intercept
  vector[I] mu;

  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {

  matrix[I,M] z_ind;

  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S_group;
  vector[M] diag_S_ind;

  if (kernel == 1) {
    diag_S_group = diagSPD_EQ(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_EQ(alpha_ind,rho_ind,M,L);
  } else if (kernel == 2) {
    diag_S_group = diagSPD_Matern12(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_Matern12(alpha_ind,rho_ind,M,L);
  } else if (kernel == 3) {
    diag_S_group = diagSPD_Matern32(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_Matern32(alpha_ind,rho_ind,M,L);
  } else if (kernel == 4) {
    diag_S_group = diagSPD_Matern52(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_Matern52(alpha_ind,rho_ind,M,L);
  } else if (kernel == 5) {
    diag_S_group = diagSPD_Periodic(alpha_group,rho_group,M);
    diag_S_ind = diagSPD_Periodic(alpha_ind,rho_ind,M);
  }

  z_ind[1:(I - 1), ] = z_ind_raw;

  for (m in 1:M) {
    z_ind[I, m] = -sum(z_ind_raw[, m]);
  }

  vector[M] beta_group = diag_S_group .* z_group;

  matrix[I, M] beta_ind = diag_post_multiply(z_ind,diag_S_ind);

}

model {

  // Priors
  z_group ~ std_normal();
  to_vector(z_ind) ~ std_normal();

  if (distributions[1] == 1) {
    mu ~ normal(params[1,1],params[1,2]);
  } else if (distributions[1] == 2) {
    mu ~ lognormal(params[1,1],params[1,2]);
  } else if (distributions[1] == 3) {
    mu ~ cauchy(params[1,1],params[1,2]);
  } else if (distributions[1] == 4) {
    mu ~ exponential(params[1,1]) T[1, ];
  }

  apply_prior_lp(alpha_group, distributions[2], params[2, 1], params[2, 2]);
  apply_prior_lp(alpha_ind, distributions[3], params[3, 1], params[3, 2]);

  apply_prior_lp(rho_group, distributions[4], params[4, 1], params[4, 2]);
  apply_prior_lp(rho_ind, distributions[5], params[5, 1], params[5, 2]);


  if (include_k != 0) apply_prior_lp(k[1], distributions[include_k],
  params[include_k, 1], params[include_k, 2]);
  if (include_shape != 0) apply_prior_lp(shape[1], distributions[include_shape],
  params[include_shape, 1], params[include_shape, 2]);
  if (include_sigma_lognormal != 0) apply_prior_lp(sigma_lognormal[1], distributions[include_sigma_lognormal],
  params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2]);


  array[3] vector[N_total] eta_quad;

  for (j in 1:3) {
    vector[N_total] f_group = PHI_quad[j] * beta_group;

    vector[N_total] f_ind = rows_dot_product(PHI_quad[j], beta_ind[ind_id, ]);


    eta_quad[j] = mu[ind_id] + f_group + f_ind;
  }



  matrix[N_total, 3] log_kernel;

  if (family == 1) {
    for (j in 1:3) {
      vector[N_total] scale = exp(eta_quad[j]);
      for (n in 1:N_total) {
        log_kernel[n, j] = log_qw[j] + exponential_lpdf(dt[n] | inv(scale[n]));
      }
    }
  } else if (family == 2) {
    for (j in 1:3) {
      vector[N_total] scale = exp(eta_quad[j]);
      for (n in 1:N_total) {
        log_kernel[n, j] = log_qw[j] + gamma_lpdf(dt[n] | k[1], inv(scale[n]));
      }
    }
  } else if (family == 3) {
    for (j in 1:3) {
      vector[N_total] scale = exp(eta_quad[j]);
      for (n in 1:N_total) {
        log_kernel[n, j] = log_qw[j] + weibull_lpdf(dt[n] | shape[1], scale[n]);
      }
    }
  } else if (family == 4) {
    for (j in 1:3) {
      vector[N_total] mu_lognormal = eta_quad[j];
      for (n in 1:N_total) {
        log_kernel[n, j] = log_qw[j] + lognormal_lpdf(dt[n] | mu_lognormal[n], sigma_lognormal[1]);
      }
    }
  } else if (family == 5) {
    for (j in 1:3) {
      vector[N_total] scale = exp(eta_quad[j]);
      for (n in 1:N_total) {
        log_kernel[n, j] = log_qw[j] + gengamma_lpdf(dt[n] | shape[1], k[1], scale[n]);
      }
    }
  }

  for (n in 1:N_total) {
    target += log_sum_exp(log_kernel[n, ]);
  }
}
