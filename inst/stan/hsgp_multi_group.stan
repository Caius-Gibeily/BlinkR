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
    vector[M] denom = 1 / rho + rho * basis_indices(M, L);
    return alpha * sqrt(2 ./ denom);
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
    return q;
  }

  // Matern squared exponential
  vector diagSPD_EQ(real alpha, real rho, int M, real L) {
    vector[M] indices = linspaced_vector(M, 1, M);
    real factor = alpha * sqrt(sqrt(2 * pi()) * rho);
    real exponent = -0.25 * (rho * pi() / 2 / L)^2;
    return factor * exp(exponent * square(indices));
  }

  real gengamma_lpdf(real x, real shape, real k, real scale) {

    if (x <= 0.0) {
      return negative_infinity();
    }
    real logdens = log(shape) - lgamma(k) + (shape * k - 1.0) * log(x) - (shape * k) * log(scale) - pow(x/scale,shape);
    return logdens;

  }

  matrix sum_to_zero_groups(int G, int I, int M, array[] int g_membership, matrix z_ind_raw, array[] int I_per_group) {
    matrix[I, M] z_ind;

    int pos = 1;
    for (g in 1:G) {

      int n = I_per_group[g];

      int idx = 0;
      int last_i = 0;


      for (i in 1:I) {
        if (g_membership[i] == g) {
          idx += 1;
          if (idx < n) {
            z_ind[i, ] = z_ind_raw[pos, ];
            pos += 1;
          } else {
            last_i = i;
          }
        }
      }

      // Constrained row
      for (m in 1:M) {
        real s = 0;

        for (i in 1:I) {
          if (g_membership[i] == g && i != last_i)
          s += z_ind[i, m];
        }

        z_ind[last_i, m] = -s;
      }
    }

    return z_ind;
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
  int<lower=1, upper=5> kernel;
  int<lower=1, upper=5> family;

  real<lower=0> L_factor;

  int<lower=0, upper=N_params> include_k;
  int<lower=0, upper=N_params> include_shape;
  int<lower=0, upper=N_params> include_sigma_lognormal;
  // data
  int<lower=1> N_total;
  int<lower=1> I;

  array[N_total] int<lower=1,upper=I> ind_id;

  int<lower=1> G;
  array[N_total] int<lower=1,upper=G> g_id;
  array[I] int<lower=1,upper=G> g_membership;
  array[G] int I_per_group;

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
  vector[M] z_global;
  matrix[G-1,M] z_group_raw;
  matrix[I-G,M] z_ind_raw;

  real<lower=0> rho_global;
  real<lower=0> alpha_global;

  real<lower=0> rho_group;
  real<lower=0> alpha_group;

  // Individual-level
  real<lower=0> rho_ind;
  real<lower=0> alpha_ind;


  // indect intercept
  vector[G] mu_group;
  vector[I] mu_ind;

  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {

  matrix[I,M] z_ind;
  matrix[G,M] z_group;

  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S_global;
  vector[M] diag_S_group;
  vector[M] diag_S_ind;

  if (kernel == 1) {
    diag_S_global = diagSPD_EQ(alpha_global,rho_global,M,L);
    diag_S_group = diagSPD_EQ(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_EQ(alpha_ind,rho_ind,M,L);
  } else if (kernel == 2) {
    diag_S_global = diagSPD_Matern12(alpha_global,rho_global,M,L);
    diag_S_group = diagSPD_Matern12(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_Matern12(alpha_ind,rho_ind,M,L);
  } else if (kernel == 3) {
    diag_S_global = diagSPD_Matern32(alpha_global,rho_global,M,L);
    diag_S_group = diagSPD_Matern32(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_Matern32(alpha_ind,rho_ind,M,L);
  } else if (kernel == 4) {
    diag_S_global = diagSPD_Matern52(alpha_global,rho_global,M,L);
    diag_S_group = diagSPD_Matern52(alpha_group,rho_group,M,L);
    diag_S_ind = diagSPD_Matern52(alpha_ind,rho_ind,M,L);
  } else if (kernel == 5) {
    diag_S_global = diagSPD_Periodic(alpha_global,rho_global,M);
    diag_S_group = diagSPD_Periodic(alpha_group,rho_group,M);
    diag_S_ind = diagSPD_Periodic(alpha_ind,rho_ind,M);
  }

  z_group[1:(G - 1), ] = z_group_raw;

  for (m in 1:M) {
    z_group[G, m] = -sum(z_group_raw[, m]);
  }

  z_ind = sum_to_zero_groups(G, I, M, g_membership, z_ind_raw, I_per_group);

  vector[M] beta_global = diag_S_global .* z_global;

  matrix[G, M] beta_group = diag_post_multiply(z_group,diag_S_group);

  matrix[I, M] beta_ind = diag_post_multiply(z_ind,diag_S_ind);

}

model {

  // Priors
  z_global ~ std_normal();
  to_vector(z_group) ~ std_normal();
  to_vector(z_ind) ~ std_normal();

  if (distributions[1] == 1) {
    mu_group ~ normal(params[1,1],params[1,2]);
    mu_ind ~ normal(params[2,1],params[2,2]);
  } else if (distributions[1] == 2) {
    mu_group ~ lognormal(params[1,1],params[1,2]);
    mu_ind ~ lognormal(params[2,1],params[2,2]);
  } else if (distributions[1] == 3) {
    mu_group ~ cauchy(params[1,1],params[1,2]);
    mu_ind ~ cauchy(params[2,1],params[2,2]);
  } else if (distributions[1] == 4) {
    mu_group ~ exponential(params[1,1]) T[1, ];
    mu_ind ~ exponential(params[2,1]) T[1, ];
  }

  apply_prior_lp(alpha_global, distributions[3], params[3, 1], params[3, 2]);
  apply_prior_lp(alpha_group, distributions[4], params[4, 1], params[4, 2]);
  apply_prior_lp(alpha_ind, distributions[5], params[5, 1], params[5, 2]);

  apply_prior_lp(rho_global, distributions[6], params[6, 1], params[6, 2]);
  apply_prior_lp(rho_group, distributions[7], params[7, 1], params[7, 2]);
  apply_prior_lp(rho_ind, distributions[8], params[8, 1], params[8, 2]);


  if (include_k != 0) apply_prior_lp(k[1], distributions[include_k],
    params[include_k, 1], params[include_k, 2]);
  if (include_shape != 0) apply_prior_lp(shape[1], distributions[include_shape],
    params[include_shape, 1], params[include_shape, 2]);
  if (include_sigma_lognormal != 0) apply_prior_lp(sigma_lognormal[1], distributions[include_sigma_lognormal],
    params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2]);

  array[3] vector[N_total] eta_quad;

  for (j in 1:3) {
    vector[N_total] f_global = PHI_quad[j] * beta_global;
    vector[N_total] f_group = rows_dot_product(PHI_quad[j], beta_group[g_id, ]);
    vector[N_total] f_ind = rows_dot_product(PHI_quad[j], beta_ind[ind_id, ]);


    eta_quad[j] = mu_group[g_id] + mu_ind[ind_id] + f_global + f_group + f_ind;
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
