functions {

  matrix phi(int N, int M, real L, vector x) {
    matrix[N, M] res;

    for (m in 1:M) {
      res[, m] = inv_sqrt(L) * sin(pi() * m * (x + L) / (2 * L));
    }

    return res;
  }

  vector diagSPD_Matern12(
      real alpha,
      real rho,
      int M,
      real L
  ) {

    vector[M] indices = linspaced_vector(M, 1, M);

    vector[M] denom =
      rho * (square(inv(rho)) + square((pi() / (2 * L)) * indices));

    return alpha * sqrt(2.0 * inv(denom));
  }
  vector matern_indices(int M, 
  real L
  ) {
    vector[M] indices = linspaced_vector(M, 1, M);
    
    return square(pi() / (2 * L) * indices);
  }
  vector diagSPD_Matern32(real alpha, 
    real rho, 
    int M, 
    real L
    ) {
    real factor = 2 * alpha * (sqrt(3) / rho)^1.5;
    vector[M] denom = 3 / square(rho) + matern_indices(M, L);
  return factor ./ denom;
  }
  
  vector diagSPD_EQ(real alpha, real rho, int M, real L) {
    vector[M] indices = linspaced_vector(M, 1, M);
    real factor = alpha * sqrt(sqrt(2 * pi()) * rho);
    real exponent = -0.25 * (rho * pi() / 2 / L)^2;
    return factor * exp(exponent * square(indices));
}
 
}

data {

  int<lower=1> S;
  int<lower=1> N_total;
  int<lower=1> M;
  real<lower=0> mu_rho_pop;
  real<lower=0> mu_rho_subj;
  real<lower=0> L_factor;

  vector[N_total] t_ev;
  vector[N_total] dt;

  array[N_total]
    int<lower=1, upper=S> s_id;
}

transformed data {

  real L =
    L_factor * max(t_ev);

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

    PHI_quad[j] =
      phi(
        N_total,
        M,
        L,
        t_quad
      );
  }
}

parameters {

  // population GP
  vector[M] z_pop;
  real<lower=0> rho_pop;
  real<lower=0> alpha_pop;

  // subject deviations
  matrix[S - 1, M] z_subj_raw;

  real<lower=0> rho_subj;
  real<lower=0> alpha_subj;

  // subject intercepts
  vector[S] mu;

  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;
  real<lower=1> k;
  real<lower=0> lambda_shape;
}

transformed parameters {

  //real k = 1.0 + exp(log_k_minus1);

  matrix[S, M] z_subj;

  vector[M] diag_S_pop = diagSPD_EQ(alpha_pop,rho_pop,M,L);

  vector[M] diag_S_subj = diagSPD_EQ(alpha_subj,rho_subj,M,L);

  // zero-sum subject constraint
  z_subj[1:(S - 1), ] = z_subj_raw;

  for (m in 1:M) {
    z_subj[S, m] = -sum(z_subj_raw[, m]);
  }

  vector[M] spectral_pop = diag_S_pop .* z_pop;

  matrix[S, M] beta_subj = diag_post_multiply(z_subj,diag_S_subj);
}

model {

  // Priors

  z_pop ~ std_normal();

  to_vector(z_subj_raw) ~ normal(0,sqrt(1.0 - 1.0 / S));

  rho_pop ~ lognormal(mu_rho_pop,0.8);

  rho_subj ~ lognormal(mu_rho_subj,0.8);

  alpha_pop ~ normal(0.5, 0.3);

  alpha_subj ~ normal(0.5, 0.3);

  mu ~ normal(2, 1);
  
  lambda_shape ~ normal(0.5,0.2);
  //log_k_minus1 ~ normal(0, 1);
  //k ~ exponential(lambda_shape) T[1, ];
  k ~ exponential(lambda_shape) T[1, ];
  array[3] vector[N_total] eta_quad;

  for (j in 1:3) {

    vector[N_total] f_pop = PHI_quad[j] * spectral_pop;

    vector[N_total] f_subj = rows_dot_product(PHI_quad[j], beta_subj[s_id, ]);

    // global centering:

    //f_pop -= mean(f_pop);

    eta_quad[j] = mu[s_id] + f_pop + f_subj;
  }


  for (n in 1:N_total) {

    vector[3] log_kernel;

    for (j in 1:3) {

      real theta = exp(-eta_quad[j][n]);

      log_kernel[j] = log_qw[j] + gamma_lpdf(dt[n] | k, theta);
    }

    target += log_sum_exp(log_kernel);
  }
}
