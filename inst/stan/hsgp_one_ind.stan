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
  int<lower=0, upper=N_params> include_mu_lognormal;
  int<lower=0, upper=N_params> include_sigma_lognormal;
  // data
  int<lower=1> N_total;
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
  vector[M] z;
  real<lower=0> rho;
  real<lower=0> alpha;

  // subject intercept
  real mu;

  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector[include_mu_lognormal ? 1 : 0] mu_lognormal;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {


  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S;

  if (kernel == 1) {
    diag_S = diagSPD_EQ(alpha,rho,M,L);
  } else if (kernel == 2) {
    diag_S = diagSPD_Matern12(alpha,rho,M,L);
  } else if (kernel == 3) {
    diag_S = diagSPD_Matern32(alpha,rho,M,L);
  } else if (kernel == 4) {
    diag_S = diagSPD_Matern52(alpha,rho,M,L);
  } else if (kernel == 5) {
    diag_S = diagSPD_Periodic(alpha,rho,M);
  }

  vector[M] beta = diag_S .* z;

}

model {

  // Priors
  z ~ std_normal();

  apply_prior_lp(mu, distributions[1], params[1, 1], params[1, 2]);
  apply_prior_lp(alpha, distributions[2], params[2, 1], params[2, 2]);
  apply_prior_lp(rho, distributions[3], params[3, 1], params[3, 2]);

  if (include_k != 0) apply_prior_lp(k[1], distributions[include_k],
    params[include_k, 1], params[include_k, 2]);
  if (include_shape != 0) apply_prior_lp(shape[1], distributions[include_shape],
    params[include_shape, 1], params[include_shape, 2]);
  if (include_mu_lognormal != 0) apply_prior_lp(mu_lognormal[1], distributions[include_mu_lognormal],
    params[include_mu_lognormal, 1], params[include_mu_lognormal, 2]);
  if (include_sigma_lognormal != 0) apply_prior_lp(sigma_lognormal[1], distributions[include_sigma_lognormal],
    params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2]);

  array[3] vector[N_total] eta_quad;

  for (j in 1:3) {

    vector[N_total] f = PHI_quad[j] * beta;

    eta_quad[j] = mu + f;
  }


  for (n in 1:N_total) {

    vector[3] log_kernel;


    if (family == 1) {
      for (j in 1:3) {
        real scale = exp(eta_quad[j][n]);
        log_kernel[j] = log_qw[j] + exponential_lpdf(dt[n] | inv(scale));
      }
    } else if (family == 2) {
      for (j in 1:3) {
        real scale = exp(eta_quad[j][n]);
        log_kernel[j] = log_qw[j] + gamma_lpdf(dt[n] | k[1], inv(scale));
      }
    } else if (family == 3) {
      for (j in 1:3) {
        real scale = exp(eta_quad[j][n]);
        log_kernel[j] = log_qw[j] + weibull_lpdf(dt[n] | shape[1], scale);
      }

    } else if (family == 4) {
      for (j in 1:3) {
        real scale = exp(eta_quad[j][n]);
        log_kernel[j] = log_qw[j] + lognormal_lpdf(dt[n] | mu_lognormal, sigma_lognormal[1]);
      }
    } else if (family == 5) {
      for (j in 1:3) {
        real scale = exp(eta_quad[j][n]);
        log_kernel[j] = log_qw[j] + gengamma_lpdf(dt[n] | shape[1], k[1], scale);
      }
    }

    target += log_sum_exp(log_kernel);
  }
}
