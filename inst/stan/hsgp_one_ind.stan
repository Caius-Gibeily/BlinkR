functions {
#include include/functions.stan
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
  real<lower=0> w0;

  int<lower=0, upper=N_params> include_k;
  int<lower=0, upper=N_params> include_shape;
  int<lower=0, upper=N_params> include_sigma_lognormal;
  // data
  int<lower=1> N_total;
  real<lower=0> duration;
  vector[N_total] t_ev;
  vector[N_total] dt;


}

transformed data {

  real L = L_factor * duration;

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

  if (kernel != 5) {
    for (j in 1:3) {

      vector[N_total] t_quad = t_ev - (1.0 - qx[j]) .* dt;

      PHI_quad[j] = phi(N_total,M,L,t_quad);
    }
  } else {
    for (j in 1:3) {

      vector[N_total] t_quad = t_ev - (1.0 - qx[j]) .* dt;

      PHI_quad[j] = phi_periodic(N_total,M,w0,t_quad);
    }
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
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {


  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S;

  diag_S = get_diagSPD(alpha,rho,M,L,kernel);

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
  if (include_sigma_lognormal != 0) apply_prior_lp(sigma_lognormal[1], distributions[include_sigma_lognormal],
    params[include_sigma_lognormal, 1], params[include_sigma_lognormal, 2]);

  array[3] vector[N_total] eta_quad;

  for (j in 1:3) {

    vector[N_total] f = PHI_quad[j] * beta;

    eta_quad[j] = mu + f;
  }

  matrix[N_total, 3] log_kernel;

  if (family == 1) log_kernel = exponential_likelihood(N_total, log_qw, eta_quad, dt);
  else if (family == 2) log_kernel = gamma_likelihood(N_total, log_qw, eta_quad, dt, k[1]);
  else if (family == 3) log_kernel = weibull_likelihood(N_total, log_qw, eta_quad, dt, shape[1]);
  else if (family == 4) log_kernel = lognormal_likelihood(N_total, log_qw, eta_quad, dt, sigma_lognormal[1]);
  else if (family == 5) log_kernel = gengamma_likelihood(N_total, log_qw, eta_quad, dt, k[1], shape[1]);

  for (n in 1:N_total) {
    target += log_sum_exp(log_kernel[n, ]);
  }


}
