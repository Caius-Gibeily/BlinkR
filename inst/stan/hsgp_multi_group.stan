functions {
#include include/functions.stan
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
  real<lower=0> w0;

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

  matrix[G,G-1] Q_R = qr_decomp(G);

}

parameters {

  // population GP
  vector[M] z_global;
  matrix[G-1,M] z_group_raw;
  matrix[I,M] z_ind_raw;

  real<lower=0> rho_global;
  real<lower=0> alpha_global;

  real<lower=0> rho_group;
  real<lower=0> alpha_group;

  // Individual-level
  real<lower=0> rho_ind;
  real<lower=0> alpha_ind;


  // indect intercept
  //vector[G] mu_group;
  //vector[I] mu_ind;
  real mu_global;
  real<lower=0> sigma_group;
  real<lower=0> sigma_ind;

  vector[G - 1] mu_raw_group;
  vector[I] mu_raw_ind;
  // renewal shape
  //real log_k_minus1;
  //vector<lower=1>[S] k;

  vector<lower=1>[include_k ? 1 : 0] k;
  vector<lower=1>[include_shape ? 1 : 0] shape;
  vector<lower=0>[include_sigma_lognormal ? 1 : 0] sigma_lognormal;


}

transformed parameters {

  // zero-constrained group mu intercepts
  vector[G] mu_raw_group_std = Q_R * mu_raw_group; //sum_to_zero_mu(G, 1, rep_array(1,G), mu_raw_group, {G});
  vector[G] mu_group = mu_global + sigma_group * mu_raw_group_std;

  // zero-constrained individual mu intercepts
  vector[I] mu_raw_ind_std = constrain_mu(I, G, g_membership, mu_raw_ind, I_per_group);
  vector[I] mu_ind = mu_group[g_membership] + sigma_ind * mu_raw_ind_std;


  matrix[I,M] z_ind;
  matrix[G,M] z_group;

  //real k = 1.0 + exp(log_k_minus1);
  vector[M] diag_S_global;
  vector[M] diag_S_group;
  vector[M] diag_S_ind;

  diag_S_global = get_diagSPD(alpha_global,rho_global,M,L,kernel);
  diag_S_group = get_diagSPD(alpha_group,rho_group,M,L,kernel);
  diag_S_ind = get_diagSPD(alpha_ind,rho_ind,M,L,kernel);


  //z_group[1:(G - 1), ] = z_group_raw;

  //for (m in 1:M) {
  //  z_group[G, m] = -sum(z_group_raw[, m]);
  //}
  z_group = Q_R * z_group_raw;

  z_ind = constrain_groups(G, I, M, g_membership, z_ind_raw, I_per_group);

  vector[M] beta_global = diag_S_global .* z_global;

  matrix[G, M] beta_group = diag_post_multiply(z_group,diag_S_group);

  matrix[I, M] beta_ind = diag_post_multiply(z_ind,diag_S_ind);

}

model {

  // Priors
  z_global ~ std_normal();
  to_vector(z_group_raw) ~ std_normal();
  to_vector(z_ind_raw) ~ std_normal();

  mu_raw_group ~ std_normal();
  mu_raw_ind ~ std_normal();

  sigma_group ~ normal(0.5,0.2);
  sigma_ind ~ normal(0.5,0.2);

  apply_prior_lp(mu_global, distributions[1], params[1, 1], params[1, 2]);

  apply_prior_lp(alpha_global, distributions[2], params[2, 1], params[2, 2]);
  apply_prior_lp(alpha_group, distributions[3], params[3, 1], params[3, 2]);
  apply_prior_lp(alpha_ind, distributions[4], params[4, 1], params[4, 2]);

  apply_prior_lp(rho_global, distributions[5], params[5, 1], params[5, 2]);
  apply_prior_lp(rho_group, distributions[6], params[6, 1], params[6, 2]);
  apply_prior_lp(rho_ind, distributions[7], params[7, 1], params[7, 2]);


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


    eta_quad[j] = mu_ind[ind_id] + f_global + f_group + f_ind;
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
