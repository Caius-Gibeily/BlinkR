// spectral basis set, phi
matrix phi(int N, int M, real L, vector x) {
  matrix[N, M] res;

  for (m in 1:M) {
    res[, m] = inv_sqrt(L) * sin(pi() * m * (x + L) / (2 * L));
  }

  return res;
}

matrix phi_periodic(int N, int M, real w0, vector x) {

  row_vector[M / 2] k = linspaced_row_vector(M / 2, 1, M / 2);

  matrix[N, M / 2] w0xk = (w0 * x) * k;
  return append_col(cos(w0xk), sin(w0xk));
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
  vector[M / 2] indices = linspaced_vector(M / 2, 1, M / 2);
  vector[M / 2] q = exp(
    log(alpha) + 0.5 *
    (log(2) - a + to_vector(log_modified_bessel_first_kind(indices, a)))
    );
    return append_row(q,q);
}

vector diagSPD_EQ(real alpha, real rho, int M, real L) {
  vector[M] indices = linspaced_vector(M, 1, M);

  real factor = alpha^2 * sqrt(2 * pi()) * rho;
  real exponent = -0.5 * (rho * pi() / (2 * L))^2;

  return factor * exp(exponent * square(indices));
}

real gengamma_lpdf(real x, real shape, real k, real scale) {

  if (x <= 0.0) {
    return negative_infinity();
  }
  real logdens = log(shape) - lgamma(k) + (shape * k - 1.0) * log(x) - (shape * k) * log(scale) - pow(x/scale,shape);
  return logdens;

}

vector get_diagSPD(real alpha, real rho, int M, real L, int kernel) {

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
  return diag_S;
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

matrix exponential_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + exponential_lpdf(dt[n] | inv(scale[n]));
    }
  }
  return log_kernel;
}

matrix gamma_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real k) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + gamma_lpdf(dt[n] | k, inv(scale[n]));
    }
  }
  return log_kernel;
}

matrix weibull_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real shape) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + weibull_lpdf(dt[n] | shape, scale[n]);
    }
  }
  return log_kernel;
}

matrix lognormal_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real sigma_lognormal) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] mu_lognormal = eta_quad[j];
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + lognormal_lpdf(dt[n] | mu_lognormal[n], sigma_lognormal);
    }
  }
  return log_kernel;
}

matrix gengamma_likelihood(int N_total, vector log_qw, array[] vector eta_quad, vector dt, real k, real shape) {

  matrix[N_total, 3] log_kernel;

  for (j in 1:3) {
    vector[N_total] scale = exp(eta_quad[j]);
    for (n in 1:N_total) {
      log_kernel[n, j] = log_qw[j] + gengamma_lpdf(dt[n] | shape, k, scale[n]);
    }
  }
  return log_kernel;
}
