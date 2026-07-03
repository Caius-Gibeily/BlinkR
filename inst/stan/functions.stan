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

  // Matern 1/2
  vector diagSPD_Matern12(real alpha, real rho, int M, real L) {

    vector[M] denom =
      rho * (square(inv(rho)) + square((pi() / (2 * L)) * basis_indices(M, L)));

    return alpha * sqrt(2.0 ./ denom));
  }

  // Matern 3/2
  vector diagSPD_Matern32(real alpha, real rho, int M, real L) {

    real factor = 2 * alpha * (sqrt(3) / rho)^1.5;
    vector[M] denom = 3 / square(rho) + basis_indices(M, L);

    return factor ./ denom;
  }

  // Matern squared exponential
  vector diagSPD_EQ(real alpha, real rho, int M, real L) {

    real factor = alpha * sqrt(sqrt(2 * pi()) * rho);
    real exponent = -0.25 * (rho * pi() / 2 / L)^2;

    return factor * exp(exponent * square(basis_indices(M, L)));
  }

}
