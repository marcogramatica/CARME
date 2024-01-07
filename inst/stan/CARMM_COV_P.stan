functions {
  real sparse_car_lpdf(vector phi,
  real alpha,
    array[,] int W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (sum(ldet_terms) - (phit_D * phi - alpha * (phit_W * phi)));
  }
  real corr_val(vector x, vector y, int n){
    int n_vec;
    real num;
    real den_1;
    real den_2;
    n_vec = n;
    num = n_vec*sum(x .* y) - sum(x)*sum(y);
    den_1 = sqrt(n_vec*sum(x .* x) - pow(sum(x), 2));
    den_2 = sqrt(n_vec*sum(y .* y) - pow(sum(y), 2));
    return num/(den_1 * den_2);
    }
}
data {
  // Number of areas
  int<lower=0> n;

  // Multiple membership
  int<lower=0> m;
  matrix[m, n] H;

  // Spatial matrix
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n; // number of adjacent region pairs

  // Covariates
  int<lower=0> k;   // number of predictors
  matrix[n, k] X_cov;   // predictor matrix

  // Outcomes
  array[m] int<lower = 0> y;

  // Offsets
  vector[m] log_offset;

  // Hyperparameters
  real mu_gamma;
  real mu_beta;
  real<lower = 0> sigma_gamma;
  real<lower = 0> sigma_beta;
  real<lower = 0> tau_shape;
  real<lower = 0> tau_rate;
}
transformed data {
   // Adjacency
  // Number of directed neighbours per area
  vector[n] n_i;
  array[W_n, 2] int W_sparse;   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:n) D_sparse[i] = sum(W[i]);
  {
    vector[n] invsqrtD;
    for (i in 1:n) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}
parameters {
  real<lower = 0, upper = 1> alpha;
  real<lower = 0> tau;
  vector[n] phi_unsc;

  // Intercept and covariates
  real gamma;
  vector[k] beta;
}
transformed parameters{
  vector[m] r_mm;
  real<lower = 0> invtausq;
  vector[n] phi;

  invtausq = inv_sqrt(tau);
  phi = invtausq*phi_unsc;

  // Relative risk
  r_mm = H*(gamma + X_cov * beta + phi);
}
model {
  phi_unsc ~ sparse_car(alpha, W_sparse, D_sparse, lambda, n, W_n);
  tau ~ gamma(tau_shape, tau_rate);
  beta ~ normal(mu_beta, sigma_beta);
  gamma ~ normal(mu_gamma, sigma_gamma);

  // Poisson
  y ~ poisson(exp(log_offset + r_mm));
}
generated quantities {
  array[m] int<lower = 0> yrep;
  vector[m] log_lik;
  vector[m] log_lik_rep;
  real sum_ll;
  real sum_ll_rep;
  real ppp;
  vector[n] l_RR;

  // Areal relative risk
   l_RR = gamma + X_cov * beta + phi;

  // Simulate from posterior
  for (i in 1:m) {
  //   // likelihood of the current parameter values (given the original data)
    log_lik[i] = poisson_lpmf(y[i] | exp(r_mm[i] + log_offset[i]));
  //   // generate new data based on current parameter values
    yrep[i] = poisson_rng(exp(r_mm[i] + log_offset[i]));
  //   // compute likelihood of the current parameter values (given the new data)
    log_lik_rep[i] = poisson_lpmf(yrep[i] | exp(r_mm[i] + log_offset[i]));
  }
  //   // sum up the likelihoods for all observations
    sum_ll = sum(log_lik);
    sum_ll_rep = sum(log_lik_rep);
  //   // check which is higher
    ppp = sum_ll > sum_ll_rep ? 1 : 0;
}
