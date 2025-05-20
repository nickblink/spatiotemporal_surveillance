functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi_star CAR prior phi minus the expected value of phi (alpha times the previous iteration)
  * @param tau2 variance parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param log_det_Q The precomputed log determinant of the precision matrix
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi_star, real tau2, real rho, int[,] W_sparse, 
  vector D_sparse, real log_det_Q, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
    
      phit_D = (phi_star .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi_star[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi_star[W_sparse[i, 1]];
      }
    
      return 0.5 * (log_det_Q - (1/tau2) * (rho * (phit_D * phi_star) - rho * (phit_W * phi_star) + (1 - rho)*dot_self(phi_star)));
  }
}data {
  int<lower=0> p; // number of variables
  int<lower=0> N;  // number of observations
  int<lower=0> N_F;  // number of facilities
  int<lower=0> N_T;  // number of time points
  int<lower=0> N_miss; // number of missing y points
  int<lower=0> N_obs; // number of observed y points
  int<lower=0, upper=N> ind_miss[N_miss]; // indices of missing y points
  int<lower=0, upper=N> ind_obs[N_obs]; // indices of observed y points
  matrix[N,p] X; // design matrix of observed values
  int<lower=0> y_obs[N_obs];  // output
  vector[p] mu; // prior mean betas
  matrix[p,p] Sigma; // prior variance betas
  int<lower=0> theta_shape; // theta prior shape
  int<lower=0> theta_rate; // theta prior rate
  //matrix[N_F,N_F] W_star; // D - W matrix for Leroux CAR covariance
  matrix<lower=0, upper = 1>[N_F, N_F] W; //adjacency matrix
  int W_n; // Number of adjacency pairs
  matrix[N_F,N_F] I; // Identity matrix
  vector[N_F] lambda; // the eigenvalues of the D - W - I matrix
}
transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N_F] D_sparse;     // diagonal of D (number of neigbors for each site)
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N_F - 1)) {
      for (j in (i + 1):N_F) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N_F) D_sparse[i] = sum(W[i]); // Compute the sparse representation of D
}
parameters {
  // int<lower=0> y_miss[N_miss];
  real<lower=0> tau2; // CAR variance parameter
  real<lower=0, upper=1> rho; // spatial correlation
  real<lower=0, upper=1> alpha; // temporal correlation
  vector<lower=0>[N_F] theta; // theta vector
  vector[N] phi; // CAR parameter in vector form
  vector[p] beta;  // log of rate parameter
}
transformed parameters {
  // variable declarations
  vector[N_obs] observed_est;
  vector[N] theta_long;
  vector[N_obs] theta_observed;
  vector[N] phi_star;
  vector[N_F + 1] ldet_vec;
  real log_detQ;
  
  // variable calculations
  observed_est = (X*beta)[ind_obs] + phi[ind_obs];
  for(t in 1:N_T){
    theta_long[((t-1)*N_F+1):(t*N_F)] = theta;
  }
  theta_observed = theta_long[ind_obs];
  phi_star = phi;
  // center the phi_star values from a random walk with temporal correlation alpha
  for (t in 2:N_T){
    phi_star[((t-1)*N_F+1):(t*N_F)] = phi[((t-1)*N_F+1):(t*N_F)] - alpha*phi[((t-2)*N_F+1):((t-1)*N_F)];
  }
  ldet_vec[N_F + 1] = -N_F*log(tau2);
  for (i in 1:N_F){
	ldet_vec[i] = log1p(rho*lambda[i]);
  }
  log_detQ = sum(ldet_vec);
}
model {
  beta ~ multi_normal(mu, Sigma); // beta prior
  theta ~ gamma(theta_shape, theta_rate); // theta prior
  y_obs ~ neg_binomial_2_log(observed_est, theta_observed); // negbin likelihood.
  phi_star[1:N_F] ~ sparse_car(tau2, rho, W_sparse, D_sparse, log_detQ, N_F, W_n); // first time point CAR prior
  // CAR prior for successive time points
  for (t in 2:N_T) {
    phi_star[((t-1)*N_F+1):(t*N_F)] ~ sparse_car(tau2, rho, W_sparse, D_sparse, log_detQ, N_F, W_n);
  }
}
generated quantities {
  vector[N] y_exp = exp(X*beta + phi);
  int y_pred[N] = neg_binomial_2_log_rng(X*beta + phi, theta_long);
}
