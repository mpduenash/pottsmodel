functions {
  real sparse_car_lpdf(vector phi, real tau, real alpha,
                       array[,] int W_sparse, vector D_sparse, vector lambda,
                       int n, int W_n) {
    row_vector[n] phit_D; // phi' * D
    row_vector[n] phit_W; // phi' * W
    vector[n] ldet_terms;
    
    phit_D = (phi .* D_sparse)';
    phit_W = rep_row_vector(0, n);
    for (i in 1 : W_n) {
      phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
      phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
    }
    
    for (i in 1 : n) {
      ldet_terms[i] = log1m(alpha * lambda[i]);
    }
    return 0.5
           * (n * log(tau) + sum(ldet_terms)
              - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}

data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> k;
  int<lower=1> n_miss; //Number of missing
  int<lower=1> n_pres; //Number of present values
  matrix[n, p] X;
  array[n] int<lower=0> y; //Responses included missing responses
  matrix<lower=0, upper=1>[n, n] W; // adjacency matrix
  int W_n; // number of adjacent region pairs
  array[n_miss] int<lower=0, upper=n> test_block;  
  array[n_pres] int<lower=0, upper=n> train_block;  
}

transformed data {
  array[W_n, 2] int W_sparse; // adjacency pairs
  vector[n] D_sparse; // diagonal of D (number of neigbors for each site)
  vector[n] lambda; // eigenvalues of invsqrtD * W * invsqrtD
  
  {
    // generate sparse representation for W
    int counter;
    counter = 1;
    // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1 : (n - 1)) {
      for (j in (i + 1) : n) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1 : n) {
    D_sparse[i] = sum(W[i]);
  }
  {
    vector[n] invsqrtD;
    for (i in 1 : n) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
  
  matrix[n_miss, n_miss] W_MM;
  matrix[n_miss, n_pres] W_MO;
  vector[n_miss] D_M;

  for (i in 1:n_miss) {

  int mi = test_block[i];

  D_M[i] = D_sparse[mi];

  for (j in 1:n_miss) {
    int mj = test_block[j];
    W_MM[i,j] = W[mi,mj];
  }

  for (j in 1:n_pres) {
    int oj = train_block[j];
    W_MO[i,j] = W[mi,oj];
  }
}
  
}
parameters {
  matrix[p, k-1] raw_beta; //beta in matrix form
  matrix[n, k-1] raw_phi; //phi in matrix form
  array[k-1] real<lower=0> tau; 
  array[k-1] real<lower=0, upper=0.999> alpha;
}

transformed parameters{
  
  vector[p] zeros = rep_vector(0, p); //fixed effects for reference color
  vector[n] zeros_phi = rep_vector(0, n); //"random" effects for reference color
  matrix[p, k] beta;
  matrix[n, k] phi;
  beta = append_col(zeros, raw_beta); 
  phi = append_col(zeros_phi, raw_phi);
  
  matrix[n, k] X_beta = X * beta;
  
}


model {
  // Prior for beta
  to_vector(raw_beta) ~ normal(0, 1);
  
  for(j in 1:(k-1)){
    tau[j] ~ gamma(2, 2); 
    col(raw_phi, j) ~ sparse_car(tau[j], alpha[j], W_sparse, D_sparse, lambda, n, W_n); 
  }
  
  // Likelihood
  
  for (i in 1:n_pres){
    int index = train_block[i];
    vector[k] eta = row(X_beta, index)' + row(phi, index)';
    y[index] ~ categorical_logit(eta);
  }
}

generated quantities {
  
    matrix[n_miss, k] probs;
    matrix[n_miss, k] phi_miss;
    vector[n_miss] y_imp;
  
  for (i in 1:n_miss) {
    phi_miss[i, 1] = 0;
  }

  for (j in 1:(k-1)) {

    matrix[n_miss, n_miss] A;
    matrix[n_miss, n_miss] L;

    vector[n_pres] phi_obs;
    vector[n_miss] mu_m;
    vector[n_miss] z;

    
    A =
      diag_matrix(D_M)
      - alpha[j] * W_MM;

    L = cholesky_decompose(A);


    for (i in 1:n_pres) {

      int oi = train_block[i];

      phi_obs[i] = raw_phi[oi, j];
    }

    mu_m =
      mdivide_left_spd(
        A,
        alpha[j] * W_MO * phi_obs
      );


    for (i in 1:n_miss) {
      z[i] = normal_rng(0, 1);
    }

   phi_miss[, j + 1] =
    mu_m
     + mdivide_right_tri_low(z', L)' / sqrt(tau[j]);
  }


  /*
   * ------------------------------------------------
   * Predict probabilities
   * ------------------------------------------------
   */

  
  for (i in 1:n_miss) {

  int index = test_block[i];

  vector[k] eta;

  eta =
    to_vector(row(X, index) * beta)
    + to_vector(phi_miss[i]);

  probs[i] =
    softmax(eta)';
    
    y_imp[i] = categorical_logit_rng(eta);
}
}
  


