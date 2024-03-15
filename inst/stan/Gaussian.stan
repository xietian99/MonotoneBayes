// Normal Prior
// The input data is a matrix of rows 'N', consisting of 'X', 'Y', and fixed "L" used for model definition
data {
  // Define N, L, Y[N], X[N], M is estimated non-zero alphas.
  int<lower=0> N;
  int<lower=1> L;
  vector<lower=0.0, upper=1.0>[L+1] nodes; // user-defined node point
  vector<lower=0.0, upper=1.0>[L] W; // user-defined length of each interval
  int<lower = 0, upper = 1> Y[N]; // outcome binary variable
  vector<lower=0.0, upper=1.0>[N] X; // dose level
  int<lower = 1, upper = L> J[N]; // which interval X falls in


  // Define the global hyperparameter of priors as to tau
  real<lower = 0.0> local_dof_stan;
  real<lower = 0.0> global_dof_stan;
  real<lower = 0.0> tau0_sq;

}

// The parameters accepted by the model. Our model
// accepts parameters 'alpha's and 'beta's.
parameters {
  // Main Parameter Estimated
  vector<lower = 0.0>[L] alpha_base;


  // Hyper-parameter
  // tau: half-cauchy , stand cauchy * scale ;
  // decomposed into chi-square (N^2) and inverse-gamma
  real<lower = 0.0> tau_base_sq;
  real<lower = 0.0> tau_scale_sq;

  real<lower = 0.0> gamma;

}

transformed parameters {
  // Interested: alphas and betas
  vector<lower = 0.0>[L] alpha;
  vector<lower = 0.0>[L] cumsum_alphaW;

  // Used parameter
  vector<lower = 0.0>[L] theta; // alpha = alpha_base * theta
  real<lower = 0.0> tau_sq;
  real<lower = 0.0> denominator;

  //phil Xi: fitted probabilities at x=0,1/L,...,1
  vector<lower = 0.0, upper = 1.0>[L+1] xi;


  // Define Bernoulli B(p)
  vector<lower = 0.0, upper = 1.0>[N] p;



  //phil note that tau_sq is *unscaled* by tau0.
  tau_sq = tau_base_sq * tau_scale_sq;

  // lambda_sq = lambda_base_sq .* lambda_scale_sq;
  for (i in 1:L){
    theta[i] = sqrt(tau0_sq * tau_sq);
  }
  alpha = alpha_base .* theta;
  denominator = sum(alpha .* W) + 1 + gamma;
  cumsum_alphaW = cumulative_sum(alpha .* W);
  for (i in 1:N){
    if (J[i] == 1)
       p[i] = ( X[i] * alpha[1] + 1 ) / denominator  ;
    else
      p[i] = ( cumsum_alphaW[J[i]-1] + (X[i]-nodes[J[i]]) * alpha[J[i]] + 1  ) / denominator;
  }
  xi[1] = 1 / denominator;
  for(i in 2:(L+1)) {
    xi[i] = ( cumsum_alphaW[i-1] + 1) / denominator;
  }

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // Priors
  alpha_base ~ normal(0.0, 1.0);
  tau_base_sq ~ chi_square(1.0);
  tau_scale_sq ~ inv_gamma(global_dof_stan/2.0, global_dof_stan/2.0);
  // c_sq ~ inv_gamma(c_sq_shape, c_sq_scale);
  // likelihood
  for (i in 1:N) {
    Y[i] ~ bernoulli(p[i]);
  }
}


