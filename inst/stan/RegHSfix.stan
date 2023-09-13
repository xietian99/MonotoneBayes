data {
  // Define N, L, Y[N], X[N]
  int<lower=0> N;
  int<lower=1> L;
  int<lower = 0, upper = 1> Y[N]; // outcome binary variable
  vector<lower=0.0, upper=1.0>[N] X; // dose level
  int<lower = 1, upper = L> J[N]; // which interval X falls in

  // C_sq for fixed c_sq
  real<lower = 0.0> c_sq;
  real<lower = 0.0> local_dof_stan;
  real<lower = 0.0> global_dof_stan;
  real<lower = 0.0> tau0_sq;
}

// The parameters accepted by the model. Our model
// accepts parameters 'alpha's and 'tau's, and gamma
parameters {
  // Main Parameter Estimated
  vector<lower = 0.0>[L] alpha_base;

  // tau: half-cauchy , stand cauchy * scale ;
  // decomposed into chi-square (N^2) and inverse-gamma
  real<lower = 0.0> tau_base_sq;
  real<lower = 0.0> tau_scale_sq;
  // lambda: half-cauchy decomposed into chi-square (N^2) and inverse-gamma
  vector<lower = 0.0>[L] lambda_base_sq;
  vector<lower = 0.0>[L] lambda_scale_sq;

  // Beta0 / Beta1: decomposed into chi-square (N^2) and inverse-gamma
  real<lower = 0.0> gamma;
}

transformed parameters {
  // Interested: alphas and gammas
  vector<lower = 0.0>[L] alpha;


  // Used parameter
  vector<lower = 0.0>[L] theta; // alpha = alpha_base * theta
  real<lower = 0.0> tau_sq;
  vector<lower = 0.0>[L] lambda_sq;
  real<lower = 0.0> denominator;

  //phil Xi: fitted probabilities at x=0,1/L,...,1
  vector<lower = 0.0, upper = 1.0>[L+1] xi;

  // Define Bernoulli B(p)
  vector<lower = 0.0, upper = 1.0>[N] p;

  // parameter transformation
  tau_sq = tau_base_sq * tau_scale_sq;
  lambda_sq = lambda_base_sq .* lambda_scale_sq;

  for (i in 1:L){
    // Below is where scaling by tau0 comes into play.
    theta[i] = 1 / sqrt( 1/(tau0_sq * tau_sq * lambda_sq[i]) + 1/c_sq );
  }
  alpha = alpha_base .* theta;
  denominator = sum(alpha)/L + 1 + gamma;
  for (i in 1:N){
    if (J[i] == 1)
      p[i] = ( X[i] * alpha[1] + 1 ) / denominator  ;
    else
      p[i] = ( sum(alpha[1:(J[i]-1)])/L + (X[i]-(J[i]-1) *1.0/L) * alpha[J[i]] + 1  ) / denominator;
  }
  xi[1] = 1 / denominator;
  for(i in 2:(L+1)) {
    xi[i] = (1 + sum(alpha[1:(i-1)]) / L) / denominator;
  }

}

// The model to be estimated.
model {
  // Priors
  alpha_base ~ normal(0.0, 1.0);
  tau_base_sq ~ chi_square(1.0);
  tau_scale_sq ~ inv_gamma(global_dof_stan/2.0, global_dof_stan/2.0);
  lambda_base_sq ~ chi_square(1.0);
  lambda_scale_sq ~ inv_gamma(local_dof_stan/2.0, local_dof_stan/2.0);
  // likelihood
  for (i in 1:N) {
    Y[i] ~ bernoulli(p[i]);
  }
}

