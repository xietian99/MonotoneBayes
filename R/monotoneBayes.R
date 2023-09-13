#' Title
#'
#' @param X input x vector
#' @param Y input binary response variable
#' @param prior choose the shrinkage prior to use, can be "Regularized HS", "Orginal HS", "Laplacian" and "Gaussian"
#' @param L input the number of user-defined interval
#' @param tau0 input hyperparameter
#' @param c_alpha input hyperparameter for c
#' @param c_beta input hyperparameter for c
#'
#' @return a list
#' $model: the model result
#' $plots: the fitted curve with original dataset
#' @export
#'
#' @examples
#' X <- runif(10, 0, 1)
#' Y <- rbinom(10, 1, 0.2)
#' monotoneBayes(X, Y)
#'
#' @importFrom rstan sampling stan

monotoneBayes = function(X, Y, L = 10, tau0 = 1e-4,
                         c_sq = 10^2, fix = F,
                         c_alpha = 0.01, c_beta = 0.01 * 4, prior = "Regularized HS"){
  N = length(Y)
  data.L = X %/% (1/L) + 1
  dt.stan = list(Y=Y, X = X, J = data.L,  L=L, N=N,
                  local_dof_stan = 1,
                  global_dof_stan = 1,
                  tau0_sq = tau0^2)
  if (prior == "Original HS"){
    model = rstan::sampling(stanmodels$OrgHS, data = dt.stan)
  } else if (prior == "Laplacian"){
    model = "Laplacian ... TBD..."
  } else if (prior == "Gaussian ... TBD..."){
    model = "Gaussian"
  } else {
    if (fix == T) {
      dt.stan$c_sq = c_sq
      model = rstan::sampling(stanmodels$RegHSfix, data = dt.stan)
    } else {
      dt.stan$c_sq_shape = c_alpha
      dt.stan$c_sq_scale = c_beta
      model = rstan::sampling(stanmodels$RegHS, data = dt.stan)
    }
  }


  return(model)
}
