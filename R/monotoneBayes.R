#' Title
#'
#' @param X input x vector
#' @param Y input binary response variable
#' @param prior choose the shrinkage prior to use, can be "Regularized HS", "Orginal HS", "Laplacian" and "Gaussian"
#' @param L input the number of user-defined interval
#' @param tau0 input hyperparameter
#' @param c_alpha input hyperparameter for c
#' @param c_beta input hyperparameter for c
#' @param c_sq fixed C_sq value for regularized HS with fixed C.
#' @param fix default F
#' @param ... chain parameter passing to stan, inclcuding iter, warmup, control ,...
#' @param nodes default 0,1/L,2/L,...L/L: user-defined nodes point
#' @param Eq.Space default T, should be F when user defined interval not equal space
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
#' @importFrom stats quantile

monotoneBayes = function(X, Y, L = 10, tau0 = 1e-2, nodes = seq(0,1,length.out = 10+1), Eq.Space = T,
                         c_sq = 10^2, fix = F,
                         c_alpha = 1, c_beta = 1 * 4, prior = "Regularized HS", ...){
  N = length(Y)
  if (Eq.Space == F) {
    #nodes = (nodes - min(X))/(max(X) - min(X))
    #X = (X - min(X))/(max(X) - min(X))
    data.J = rowSums(outer(X, nodes, "-")>=0)
    data.J[data.J>=L] = L
    data.W = nodes[2:(L+1)]- nodes[1:L]
  } else {
    nodes = seq(0, 1, length.out = L+1)
    #X = (X - min(X))/(max(X) - min(X))
    data.J = X %/% (1/L) + 1
    data.W = 1/L
  }


  dt.stan = list(Y=Y, X = X, J = data.J,  W = data.W, nodes = nodes, L=L, N=N,
                  local_dof_stan = 1,
                  global_dof_stan = 1,
                  tau0_sq = tau0^2)
  if (prior == "Original HS"){
    model = rstan::sampling(stanmodels$OrgHS, data = dt.stan, ...)
  } else if (prior == "Laplacian"){
    model = "Laplacian ... TBD..."
  } else if (prior == "Gaussian ... TBD..."){
    model = "Gaussian"
  } else {
    if (fix == T) {
      dt.stan$c_sq = c_sq
      model = rstan::sampling(stanmodels$RegHSfix, data = dt.stan, ...)
    } else {
      dt.stan$c_sq_shape = c_alpha
      dt.stan$c_sq_scale = c_beta
      model = rstan::sampling(stanmodels$RegHS, data = dt.stan, ...)
    }
  }


  return(model)
}
