#' Title
#'
#' @param X input x vector
#' @param Y input binary response variable
#' @param prior choose the shrinkage prior to use, can be "Regularized HS", "Orginal HS", "Laplacian" and "Gaussian"
#' @param L input the number of user-defined interval
#' @param tau0_sq Input hyper-parameter for tau
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

monotoneBayes = function(X, Y, Z=NULL, L = 10, tau0_sq = 1e-2, nodes = seq(0,1,length.out = 10+1), Eq.Space = T,
                         c_sq = 10^2, fix = F,
                         c_alpha = 1, c_beta = 1 * 200, prior = "Regularized HS", ...){
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
    data.J[data.J>=L] = L
    data.W = rep(1/L,L)
  }

  if (length(data.W)==1) data.W = array(data.W, dim = 1)
  if (is.null(Z)){
    dt.stan = list(Y=Y, X = X, J = data.J,  W = data.W, nodes = nodes, L=L, N=N,
                   local_dof_stan = 1,
                   global_dof_stan = 1,
                   tau0_sq = tau0_sq)
    if (prior == "Original HS"){
      model = rstan::sampling(stanmodels$OrgHS, data = dt.stan, ...)
    } else if (prior == "Laplacian"){
      model = rstan::sampling(stanmodels$Laplacian, data = dt.stan, ...)
    } else if (prior == "Gaussian"){
      model = rstan::sampling(stanmodels$Gaussian, data = dt.stan, ...)
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
  } else {
    dt.stan = list(Y=Y, X = X, Z=Z, J = data.J,  W = data.W, nodes = nodes, L=L, N=N,
                   local_dof_stan = 1,
                   global_dof_stan = 1, c_sq_shape = c_alpha, c_sq_scale = c_beta,
                   tau0_sq = tau0_sq)

    model = rstan::sampling(stanmodels$CovModel_HS, data = dt.stan, ...)
  }



  return(model)
}
