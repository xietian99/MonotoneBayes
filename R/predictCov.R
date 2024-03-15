#' Currently this predict function works only for cases with new coavriates.
#'
#' @param model fitted stan model
#' @param method mean or median
#' @param grid.eq equal grid, if F need to provide node.t below
#' @param node.t default equal nodes on 0-1
#' @param newX input new dose, should be transformed the same way with previous dose, could be extrapolation
#' @param newZ input new covariates
#' @param returnCI if T will return 95% Crediable Interval, otherwise return posterior mean/median only.
#'
#' @return a list containing posterior mean/median/95% CI of fitted prob at given x
#' @export
#'

predictCov = function(model, grid.eq = T, node.t, newX, newZ, returnCI = T){

  rs.sampler = rstan::extract(model, pars = c("alpha", "gamma", "Gamma","denominator","cumsum_alphaW"))
  L = ncol(rs.sampler$alpha)
  M.samples = nrow(rs.sampler$alpha)
  if (grid.eq == T){
    node.t = seq(0,1, length.out = L + 1)
  }
  x.J = rowSums(outer(newX, node.t, "-")>=0)
  x.J[x.J > L] = L # Extrapolation for new X greater than the maximum X for model fitting
  x.J[x.J == 0] = 1  # Extrapolation for new X less than the minimum X for model fitting

  y = matrix(rep(newX, each = M.samples), nrow = M.samples, ncol = length(newX))
  for (i in 1:length(newX)){
    j = x.J[i]
    if (j == 1){
      y[,i] = ( newX[i] * rs.sampler$alpha[,1] + 1 ) / rs.sampler$denominator
    } else {
      y[,i] = ( rs.sampler$cumsum_alphaW[,j-1] + (newX[i]-node.t[j]) * rs.sampler$alpha[,j] + 1  ) / rs.sampler$denominator
    }
    y[,i] = logit(y[,i])
    y[,i] = y[,i] + newZ[i] * rs.sampler$Gamma
    y[,i] = expit(y[,i])
  }

  y.mean = colMeans(y)
  y.CI = apply(y, 2, quantile,probs = c(0.5, 0.025, 0.975), na.rm =T)
  y.all <-list("mean" = y.mean, "median" = y.CI[1,],
               "l.ci" = y.CI[2,], "u.ci" = y.CI[3,])

  return(y.all)
}

