#' Title
#'Mar 6: this is corrected version of predict function: return posterior mean/median and 95% crediable interval
#'
#' @param model Fitted Model
#' @param grid.eq equal grid, if F need to provide node.t below
#' @param node.t default equal nodes on 0-1
#' @param newdata newdata
#'
#' @return a list containing posteiror mean, median and 95% C.I. of fitted prob at given x
#' @export
#'

predictModel = function(model, grid.eq = T, node.t, newdata){
  rs.sampler = rstan::extract(model, pars = c("alpha", "gamma", "denominator","cumsum_alphaW"))
  L = ncol(rs.sampler$alpha)
  M.samples = nrow(rs.sampler$alpha)
  if (grid.eq == T){
    node.t = seq(0,1, length.out = L + 1)
  }
  x.J = rowSums(outer(newdata, node.t, "-")>=0)
  x.J[x.J > L] = L # Extrapolation for new X greater than the maximum X for model fitting
  x.J[x.J == 0] = 1  # Extrapolation for new X less than the minimum X for model fitting

  y = matrix(rep(newdata, each = M.samples), nrow = M.samples, ncol = length(newdata))
  for (i in 1:length(newdata)){
    j = x.J[i]
    if (j == 1){
      y[,i] = ( newdata[i] * rs.sampler$alpha[,1] + 1 ) / rs.sampler$denominator
    } else {
      y[,i] = ( rs.sampler$cumsum_alphaW[,j-1] + (newdata[i]-node.t[j]) * rs.sampler$alpha[,j] + 1  ) / rs.sampler$denominator
    }
  }

  y.mean = colMeans(y)
  y.CI = apply(y, 2, quantile,probs = c(0.5, 0.025, 0.975), na.rm =T)
  y.all <-list("mean" = y.mean, "median" = y.CI[1,],
               "l.ci" = y.CI[2,], "u.ci" = y.CI[3,])

  return(y.all)
}
