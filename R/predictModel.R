#' Title
#'Mar 6: this is corrected version of predict function: return posterior mean/median and 95% crediable interval
#'
#' @param model Fitted Model
#' @param L the pre-defined L for model
#' @param method posteiror mean or median used
#' @param grid.eq equal grid, if F need to provide node.t below
#' @param node.t default equal nodes on 0-1
#' @param newdata newdata
#'
#' @return a vector of fitted prob at given x
#' @export
#'

predictModel = function(model,L,method = "mean", grid.eq = T, node.t, newdata){
  y = x = newdata
  if (method == "mean") {
    yi = rstan::summary(model, par =  "xi")$summary[,"mean"]
  } else if (method == "median") {
    yi = rstan::summary(model, par =  "xi")$summary[,"50%"]
  }
  M = length(yi) - 1
  if (grid.eq == T){
    node.t = seq(0,1, length.out = M + 1)
  }
  x.J = x %/% (1/M) + 1
  for (i in 1:length(x)){
    j = x.J[i]
    y[i] = yi[j] + ( yi[j+1] - yi[j] )/(node.t[j+1] - node.t[j]) * (x[i] - node.t[j])
  }
  return(y)
}
