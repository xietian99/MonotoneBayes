#' sumStats
#'
#' @param model rstan object, including the result from monotoneBayes()
#' @param AllPara default F, if T then return all parameter in the model including c,tau,and lambda
#' @param ndigit default 3, the digits of results
#' @param orgHS default F, only T when the model is fitted via original HS method
#'
#' @return The summary statistics of fitted model
#' @export
#'
#'@importFrom stats median

sumStats = function(model, AllPara = F, ndigit = 3, orgHS = F){ # find the posterior distribution of Gamma, IntAlpha, alpha , lambda, and inte
  if (AllPara == T){
    if (orgHS == T){ # then no c_sq parameter
      outcome = round(rstan::summary(model, pars =c("gamma","tau_sq","alpha","lambda_sq"))$summary, ndigit)
    } else {
      outcome = round(rstan::summary(model, pars =c("gamma","c_sq","tau_sq", "alpha","lambda_sq"))$summary, ndigit)
    }
  } else {
    outcome.sum = round(rstan::summary(model, probs = c(0.5, 0.05, 0.25, 0.75, 0.95), pars =c("gamma", "alpha", "xi"))$summary, ndigit)
    outcome = rbind(outcome.sum[1,], 0, outcome.sum[-1,])
    outcome[2,] = colSums(outcome.sum[-1,])
    outcome[2,c(2,3,5:10)] = NA
    rownames(outcome)[1:2] = c("gamma","sumAlpha")
  }
  return(outcome)
}



#' plotModel
#'
#' @param model rstan object, including the result from monotoneBayes()
#' @param L user-defined interval
#' @param method fitted curve use either posteiror "mean" or "median"
#' @param grid.eq default T for equal spacing
#' @param node.t if grid.eq == T, equally spaced nodes at 0-1, or when grid.eq = F, this should be user-defined node point
#'
#' @return a list
#' $FittedProb: a data frame list the estimated value at specific points of [0,1]
#' $plots: the fitted model curve
#' @export
#'
#'
#' @importFrom rstan extract
#' @importFrom ggplot2 ggplot geom_line ylim
plotModel = function(model,L,method = "mean", grid.eq = T, node.t){
  y = x = seq(0, 1, length.out = 100)
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
  df = data.frame("x" = x, "EstPr.Y" = y)
  plt = ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x = x, y = y)) + ggplot2::ylim(0,1)
  rs = list(df, plt)
  names(rs) = c("FittedProb", "Curve")
  return(rs)
}


