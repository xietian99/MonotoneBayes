#' logit
#'
#' @param x value
#' @noRd
logit = function(x){
  return(log(x/(1-x)))
}

#' expit
#'
#' @param x value
#' @noRd
expit = function(x){
  return(exp(x)/(1+exp(x)))
}
