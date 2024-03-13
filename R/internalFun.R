#' logit
#'
#' @param x value
#' @noRd
logit = function(x){
  return(qlogis(x))
}

#' expit
#'
#' @param x value
#' @noRd
expit = function(x){
  return(plogis(x))
}
