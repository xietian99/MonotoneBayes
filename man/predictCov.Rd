% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictCov.R
\name{predictCov}
\alias{predictCov}
\title{Currently this predict function works only for cases with new coavriates.}
\usage{
predictCov(model, grid.eq = T, node.t, newX, newZ, returnCI = T)
}
\arguments{
\item{model}{fitted stan model}

\item{grid.eq}{equal grid, if F need to provide node.t below}

\item{node.t}{default equal nodes on 0-1}

\item{newX}{input new dose, should be transformed the same way with previous dose, could be extrapolation}

\item{newZ}{input new covariates}

\item{returnCI}{if T will return 95\% Crediable Interval, otherwise return posterior mean/median only.}

\item{method}{mean or median}
}
\value{
a list containing posterior mean/median/95\% CI of fitted prob at given x
}
\description{
Currently this predict function works only for cases with new coavriates.
}
