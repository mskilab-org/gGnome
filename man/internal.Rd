% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rcplex.R
\name{Rcplex2}
\alias{Rcplex2}
\title{Rcplex2}
\usage{
Rcplex2(
  cvec,
  Amat,
  bvec,
  Qmat = NULL,
  lb = 0,
  ub = Inf,
  control = list(),
  objsense = c("min", "max"),
  sense = "L",
  vtype = NULL,
  n = 1,
  tuning = FALSE
)
}
\description{
Modification of Rcplex which takes in mipcontrol parameters.
}
