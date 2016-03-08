#' MultiBD: Cyclic coordinate descent for logistic, Poisson and survival analysis
#'
#' The MultiBD package incorporates cyclic coordinate descent and
#' majorization-minimization approaches to fit a variety of regression models
#' found in large-scale observational healthcare data.  Implementations focus
#' on computational optimization and fine-scale parallelization to yield
#' efficient inference in massive datasets.
#' 
#' @docType package
#' @name MultiBD
#' @import Rcpp RcppParallel
#' @useDynLib MultiBD
NULL