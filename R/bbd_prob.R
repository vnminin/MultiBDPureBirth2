###########################################################
### Transition probabilities of a birth/birth-death process
###########################################################

#' Transition probabilities of a birth/birth-death process
#'
#' Computes the transition pobabilities of a birth/birth-death process 
#' using the continued fraction representation of its Laplace transform
#' @param t time
#' @param a0 total number of type 1 particles at \code{t = 0}
#' @param b0 total number of type 2 particles at \code{t = 0}
#' @param lambda1 birth rate of type 1 particles (a two variables function)
#' @param lambda2 birth rate of type 2 particles (a two variables function)
#' @param mu2 death rate function of type 2 particles (a two variables function)
#' @param gamma transition rate from type 2 particles to type 1 particles (a two variables function)
#' @param A upper bound for the total number of type 1 particles
#' @param B upper bound for the total number of type 2 particles
#' @param nblocks number of blocks
#' @param tol tolerance
#' @param computeMode computation mode
#' @param nThreads number of threads
#' @param maxdepth maximum number of iterations for Lentz algorithm
#' @export
bbd_prob <- function(t, a0, b0, lambda1, lambda2, mu2, gamma, A, B,
                     nblocks=256, tol=1e-12, computeMode=0, nThreads=4,
                     maxdepth=400) {
  
#   setThreadOptions(numThreads = nThreads)
  
  ###################
  ### Input checking
  ###################
  
  if (a0 < 0) stop("a0 cannot be negative.")
  if (a0 > A) stop("a0 cannot be bigger than A.")
  if (b0 < 0) stop("b0 cannot be negative.")
  if (b0 > B) stop("b0 cannot be bigger than B.")
  if (t < 0) stop("t cannot be negative.")
  
  ################################
  ### t is too small
  ### set probability 1 at a0, b0
  ################################
  
  if (t < tol) {
    res = matrix(0, nrow=A-a0+1, ncol=B+1)
    res[1,b0+1] = 1
    colnames(res) = 0:B
    rownames(res) = a0:A
    return(res)
  }
  
  ##################################
  ### store rate values in matrices
  ##################################
  
  grid  = expand.grid(a0:A,0:(B+maxdepth))
  l1 = matrix(mapply(lambda1, grid[,1], grid[,2]), ncol=B+1+maxdepth)
  l2 = matrix(mapply(lambda2, grid[,1], grid[,2]),ncol=B+1+maxdepth)
  m2 = matrix(mapply(mu2, grid[,1], grid[,2]), ncol=B+1+maxdepth)
  g = matrix(mapply(gamma, grid[,1], grid[,2]), ncol=B+1+maxdepth)
  
  ##########################################################
  ### store x and y values (continued fraction) in matrices
  ##########################################################
  
  xf <- function(a,b){
    if (b==0) x = 1
    else x = - l2[a-a0+1,b]*m2[a-a0+1,b+1]   
    return(x)
  }
  x = matrix(mapply(xf, grid[,1], grid[,2]), nrow=B+1+maxdepth, byrow=TRUE)
  
  yf <- function(a,b) {
    y = l1[a-a0+1,b+1] + l2[a-a0+1,b+1] + m2[a-a0+1,b+1] + g[a-a0+1,b+1]
    return(y)
  }
  y = matrix(mapply(yf, grid[,1], grid[,2]), nrow=B+1+maxdepth, byrow=TRUE)	
  
  #############################
  ### Call C function via Rcpp
  #############################
  
  res = matrix(bbd_lt_invert_Cpp(t, a0, b0, l1, l2, m2, g, x, y, A, B+1,
                                 nblocks, tol, computeMode, nThreads, maxdepth),
               nrow=(A-a0+1),byrow=T)
  
  colnames(res) = 0:B
  rownames(res) = a0:A
  
  return(abs(res))
}
