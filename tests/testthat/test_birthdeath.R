library(MultiBD)
library(compiler)

test_that("Birth death process", {
  tolerance = 1E-8
  
  states = 0:50
  print("Forrest's code:")
  system.time(p <- sapply(states, bd_prob, m=3, t=1, brates=function(k){return(0.5*k)}, drates=function(k){0.3*k}))
  print("non-parallel bbd_prob:")
  system.time(p1 <- bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},
                            mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50))
  print("parallel bbd_prob:")
  system.time(p2 <- bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},
                             mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50,computeMode=1))
  system.time(p3 <- bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},
                             mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50,computeMode=2))
  
  expect_equal(0.0, sum(abs(p-p1)), tolerance)
  expect_equal(0.0, sum(abs(p1-p2)), tolerance)
  expect_equal(0.0, sum(abs(p2-p3)), tolerance)
  
})