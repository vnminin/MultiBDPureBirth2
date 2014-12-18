library(BirthDeathBirth)
library(compiler)
library(expoRkit)
library(Matrix)

test_that("SIR", {
  tolerance = 1E-7
  
  alpha = 2.73
  beta = 0.0178
  brates1=function(a,b){0}
  drates1=function(a,b){0}
  brates2=function(a,b){0}
  drates2=function(a,b){alpha*b}
  trans=function(a,b){beta*a*b}
  
  p <- dbd_expM(t=15,a0=235,b0=15,drates1,brates2,drates2,trans,a=201,B=49)
  
  print("non-parallel bbd_prob:")
  p1 <- dbd_prob(t=15,a0=235,b0=15,drates1,brates2,drates2,trans,a=201,B=49)
  
  print("parallel bbd_prob:")
  p2 <- dbd_prob(t=15,a0=235,b0=15,drates1,brates2,drates2,trans,a=201,B=49,computeMode=1)
  p3 <- dbd_prob(t=15,a0=235,b0=15,drates1,brates2,drates2,trans,a=201,B=49,computeMode=2)
  
  expect_equal(0.0, sum(abs(p-p1)), tolerance)
  expect_equal(0.0, sum(abs(p2-p1)), tolerance)
  expect_equal(0.0, sum(abs(p3-p2)), tolerance)
})