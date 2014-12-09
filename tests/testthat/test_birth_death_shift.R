library(BirthDeathBirth)
library(compiler)
library(deSolve)

test_that("Birth death process", {
  tolerance = 1E-7
  
  tList = 1;  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
  gridLength = 51
  s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
  s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
  print("Vladimir's code:")
  p0 <- getTrans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)[[1]]
  p1 = p0[1:11,]
  
  brates1=function(a,b){0}
  drates1=function(a,b){mu*a}
  brates2=function(a,b){lam*(a+b)}
  drates2=function(a,b){mu*b}
  trans=function(a,b){v*a}
  
  a0 = 10
  b0 = 0
  A = 0
  B = 50
  
  print("non-parallel bbd_prob:")
  p2 <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B)
  print("parallel bbd_prob:")
  p3 <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1)
  p4 <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=2)
  
  expect_equal(0.0, sum(abs(p1-p2)), tolerance)
  expect_equal(0.0, sum(abs(p2-p3)), tolerance)
  expect_equal(0.0, sum(abs(p3-p4)), tolerance)
  
})