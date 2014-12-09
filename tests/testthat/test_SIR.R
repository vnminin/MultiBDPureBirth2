library(BirthDeathBirth)
library(compiler)
library(expm)

test_that("SIR", {
  tolerance = 1E-8
  
  a0 = 10
  A = 3
  N = 20
  B = N
  b0 = N-a0
  t = 1
  
  alpha = 2.73
  beta = 0.0178
  brates1=function(a,b){0}
  drates1=function(a,b){0}
  brates2=function(a,b){0}
  drates2=function(a,b){alpha*b}
  trans=function(a,b){beta*a*b}
  
  Q = matrix(0,nrow = (a0-A+1)*(B+1), ncol = (a0-A+1)*(B+1))
  P0 = rep(0,(a0-A+1)*(B+1))
  
  for (i in A:a0) {
    for (j in 0:B) {
      tmp = (i-A)*(B+1) + j+1
      Q[tmp,tmp] = - alpha*j - beta*i*j
      if (j>0) {
        tmp1 = (i-A)*(B+1) + j
        Q[tmp,tmp1] = alpha*j
      }
      if ((i>A)&&(j<B)) {
        tmp2 = (i-A-1)*(B+1) + j+2
        Q[tmp,tmp2] = beta*i*j
      }
    }
  }
  tmp0 = (a0-A)*(B+1) + b0+1
  P0[tmp0] = 1
  
  expM = expAtv(t(Q),P0,t)$eAtv
  P = matrix(0,nrow=a0-A+1, ncol=B+1)
  
  for (i in A:a0) {
    for (j in 0:B) {
      P[i-A+1,j+1] = expM[(i-A)*(B+1) + j+1]
    }
  }
    
  print("non-parallel bbd_prob:")
  p1 <- dbd_prob(t,a0,b0,drates1,brates2,drates2,trans,a=A,B)
  
  print("parallel bbd_prob:")
  p2 <- dbd_prob(t,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1)
  p3 <- dbd_prob(t,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=2)
  
  expect_equal(0.0, sum(abs(P-p1)), tolerance)
  expect_equal(0.0, sum(abs(p2-p1)), tolerance)
  expect_equal(0.0, sum(abs(p3-p2)), tolerance)
})