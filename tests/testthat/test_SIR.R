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
  
  alpha = 2.73
  beta = 0.0178
  brates1=function(a,b){0}
  drates1=function(a,b){0}
  brates2=function(a,b){0}
  drates2=function(a,b){alpha*b}
  trans=function(a,b){beta*a*b}
  
  grid = expand.grid(A:a0,0:B)
  Q = matrix(0,nrow = nrow(grid), ncol = nrow(grid))
  for (i in 1: nrow(Q))
    for (j in 1: ncol(Q))
    {
      if ((grid[i,1]==grid[j,1])&&(grid[i,2] == grid[j,2]+1)) Q[i,j] = alpha*grid[i,2]
      if ((grid[i,1]==grid[j,1]+1)&&(grid[i,2] == grid[j,2]-1)) Q[i,j] = beta*grid[i,1]*grid[i,2]
      if ((grid[i,1]==grid[j,1])&&(grid[i,2] == grid[j,2])) Q[i,j] = - alpha*grid[i,2] - beta*grid[i,1]*grid[i,2]
    }
  
  P0 = rep(0,nrow(grid))
  for (i in 1:nrow(grid)) {if ((grid[i,1]==a0)&&(grid[i,2]==b0)) P0[i]=1}
  tmp = expAtv(t(Q),P0,t=15)$eAtv
  P = matrix(0,nrow=a0-A+1, ncol=B+1)
  for (i in 1:nrow(grid)) {P[grid[i,1]-A+1,grid[i,2]+1]=tmp[i]}
  
  print("non-parallel bbd_prob:")
  p1 <- dbd_prob(t=15,a0,b0,drates1,brates2,drates2,trans,a=A,B)
  
  print("parallel bbd_prob:")
  p2 <- dbd_prob(t=15,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1)
  
  expect_equal(0.0, sum(abs(P-p1)), tolerance)
  expect_equal(0.0, sum(abs(p2-p1)), tolerance)
})