library(BirthDeathBirth)
library(compiler)

test_that("Sum of conditional probabilities", {
  tolerance = 1E-8
  
  a0 = 10
  b0 = 0
  A = 0
  B = a0
  
  muL = runif(1,0,1)
  muM = 0.0015
  eta = runif(1,0,1)
  gamma = 0.04
  
  drates1=function(a,b){muL*a+eta*a^2}
  brates2=function(a,b){0}
  drates2=function(a,b){muM*b}
  trans=function(a,b){gamma*a} # a -> b
  
  #Rprof("func.out",memory.profiling=T)
#   system.time(
    p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B)
#     )
  #p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B)
  #Rprof(NULL)
  #summaryRprof("func.out",memory="both")
  
  
  expect_equal(1.0, sum(p), tolerance)
  
})