source("bbd_prob.r")
source("bbd_lt.r")
source("bd_series_accel.r")
source("contfr.r")

a0 = 3
b0 = 3
A = 20
B = 20
p = matrix(0,ncol=A+1,nrow=B+1)
N = 20
gamma = 1
beta = 1

brates1=function(x,y){0}
drates1=function(x,y){0}
brates2=function(x,y){0}
drates2=function(x,y){gamma*y}
trans=function(x,y){beta*x*y/N}

system.time(p <- bbd_prob(t=1,a0,b0,brates1,brates2,drates2,trans,A,B))

system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,B))

bbd_phi(s=1,a0,b0,brates1,brates2,drates2,A,B)
  
  source("bd_prob.r")
  source("bd_lt.r")
  
  states = 0:50
  p1 = sapply(states, bd_prob, m=3, t=1, brates=function(k){0.5*k}, drates=function(k){0.3*k})
  
  p = bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50)
  
  