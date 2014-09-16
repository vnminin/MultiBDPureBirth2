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

brates1=function(x,y){beta*x*y/N}
brates2=function(x,y){beta*x*y/N}
drates2=function(x,y){gamma*y}

system.time(p <- bbd_prob(t=1,a0,b0,brates1,brates2,drates2,A,B))

bbd_phi(s=1,a0,b0,brates1,brates2,drates2,A,B)
  
  source("bd_prob.r")
  source("bd_lt.r")
  
  states = 0:50
  p1 = sapply(states, bd_prob, m=3, t=1, brates=function(k){0.5*k}, drates=function(k){0.3*k})
  plot(states,p1)
  
  p = sapply(states,bbd_prob,t=1,a0=0,b0=3,brates1=function(a,b){0.3*a+0.1},brates2=function(a,b){0.5*b},drates2=function(a,b){0.3*b},1,5,0) 
