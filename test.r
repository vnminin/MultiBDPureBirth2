source("bd_prob.r")
source("bbd_prob.r")
source("bbd_lt.r")
source("bd_estep.r")
source("bd_lt.r")
source("bd_series_accel.r")
source("contfr.r")

p = matrix(0,ncol=4,nrow=4)
 
for (i in 0:3)
for (j in 0:3) 
p[i+1,j+1] = bbd_prob(t=1,a0=0,b0=0,brates1=function(a,b){0.3*a+0.1},brates2=function(a,b){0.3*b+0.1},drates2=function(a,b){0.5*b+0.1},3,3,i,j) 
  #plot(states,p)
  
  
  states = 0:5
  p1 = sapply(states, bd_prob, m=3, t=1, brates=function(k){0.5*k}, drates=function(k){0.3*k})
  plot(states,p1)
  
  p = sapply(states,bbd_prob,t=1,a0=0,b0=3,brates1=function(a,b){0.3*a+0.1},brates2=function(a,b){0.5*b},drates2=function(a,b){0.3*b},1,5,0) 
