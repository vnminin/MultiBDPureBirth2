source("bbd_prob.r")
source("bbd_lt.r")
source("bd_series_accel.r")
source("contfr.r")

a0 = 0
b0 = 0
A = 10
B = 10
p = matrix(0,ncol=A+1,nrow=B+1)

brates1=function(a,b){0.3*a+0.01}
brates2=function(a,b){0.3*a+0.01}
drates2=function(a,b){0.5*b+0.01}

for (i in 0:A)
for (j in 0:B)
p[i+1,j+1] = bbd_prob(t=1,a0,b0,brates1,brates2,drates2,A,B,a=i,b=j)

bbd_phi(s=1,a0,b0,brates1,brates2,drates2,A,B)
  
  states = 0:5
  p1 = sapply(states, bd_prob, m=3, t=1, brates=function(k){0.5*k}, drates=function(k){0.3*k})
  plot(states,p1)
  
  p = sapply(states,bbd_prob,t=1,a0=0,b0=3,brates1=function(a,b){0.3*a+0.1},brates2=function(a,b){0.5*b},drates2=function(a,b){0.3*b},1,5,0) 
