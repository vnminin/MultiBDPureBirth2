# Functions for simulating a BDP starting at a certain state
# and ending at a certain time.  
#
# Author: Forrest W. Crawford

########################################

bd_simulate_continuous = function(start,tmax,brates,drates) {

  x = rep(NA,100)
  t = rep(NA,100)
  tcur = 0
  x[1] = start
  i = 1

  while(1) {
    u = runif(1)
    lambda = brates(x[i])
    mu = drates(x[i])

    if((mu == 0)&&(lambda==0)) {
      if(i>1) t[i] = tmax-tcur
      else t[i] = tmax
      break
    }

    if(u < lambda/(lambda+mu))  xnew = 1 
    else xnew = -1 

    tnew = rexp(1,lambda+mu)

    if(tcur+tnew > tmax) {
      t[i] = tmax - tcur
      break
    } else {

      t[i] = tnew
      tcur = tcur + tnew
      x[i+1] = x[i] + xnew
      i = i + 1
    }
  }
  return(cbind(x[1:i], t[1:i]))
}


########################################################

bd_simulate_discrete = function(start,tmax,brates,drates) {
  x = bd_simulate_continuous(start,tmax,brates,drates)
  return(x[length(x[,1]),1])
}


