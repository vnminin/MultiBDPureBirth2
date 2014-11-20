library(BirthDeathBirth)
library(parallel)
library(compiler)
library(Rcpp)

### Within-host macroparasite population

a0 = 20
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
system.time(p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B))
p[1:10,1:5]
#p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B)
#Rprof(NULL)
#summaryRprof("func.out",memory="both")
sum(p)
#print(c(muL,muM))


### SIR

a0 = 1
A = 0
B = 3
N = 3
b0 = N-a0

gamma = 1
beta = pi

brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){gamma*b}
trans=function(a,b){beta*a*b/N}

# system.time(p <- bbd_prob(t=1,a0,b0,brates1,brates2,drates2,trans,A,B))
Rprof("func.out",memory.profiling=T)
#system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,B))
p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B)
Rprof(NULL)
summaryRprof("func.out",memory="both")

#source("bbd_prob0.r")
#source("bbd_lt0.r")
#source("contfr0.r")

#fun <-function(x,y) {return(dbd_prob0(t=1,a0,b0,drates1,brates2,drates2,trans,x,y,B))}
#grid = expand.grid(0:a0,0:B)
#system.time(p0 <- matrix(mapply(fun,grid[,1],grid[,2]),ncol=B+1))

source("sir.r")
#sir(t=1,n=a0,a=N-a0,f=trans,mu=gamma,0,0)

funsir <- function(x,y) {return(sir(t=1,n=a0,a=N-a0,f=trans,mu=gamma,x,y))}
grid = expand.grid(0:a0,0:N)
system.time(p1 <- matrix(mapply(funsir,grid[,1],grid[,2]),nrow = a0+1))


bbd_phi(s=1,a0,b0,brates1,brates2,drates2,trans,A,B)
  
  source("bd_prob.r")
  source("bd_lt.r")
  
  states = 0:50
  #Rprof("func.out",memory.profiling=T)
  system.time(p1 <- sapply(states, bd_prob, m=3, t=1, brates=function(k){return(0.5*k)}, drates=function(k){0.3*k}))
  #Rprof(NULL)
#summaryRprof("func.out",memory="both")

 #Rprof("func.out",memory.profiling=T)
 system.time(p <- bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50))
  sum(p)
 #Rprof(NULL)
#summaryRprof("func.out",memory="both")

 
###########################################

# Data: Spread of Smallpox in a Nigerian Village (Yip 1989 - Theoretical Pop Bio)
# This data is not appropriate for SIR model

t=0:83
i = c(rep(1,23),rep(2,3),5,6,5,5,4,5,5,2,1,1,2,2,1,2,2,4,
	4,5,5,5,4,4,3,3,1,2,3,3,3,2,4,5,5,5,5,7,8,6,5,4,3,5,3,2,
	2,2,3,3,1,1,1,2,2,1,1,1,1,1)
s = c(119,rep(118,7),117,117,116,116,116,113,rep(112,4),rep(111,5),
	110,110,110,109,109,107,107,rep(105,5),104,104,104,103,rep(102,4),
	100,99,98,97,97,95,rep(94,5),rep(92,5),rep(91,5),rep(90,20))
	
N = 120

gamma = 10
beta =  10
brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){beta*b}
trans=function(a,b){alpha*a*b}

loglik <- function(i,s,drates1,brates2,drates2,trans) {
	loglik = 0
	n = length(i)
	for (k in 1:(n-1)) {
		p = dbd_prob(t=1,a0=s[k],b0=i[k],drates1,brates2,drates2,trans,
			a=s[k+1],B=s[k]+i[k]-s[k+1])
		loglik = loglik + log(p[s[k]-s[k+1]+1,i[k+1]])	
		print(c(loglik,s[k],i[k]))
	}
	return(loglik)
}

loglik(i,s,drates1,brates2,drates2,trans)


lik = matrix(NA,nrow=10,ncol=10)
for (x in 1:10)
for (y in 1:10){
	alpha = x/10
	beta =  y/10
	brates1=function(a,b){0}
	drates1=function(a,b){0}
	brates2=function(a,b){0}
	drates2=function(a,b){beta*b}
	trans=function(a,b){alpha*a*b}
	lik[x,y] = loglik(N,i,s,drates1,brates2,drates2,trans)
	print(c(x,y))
}  

#########################################################
# The Great Plague in Eyam 
# Mathemetical Epidemiology (2008) - Fred Brauer et al.

# t = c(0,1.5,2,2.5,3,3.5,4,4.5)*30
# s = c(254,235,201,154,121,108,97,83)
# i = c(7,15,22,29,21,8,8,0)

t = c(0,16,17,16,17,32)
s = c(201,154,121,108,97,83)
i = c(22,29,21,8,8,0)

## Independence sorce
## Date - S - I
## June 18 - 254 - 7
## July 3/4 - 235 - 14.5
## July 19 - 201 - 22
## August 3/4 - 153.5 - 29
## August 19 - 121 - 20
## Sep 3/4 - 110 - 8
## Sep 19 - 97 - 8
## Oct 4/5 - Unknown - Unknown
## Oct 20 - 83 - 0


### Likelihood
loglik <- function(param) {

  alpha = param[1]
  beta = param[2]
  
  if ((alpha<0)||(beta<0)) return(-Inf)
  
  brates1=function(a,b){0}
  drates1=function(a,b){0}
  brates2=function(a,b){0}
  drates2=function(a,b){alpha*b}
  trans=function(a,b){beta*a*b}
  
	n = length(i)
  fun <- function(k){return(log(dbd_prob(t=t[k+1]-t[k],a0=s[k],b0=i[k],drates1,brates2,drates2,trans,
                                     a=s[k+1],B=s[k]+i[k]-s[k+1]))[1,i[k+1]+1])}
  #tmp = mclapply(1:(n-1),fun,mc.cores=3)
  tmp = sapply(1:(n-1),fun)
  
  #loglik = sum(unlist(tmp))
  loglik = sum(tmp)
	return(loglik)
}


alpha = 2.73
beta = 0.0178
brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){alpha*b}
trans=function(a,b){beta*a*b}

Rprof("func.out",memory.profiling=T)
p <- dbd_prob(t=15,a0=235,b0=15,drates1,brates2,drates2,trans,a=201,B=49)  
#p <- dbd_prob(t=15,a0=235,b0=5,drates1,brates2,drates2,trans,a=220,B=20)  
#sum(p)
Rprof(NULL)
summaryRprof("func.out",memory="both")

# system.time(l<-loglik(c(alpha,beta)))
# print(c(l,alpha,beta))


### Prior
logprior <- function(param){
  alpha = param[1]
  beta = param[2]
  aprior = dgamma(alpha, shape = 2.73, log = T)
  bprior = dgamma(beta, shape = 0.0178, log = T)
  return(aprior+bprior)
}

posterior <- function(param){
  return (loglik(param) + logprior(param))
}

######## Metropolis algorithm ################

proposalfunction <- function(param){
  return(rnorm(2,mean = param, sd= c(0.05,0.005)))
  # small sd, more acceptance
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    print(probab)
    print("--")
    print(chain[i,])
    print(proposal)
    print("--")
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    print(i)
  }
  return(chain)
}

#alpha = 2.73
#beta =  0.0178

alpha = 2.93516604664236
beta = 0.0174758550488344
startvalue = c(alpha,beta)

#Rprof("func.out",memory.profiling=T)
chain <- run_metropolis_MCMC(startvalue, 200)
#Rprof(NULL)
#summaryRprof("func.out",memory="both")

burnIn = 20
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
acceptance
plot(chain[,1],type="l")
plot(chain[,2],type="l")
hist(chain[,1],breaks=20)
hist(chain[,2],breaks=20)

write.table(chain, "tests/chain.txt", col.names=F, row.names=F, append = T)

# minus.loglik <- function(parameter){return(-loglik(parameter))}
# optim(c(2.73,0.0178), minus.loglik, method="L-BFGS-B", lower=c(2,0.01), upper=c(3,0.03))

dat = read.table("tests/chain.txt",header=F)
plot(dat[,1],type="l")
plot(dat[,2],type="l")
hist(dat[,1],breaks=20)
hist(dat[,2],breaks=20)
