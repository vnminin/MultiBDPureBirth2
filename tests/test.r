library(BirthDeathBirth)
library(parallel)
library(compiler)
library(Rcpp)
library(Matrix)
library(expoRkit)

#### Time comparision

a0 = 100
b0 = 0
A = 0
B = a0

muL = 6.82
muM = 0.0015
eta = 0.09
gamma = 0.04

drates1=function(a,b){muL*a+eta*a^2}
brates2=function(a,b){0}
drates2=function(a,b){muM*b}
trans=function(a,b){gamma*a} # a -> b

system.time(p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B,
                          computeMode=4,nblocks=4, nThreads = 4))
sum(p)
system.time(p1 <- dbd_expM(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B))
sum(abs(p-p1))

sim = 100
t = rep(NA,sim)
for (i in 1:sim) {
  t[i] = system.time(p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=0))[3]    
}
mean(t)
sd(t)

###########
### Birth-death-shift
library(compiler)
library(deSolve)
library(expoRkit)
library(Matrix)

initNum = 30
gridLength = 81

tList = 1;  dt = 1; lam = .5; v = .2; mu = .4; 
s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)

brates1=function(a,b){0}
drates1=function(a,b){mu*a}
brates2=function(a,b){lam*(a+b)}
drates2=function(a,b){mu*b}
trans=function(a,b){v*a}

a0 = initNum
b0 = 0
A = 0
B = gridLength - 1

# system.time(p0 <- getTrans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)[[1]])
# p1 = p0[1:(initNum+1),]
system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B))
sum(p)
# sum(abs(p-p1))
system.time(p2 <- dbd_expM(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B))
# sum(p2)
sum(abs(p-p2))


### Within-host macroparasite population

expMtime = rep(0,10)
dbdtime = rep(0,10)

for (sim in 1:10) {
a0 = 100
# a0 = 10*sim
b0 = 0
A = 0
B = a0

# nsim = 100
# ser = rep(0,nsim)
# par2 = rep(0,nsim)
# par4 = rep(0,nsim)

# for (i in 1:100) {
  muL = runif(1,0,1)
  muM = 0.0015
  eta = runif(1,0,1)
  gamma = 0.04
  
  drates1=function(a,b){muL*a+eta*a^2}
  brates2=function(a,b){0}
  drates2=function(a,b){muM*b}
  trans=function(a,b){gamma*a} # a -> bs

Q = matrix(0,nrow = (a0-A+1)*(B+1), ncol = (a0-A+1)*(B+1))
P0 = rep(0,(a0-A+1)*(B+1))

for (i in A:a0) {
  for (j in 0:B) {
    tmp = (i-A)*(B+1) + j+1
    Q[tmp,tmp] = - muL*i - eta*i^2 - muM*j - gamma*i
    if (i>A) {
      tmp1 = (i-A-1)*(B+1) + j+1
      Q[tmp,tmp1] =  muL*i+ eta*i^2
    }
    if (j>0) {
      tmp1 = (i-A)*(B+1) + j
      Q[tmp,tmp1] =  muM*j
    }
    if ((i>A)&&(j<B)) {
      tmp2 = (i-A-1)*(B+1) + j+2
      Q[tmp,tmp2] =  gamma*i
    }
  }
}
tmp0 = (a0-A)*(B+1) + b0+1
P0[tmp0] = 1

t = 1
expMtime[sim] = system.time(expM <- expAtv(t(Q),P0,t)$eAtv)[3]
P = matrix(0,nrow=a0-A+1, ncol=B+1)
for (i in A:a0) {
  for (j in 0:B) {
    P[i-A+1,j+1] = expM[(i-A)*(B+1) + j+1]
  }
}
dbdtime[sim] = system.time(p <- dbd_prob(t,a0,b0,drates1,brates2,drates2,trans,a=A,B))[3]
print(sum(abs(p-P)))
}

pdf("expM_vs_bbd.pdf")
N = (1:10)*10
plot(expMtime~N,type="b",col="red",ylab="Time (seconds)",log="y",xlab="Number of lavae",ylim=c(1e-3,1000))
points(dbdtime~N,type="b",col="blue")
legend("bottomright",c("expM","bbd"),col=c("red","blue"),lty=1)
dev.off()
  #Rprof("func.out",memory.profiling=T)
#   ser[i] = system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=0))[3]
#   par2[i] = system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1,nThreads=2))[3]
#   par4[i] = system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=1,nThreads=4))[3]
  #p[1:10,1:5]
  system.time(p <- dbd_prob(t=400,a0,b0,drates1,brates2,drates2,trans,a=A,B,computeMode=2))
  #Rprof(NULL)
  #summaryRprof("func.out",memory="both")
  #sum(p)
  #print(c(muL,muM))
# }

# dat = data.frame(
#   ser = ser,
#   par2 = par2,
#   par4 = par4
#   )
# boxplot(dat, at = c(1,2,4),ylab="Time",xlab="nThreads")

### SIR

a0 = 50
A = 0
N = 100
B = N
b0 = N-a0

gamma = 1
beta = 1

brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){gamma*b}
trans=function(a,b){beta*a*b}

# system.time(p <- bbd_prob(t=1,a0,b0,brates1,brates2,drates2,trans,A,B))
#Rprof("func.out",memory.profiling=T)
system.time(p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B))
sum(p)
system.time(p1 <- dbd_expM(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B))
sum(abs(p-p1))
# p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,a=A,B)
#Rprof(NULL)
#summaryRprof("func.out",memory="both")

#source("bbd_prob0.r")
#source("bbd_lt0.r")
#source("contfr0.r")

# fun <-function(x,y) {return(dbd_prob0(t=1,a0,b0,drates1,brates2,drates2,trans,x,y,B))}
# grid = expand.grid(0:a0,0:B)
# system.time(p0 <- matrix(mapply(fun,grid[,1],grid[,2]),ncol=B+1))

#source("sir.r") 
# bbd_phi(s=1,a0,b0,brates1,brates2,drates2,trans,A,B)
  
#   source("bd_prob.r")
#   source("bd_lt.r")
#   
  states = 0:50
  #Rprof("func.out",memory.profiling=T)
  system.time(p1 <- sapply(states, bd_prob, m=3, t=1, brates=function(k){return(0.5*k)}, drates=function(k){0.3*k}))
  #Rprof(NULL)
#summaryRprof("func.out",memory="both")

 #Rprof("func.out",memory.profiling=T)
  system.time(p <- bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50))
  sum(p)
  system.time(p2 <- bbd_expM(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50))
  sum(abs(p-p2))
 #Rprof(NULL)
#summaryRprof("func.out",memory="both")

system.time(p <- bbd_prob(t=1,0,0,lambda1=function(a,b){return(0.5*(max(3+a-b,0)))},lambda2=function(a,b){return(0.3*(max(3+a-b,0)))},mu2=function(a,b){return(0)},gamma=function(a,b){0},A=50,B=50))
tmp = rep(0,51)
for (i in 1:51)
  for (j in 1:51) {
    if ((3+i-j >= 0)&&(3+i-j <= 50)) tmp[3+i-j+1] = tmp[3+i-j+1] + p[i,j]
  }
sum(abs(p1-tmp))

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

# t = c(0,16,33,49,66,98)
t = c(0,0.5,1,1.5,2,3)
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
loglik <- function(param,method=1) {

  alpha = exp(param[1])
  beta = exp(param[2])
  
  #if ((alpha<0)||(beta<0)) return(-Inf)
  
  brates1=function(a,b){0}
  drates1=function(a,b){0}
  brates2=function(a,b){0}
  drates2=function(a,b){alpha*b}
  trans=function(a,b){beta*a*b}
  
	n = length(i)
  if (method == 1)
    fun <- function(k){return(log(dbd_prob(t=t[k+1]-t[k],a0=s[k],b0=i[k],drates1,brates2,drates2,trans,
                                     a=s[k+1],B=s[k]+i[k]-s[k+1], computeMode = 0, nblocks = 20))[1,i[k+1]+1])}
  else 
    fun <- function(k){return(log(dbd_expM(t=t[k+1]-t[k],a0=s[k],b0=i[k],drates1,brates2,drates2,trans,
                                              a=s[k+1],B=s[k]+i[k]-s[k+1]))[1,i[k+1]+1])}
  #tmp = mclapply(1:(n-1),fun,mc.cores=3)
  tmp = sapply(1:(n-1),fun)
#   tmp1 = sapply(1:(n-1),fun1)
#   print(abs(tmp-tmp1))
  
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
system.time(p <- dbd_prob(t=15,a0=254,b0=7,drates1,brates2,drates2,trans,a=235,B=28))
system.time(p1 <- dbd_expM(t=15,a0=254,b0=7,drates1,brates2,drates2,trans,a=235,B=28))  
sum(abs(p-p1))
#### Example where exponential matrix method fails
p <- dbd_prob(t=15,a0=235,b0=5,drates1,brates2,drates2,trans,a=220,B=20)  
p1 <- dbd_expM(t=15,a0=235,b0=5,drates1,brates2,drates2,trans,a=220,B=20)  
#sum(p)
Rprof(NULL)
summaryRprof("func.out",memory="both")

system.time(l<-loglik(c(log(alpha),log(beta))))
system.time(l1<-loglik(c(log(alpha),log(beta)),method=0))
print(c(l,l1,alpha,beta))


### Prior
logprior <- function(param){
  alpha = param[1]
  beta = param[2]
  aprior = dnorm(alpha, mean = log(2.73), sd = 0.1, log = TRUE)
  bprior = dnorm(beta, mean = log(0.0178), sd = 0.1, log = TRUE)
  return(aprior+bprior)
}

posterior <- function(param){
  return (loglik(param) + logprior(param))
}

######## Metropolis algorithm ################

proposalfunction <- function(param){
  return(rnorm(2,mean = param, sd= c(0.1,0.1)))
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

alpha = 2.73
beta =  0.0178
startvalue = c(log(alpha),log(beta))

#Rprof("func.out",memory.profiling=T)
chain <- run_metropolis_MCMC(startvalue, 20)
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

minus.loglik <- function(parameter){return(-loglik(parameter))}
system.time(opt <- optim(c(log(2.73),log(0.0178)),minus.loglik), hessian = TRUE)

dat = read.table("tests/chain.txt",header=F)
plot(dat[,1],type="l")
plot(dat[,2],type="l")
hist(dat[,1],breaks=20)
hist(dat[,2],breaks=20)


####
# SIR - Matrix exponential

library(expm)
library(expoRkit)
library(Matrix)
a0 = 20
A = 0
B = 50
b0 = 30

alpha = 2.73
beta = 0.0178
brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){alpha*b}
trans=function(a,b){beta*a*b}

# grid = expand.grid(A:a0,0:B)
# Q = matrix(0,nrow = nrow(grid), ncol = nrow(grid))
# for (i in 1: nrow(Q))
#   for (j in 1: ncol(Q))
#   {
#     if ((grid[i,1]==grid[j,1])&&(grid[i,2] == grid[j,2]+1)) Q[i,j] = alpha*grid[i,2]
#     if ((grid[i,1]==grid[j,1]+1)&&(grid[i,2] == grid[j,2]-1)) Q[i,j] = beta*grid[i,1]*grid[i,2]
#     if ((grid[i,1]==grid[j,1])&&(grid[i,2] == grid[j,2])) Q[i,j] = - alpha*grid[i,2] - beta*grid[i,1]*grid[i,2]
#   }
# 
# P0 = rep(0,nrow(grid))
# for (i in 1:nrow(grid)) {if ((grid[i,1]==a0)&&(grid[i,2]==b0)) P0[i]=1}

Q = matrix(0,nrow = (a0-A+1)*(B+1), ncol = (a0-A+1)*(B+1))
P0 = rep(0,(a0-A+1)*(B+1))

for (i in A:a0) {
  for (j in 0:B) {
    tmp = (i-A)*(B+1) + j+1
    Q[tmp,tmp] = - alpha*j - beta*i*j
    if (j>0) {
      tmp1 = (i-A)*(B+1) + j
      Q[tmp,tmp1] =  alpha*j
    }
    if ((i>A)&&(j<B)) {
      tmp2 = (i-A-1)*(B+1) + j+2
      Q[tmp,tmp2] =  beta*i*j
    }
  }
}
tmp0 = (a0-A)*(B+1) + b0+1
P0[tmp0] = 1

t = 1
QQ = as(Q, "TsparseMatrix")
system.time(expM <- expAtv(t(Q),P0,t)$eAtv)
system.time(expMM <- expv(x = QQ,v = P0,t = t,transpose=TRUE))
#system.time(expM0 <- expm(t(Q)*t)%*%P0)
P = matrix(0,nrow=a0-A+1, ncol=B+1)
# for (i in 1:nrow(grid)) {P[grid[i,1]-A+1,grid[i,2]+1]=tmp[i]}
for (i in A:a0) {
  for (j in 0:B) {
    P[i-A+1,j+1] = expM[(i-A)*(B+1) + j+1]
  }
}

system.time(p <- dbd_prob(t,a0,b0,drates1,brates2,drates2,trans,a=A,B))
sum(abs(P-p))

####
# Birth-death-birth

a0 = 3
A = 20
B = 100
b0 = 3

alpha = 1
beta = 1
gamma = 1
brates1=function(a,b){alpha*a}
drates1=function(a,b){0}
brates2=function(a,b){alpha*b}
drates2=function(a,b){beta*b}
trans=function(a,b){gamma*(a+b)}

Q = matrix(0,nrow = (A-a0+1)*(B+1), ncol = (A-a0+1)*(B+1))
P0 = rep(0,(A-a0+1)*(B+1))

for (i in a0:A) {
  for (j in 0:B) {
    tmp = (i-a0)*(B+1) + j+1
    Q[tmp,tmp] = - alpha*(i+j) - beta*j - gamma*(i+j)
    if (i<A) {
      tmp1 = (i-a0+1)*(B+1) + j+1
      Q[tmp,tmp1] =  alpha*i
    }
    if (j<B) {
      tmp1 = (i-a0)*(B+1) + j+2
      Q[tmp,tmp1] =  alpha*j
    }
    if (j>0) {
      tmp1 = (i-a0)*(B+1) + j
      Q[tmp,tmp1] = beta*j
    }
    if ((i<A)&&(j>0)) {
      tmp2 = (i-a0+1)*(B+1) + j
      Q[tmp,tmp2] = gamma*(i+j)
    }
  }
}
tmp0 = b0+1
P0[tmp0] = 1

t = 1e-6
system.time(expM <- expAtv(t(Q),P0,t)$eAtv)
P = matrix(0,nrow=A-a0+1, ncol=B+1)
for (i in a0:A) {
  for (j in 0:B) {
    P[i-a0+1,j+1] = expM[(i-a0)*(B+1) + j+1]
  }
}

system.time(p <- bbd_expM(t,a0,b0,brates1,brates2,drates2,trans,A,B))

print("non-parallel bbd_prob:")
system.time(p1 <- bbd_prob(t,a0,b0,brates1,brates2,drates2,trans,A,B))
sum(abs(p1-p))
sum(abs(p-P))
# print("parallel bbd_prob:")
# system.time(p2 <- bbd_prob(t,a0,b0,brates1,brates2,drates2,trans,A,B,computeMode=1))


