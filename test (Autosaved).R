source("bbd_prob.r")
source("bbd_lt.r")
source("bd_series_accel.r")
source("contfr.r")
library(parallel)

a0 = 30
A = 50
B = 50
N = 50
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
p <- dbd_prob(t=1,a0,b0,drates1,brates2,drates2,trans,B)
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
  Rprof("func.out",memory.profiling=T)
  system.time(p1 <- sapply(states, bd_prob, m=3, t=1, brates=function(k){0.5*k}, drates=function(k){0.3*k}))
  Rprof(NULL)
summaryRprof("func.out",memory="both")

 Rprof("func.out",memory.profiling=T)
 system.time(p <- bbd_prob(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=0,B=50))
 Rprof(NULL)
summaryRprof("func.out",memory="both")

 
 source("bbd_prob0.r")
 source("bbd_lt0.r")
 source("contfr0.r")
 
 fun <-function(x,y) {return(bbd_prob0(t=1,0,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},x,y,B=50))}
 grid = expand.grid(0,0:50)
 system.time(p0 <- mapply(fun,grid[,1],grid[,2]))

 
  	system.time(phi<-bbd_phi(s=1,3,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=20,B=20))

	system.time(phi0<-bbd_phi0(s=1,3,3,lambda1=function(a,b){0},lambda2=function(a,b){return(0.5*b)},mu2=function(a,b){return(0.3*b)},gamma=function(a,b){0},A=20,B=20))
	
j = 5
l = 0
h =	10
	count = start_count(j,l,h)
	while (!(is.na(count[1]))) {
		print(count)
		count = add_count(count,j,l,h)
	}

  
  