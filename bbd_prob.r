### Sep 10 2014
### Lam Ho
### Transition probability of a birth/birth-death process


bbd_prob <- function(t,a0,b0,lambda1,lambda2,mu2,gamma,A,B) {
	maxdepth = 400
	enableJIT(1)
	dyn.load("cf_BidBj.so")
	dyn.load("prod_vec.so")
	dyn.load("phi_routine.so")
	#l1 <- function(a,b){return(lambda1(a,b))}
	#l2 <- function(a,b){return(lambda2(a,b))}
	#m2 <- function(a,b){return(mu2(a,b))}
	#g <- function(a,b){return(gamma(a,b))}
	
	grid  = expand.grid(a0:A,0:(B+maxdepth))
	l1 = matrix(mapply(lambda1,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	l2 = matrix(mapply(lambda2,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	m2 = matrix(mapply(mu2,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	g = matrix(mapply(gamma,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	
	xf <- function(a,b){
		if (b==0) x = 1
			else x = - l2[a-a0+1,b]*m2[a-a0+1,b+1] 	
		return(x)
	}
	x = matrix(mapply(xf,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	
	yf <- function(a,b) {
		y = l1[a-a0+1,b+1] + l2[a-a0+1,b+1] + m2[a-a0+1,b+1] + g[a-a0+1,b+1]
		return(y)
	}
	y = matrix(mapply(yf,grid[,1],grid[,2]),ncol=B+1+maxdepth)	
	
	#pmu2 = array(0,dim=c(A-a0+1,B+1,B+1))
	#pl2 = array(0,dim=c(A-a0+1,B+1,B+1))
	#for (a in a0:A) {
	#	pmu2[a-a0+1,,] = prod_vec(a-a0+1,B,m2)
	#	pl2[a-a0+1,,] = prod_vec(a-a0+1,B,l2)
	#}		
				
	res = bbd_lt_invert(t,f=function(s) {
		return(bbd_lt(s,a0,b0,l1,l2,m2,g,x,y,A,B))
		#return(bbd_lt(s,a0,b0,l1,l2,m2,g,x,y,pmu2,pl2,A,B))
		})
	#if(any(is.na(res))) cat("bbd_prob(",a0,",",b0,",",t,") failed\n")	
	
	enableJIT(0)
	return(res)
}

dbd_prob <-function(t,a0,b0,mu1,lambda2,mu2,gamma,B) {
	l1 <- function(a,b){return(mu1(a0-a,B-b))}
	l2 <- function(a,b){return(mu2(a0-a,B-b))}
	m2 <- function(a,b){return(lambda2(a0-a,B-b))}
	g <- function(a,b){return(gamma(a0-a,B-b))}
	res = matrix(NA,nrow=a0+1,ncol=B+1)
	res[(a0+1):1,(B+1):1] = bbd_prob(t,0,B-b0,l1,l2,m2,g,A=a0,B)
	return(res)
}
