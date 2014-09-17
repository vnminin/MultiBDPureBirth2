### Sep 10 2014
### Lam Ho
### Transition probability of a birth/birth-death process

bbd_prob <- function(t,a0,b0,lambda1,lambda2,mu2,gamma,A,B) {
	l1 <- function(a,b){return(lambda1(a,b))}
	l2 <- function(a,b){return(lambda2(a,b))}
	m2 <- function(a,b){return(mu2(a,b))}
	g <- function(a,b){return(gamma(a,b))}
	res = lapply(t,bbd_lt_invert,f=function(s){return(bbd_lt(s,a0,b0,l1,l2,m2,g,A,B))})
	if(any(is.na(res))) cat("bbd_prob(",a0,",",b0,",",t,") failed\n")		
	return(res)
}
