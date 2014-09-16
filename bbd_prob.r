### Sep 10 2014
### Lam Ho
### Transition probability of a birth/birth-death process

bbd_prob <- function(t,a0,b0,brates1,brates2,drates2,A,B) {
	br1 <- function(a,b){return(brates1(a,b))}
	br2 <- function(a,b){return(brates2(a,b))}
	dr2 <- function(a,b){return(drates2(a,b))}
	res = lapply(t,bbd_lt_invert,f=function(s){return(bbd_lt(s,a0,b0,br1,br2,dr2,A,B))})
	if(any(is.na(res))) cat("bbd_prob(",a0,",",b0,",",t,") failed\n")		
	return(res)
}
