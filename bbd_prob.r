### Sep 10 2014
### Lam Ho
### Transition probability of a birth/birth-death process

bbd_prob <- function(t,a0,b0,a,b,brates1,brates2,drates2,A,B) {
	fhat = function(s){
		return(bbd_lt(s,a0,b0,brates1,brates2,drates2,A,B))
		}
	res = sapply(t,lt_invert,f=fhat)
	if(any(is.na(res))) cat("bbd_prob(",a0,",",b0,",",a,",",b,",",t,") failed\n")
	return(res)
}
