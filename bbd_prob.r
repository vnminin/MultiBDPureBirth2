### Sep 10 2014
### Lam Ho
### Transition probability of a birth/birth-death process

bbd_prob <- function(t,a0,b0,brates1,brates2,drates2,A,B,a,b) {
	# use lt_invert function in bd_lt.r
	res = sapply(t,lt_invert,f=function(s){return(bbd_lt(s,a0,b0,brates1,brates2,drates2,A,B,a,b))})
	if(any(is.na(res))) cat("bbd_prob(",a0,",",b0,",",t,") failed\n")		
	return(res)
}
