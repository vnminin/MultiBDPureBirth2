### Sep 12 2014
### Lam Ho
### Computation with continued fraction
### Adapted from Forrest's code

cf_lentz_m <- function(m,xvec,yvec,maxdepth=400) {
	tiny = 1e-30
	eps = 1e-8
	j = m
	
	fj1 = tiny
	Cj = 0
	Cj1 = tiny
	Dj = 0
	Dj1 = 0
	jdiff = 2
	jbound = 1
	
	while(any(jbound > eps)) {
		aj = xvec[j+1]
		bj = yvec[j+1]
    	Dj = bj + aj*Dj1 
    	if (Dj == 0) Dj = tiny 
    	Cj = bj + aj/Cj1
    	if (Cj == 0) Cj = tiny
    	Dj = 1/Dj 
    	jdiff = Cj*Dj 
    	fj = fj1*jdiff
    	
    	truncerr = abs(fj-fj1)
    	if (truncerr==0) truncerr = tiny
    	jbound = ( abs(1/Dj)/abs(Im(1/Dj)) ) * truncerr
    	if (Im(Dj)==0) jbound = abs(jdiff-1)

	    j = j + 1	
    	fj1 = fj 
    	Dj1 = Dj
    	Cj1 = Cj

    	if((j-m) > maxdepth) {
      		cat("\n\n***************************\nlentz: maxdepth =",maxdepth," reached.\n")
      		cat("m =", m, "\n")
      		cat("fj =", fj, "\n")
      		return(NA)
    	}
  	}
  	return(fj)
}

#cf_BidBj <- function(i,j,xvec,yvec,Bk1dBk) {
#	if(i==j) return(1)
#	a = min(i,j)
#	b = max(i,j)
#	if(b==(a+1)) return(Bk1dBk[b])
#	
#	ans = rep(0,b+1)
#	ans[1] = 1 # Ba/Ba
#	ans[2] = 1/Bk1dBk[a+1] # Ba+1/Ba
#	
#	idx = 3
#	k = a + 2
#	while(k<=b) {
#		ak = xvec[k]
#	   	bk = yvec[k]
#   	ans[idx] = bk*ans[idx-1] + ak*ans[idx-2]
#
#   	idx = idx + 1
#  		k   = k + 1
#  	}
# 
#  	return(1/ans[idx-1])
#}

cf_BidBj <- function(B,xvec,yvec) {
	res = matrix(NA,nrow=B+1,ncol=B+1)
	Bk1dBk = cf_Bk1dBk(B,xvec,yvec)
	for (i in 0:B) {
		ans = rep(1,(B+1))
		for (j in i:B) {
			if (j==i) ans[j+1] = 1	
			else {
				if (j==(i+1)) ans[j+1] = 1/Bk1dBk[j]
				else ans[j+1] = yvec[j]*ans[j] + xvec[j]*ans[j-1]
				}
			res[i+1,j+1] = 1/ans[j+1]	
		}					
	}
	return(res)			
}

#cf_Bk1dBk <- function(k,xvec,yvec) {
#	if(k==0) {cat("k=0 not allowed!\n"); return(NA)}
#  
#  	tiny = 1e-30
#  	j = 0  
#  	Dj = 0
#  	Dj1 = 0
#	
#	while(j<k) {
#    	aj = xvec[j+1]
#    	bj = yvec[j+1]
#
#    	Dj = bj + aj*Dj1 
#    	if (Dj==0) Dj = tiny 
#    	Dj1 = 1/Dj 
#    	j = j + 1
#	}	
#  return(Dj1)
#}

cf_Bk1dBk <- function(B,xvec,yvec) {
	res = rep(NA,B+1)
  	tiny = 1e-30
   	Dj = 0
  	Dj1 = 0 
  	for (j in 0:B) {
  		aj = xvec[j+1]
  		bj = yvec[j+1]
  		Dj = bj + aj*Dj1
  		if (Dj==0) Dj = tiny 
  		Dj1 = 1/Dj
  		res[j+1] = Dj1
  	}
  	return(res)
}

