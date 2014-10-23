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

cf_BidBj <- function(B,xvec,yvec,Bk1dBk) {
	res = matrix(0,nrow=B+1,ncol=B+1)
	ans = rep(1,(B+1))
	
    tmp = .C("cf_BidBj",as.integer(B),as.double(xvec),as.complex(yvec),
		as.complex(Bk1dBk),as.complex(res),as.complex(ans))
    res = matrix(tmp[[5]], nrow = B+1, byrow = T)
	return(res)
}

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

