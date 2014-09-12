### Sep 10 2014
### Lam Ho
### Laplace transform of a birth/birth-death process

bbd_lt <- function(s,a0,b0,brates1,brates2,drates2,A,B,a,b) {
	# initialize phi
	phi = bbd_phi(s,a0,b0,brates1,brates2,drates2,A,B)
	# phi[a,b,m]
	size = length(s)
	f = array(0,dim=c(size,A+1,B+1))
	for (j in 0:B) {
		f[,a0+1,j+1] = phi[,a0+1,j+1,b0+1]
	}
	for (i in (a0+1):A) {	
		br1 = rep(NA,(B+1))
		for (k in 0:B) br1[k+1] = brates1(k)
		for (j in 0:B) {
			v1 = as.vector(br1)
			v2 = as.vector(f[,i,])
			v3 = as.vector(phi[,i+1,j+1,])
			f[,i+1,j+1] = sum(v1*v2*v3)
		}
	}
	return(f[,a+1,b+1])
}

bbd_phi <-function(s,a0,b0,brates1,brates2,drates2,A,B) {
	if ((a0<0)||(b0<0)||(a0>A)||(b0>B)) return(0)
	size = length(s)
	phi = array(0,dim=c(size,A+1,B+1,B+1))
	for (a in a0:A)
		for (n in 0:B) {
			xf <- function(k){
					if (k==1) x = -1/drates2(a,1)
					else x = -brates2(a,k-2)/drates2(a,k) 	
					return(x)
				}
			yf <- function(k) {
					y = (s + brates1(a,k-1) + brates2(a,k-1) + drates2(a,k-1))/drates2(a,k)
					return(y)
				}		
			for (m in 0:B) {				
				if(n<=m) {
    				if(n==m) { fac = -1/drates2(a,n+1) } 
    				else { fac = (-1)^(m-n+1)/drates2(a,n+1)/prod(xf((n+2):(m+1)))}
    				B1 = cf_BidBj(n,m,xf,yf)
    				B2 = 1/cf_Bk1dBk(m+1,xf,yf)
    				v = fac * B1 / (B2 + cf_lentz_vec_m(m+1,xf,yf))
    				} else {
    					fac = -prod(xf((m+2):n))/drates2(a,n+1)
    					B1 = cf_BidBj(m,n,xf,yf)
    					B2 = 1/cf_Bk1dBk(n+1,xf,yf)
    					v = fac * B1 / (B2 + cf_lentz_vec_m(n+1,xf,yf)) 
  						}
  				phi[,a+1,n+1,m+1] = v				
				}
			}
	return(phi)
}

