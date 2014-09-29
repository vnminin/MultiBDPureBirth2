### Sep 10 2014
### Lam Ho
### Laplace transform of a birth/birth-death process

bbd_lt <- function(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B) {
	# initialize phi
	phi = bbd_phi(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B)
	# phi[a,b,m]
	f = matrix(0,nrow=A+1,ncol=B+1)
	# f[a,b]
	f[a0+1,] = phi[a0+1,,b0+1]
	if (a0 < A) {
		for (i in (a0+1):A) {	
			#br1 = sapply(0:B,lambda1,a=i-1)
			br1 = lambda1[i-a0,1:(B+1)]
			#g = c(sapply(1:B,gamma,a=i-1),0)
			g = c(gamma[i-a0,2:(B+1)],0)
			for (j in 0:B) {
				v1 = as.vector(br1)
				v2 = as.vector(f[i,])
				v3 = as.vector(phi[i+1,j+1,])
				v4 = as.vector(g)
				v5 = as.vector(c(f[i,2:(B+1)],0))
				f[i+1,j+1] = sum((v1*v2+v4*v5)*v3)
			}
		}
	}
	return(f)
}

bbd_phi <- function(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B) {
	if ((a0<0)||(b0<0)||(a0>A)||(b0>B)) return(0)
	phi = array(0,dim=c(A+1,B+1,B+1))
	# phi[a,b,m]
	for (a in a0:A) {
		xf <- function(k){
					if (k<1) stop("error at xf")
					if (k==1) x = 1
						#else x = - lambda2(a,k-2)*mu2(a,k-1) 	
						else x = - lambda2[a-a0+1,k-1]*mu2[a-a0+1,k] 	
					return(x)
				}
		yf <- function(k) {
					if (k<1) stop("error at yf")
					#y = s + lambda1(a,k-1) + lambda2(a,k-1) + mu2(a,k-1) + gamma(a,k-1)
					y = s + lambda1[a-a0+1,k] + lambda2[a-a0+1,k] + mu2[a-a0+1,k] + gamma[a-a0+1,k]
					return(y)
				}
		lentz = sapply(1:(B+1),cf_lentz_m,xf,yf)
		for (b in 0:B) {
			for (m in 0:B) {				
					if(b<=m) {						
    					#if (b==m) fac = 1 else fac = prod(sapply((b+1):m,mu2,a=a))
    					if (b==m) fac = 1 else fac = prod(as.vector(mu2[a-a0+1,(b+2):(m+1)]))
    					B1 = cf_BidBj(b,m,xf,yf)
    					B2 = 1/cf_Bk1dBk(m+1,xf,yf)
    					v = fac * B1 / (B2 + lentz[m+1])
    					} else {
    						#fac = prod(sapply(m:(b-1),lambda2,a=a))	
    						fac = prod(as.vector(lambda2[a-a0+1,(m+1):b]))	
    						B1 = cf_BidBj(m,b,xf,yf)
    						B2 = 1/cf_Bk1dBk(b+1,xf,yf)
    						v = fac * B1 / (B2 + lentz[b+1]) 
  							}
  					phi[a+1,b+1,m+1] = v				
				}
			}
		}
	return(phi)
}

bbd_lt_invert = function(f,t,A=20) {
	cores = 5
	kmax = cores
  	ig = lapply(complex(real=A, imaginary=2*pi*(1:kmax))/(2*t),f) 	
  	#ig = mclapply(complex(real=A, imaginary=2*pi*(1:kmax))/(2*t),f,mc.cores = kmax)  	
  	psum0 = Re(f(A/(2*t))) / (2*t)
  	nr = nrow(ig[[1]])
  	nc = ncol(ig[[1]])
  	result = matrix(NA,nrow=nr,ncol=nc)
  	tol = 1e-12
  		
  	for (i in 1:nr)
  		for (j in 1:nc) {			
  			levin_init(20,tol) # initialize the levin acceleration method.
  			term = 10e30
  			sdiff = 10e30
  			k = 1
  			psum = psum0[i,j]
  			while((abs(sdiff) > 1e-16)||(abs(term) > 1e-3)) {
    			term = (-1)^k * Re(ig[[k]][i,j]) / t
    			psum = psum + term
    			omega = k*term
				sk = next_approx(psum,omega)
	    		if (is.na(sk)) return(NA)
	    		if (k>1) {
    				sdiff = sk - sk1
    				} else {sdiff = 10e30}
    			k = k+1	
    			sk1 = sk
	    		if (k>kmax) {
      			ig[(kmax+1):(kmax+cores)] = lapply(complex(real=A,imaginary=2*pi*((kmax+1):(kmax+cores)))/(2*t),f)
      			#ig[(kmax+1):(kmax+cores)] = mclapply(complex(real=A,imaginary=2*pi*((kmax+1):(kmax+cores)))/(2*t),f,mc.cores = cores)
	    		kmax = kmax + cores
	    		#print(c(i,j,kmax))
    			}   		
  			}
		result[i,j] = sk1*exp(A/2)
  		} 
  	return(result)
}
