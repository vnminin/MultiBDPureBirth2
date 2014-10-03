### Sep 10 2014
### Lam Ho
### Laplace transform of a birth/birth-death process

bbd_lt <- function(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B) {
	# initialize phi
	phi = bbd_phi(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B)
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
				#v1 = as.vector(br1)
				v1 = br1
				#v2 = as.vector(f[i,])
				v2 = f[i,]
				#v3 = as.vector(phi[i+1,j+1,])
				v3 = phi[i+1,j+1,]
				#v4 = as.vector(g)
				v4 = g
				#v5 = as.vector(c(f[i,2:(B+1)],0))
				v5 = c(f[i,2:(B+1)],0)
				f[i+1,j+1] = sum((v1*v2+v4*v5)*v3)
			}
		}
	}
	return(f)
}

bbd_phi <- function(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B) {
	if ((a0<0)||(b0<0)||(a0>A)||(b0>B)) return(0)
	phi = array(0,dim=c(A+1,B+1,B+1))
	# phi[a,b,m]
	for (a in a0:A) {
		xvec = x[a-a0+1,]
		yvec = s+y[a-a0+1,]
		lentz = sapply(1:(B+1),cf_lentz_m,xvec,yvec)
		Bk1dBk = cf_Bk1dBk(B,xvec,yvec)
		BidBj = cf_BidBj(B,xvec,yvec,Bk1dBk)
		prod_mu2 = prod_vec(a-a0+1,B,mu2)
		prod_lambda2 = prod_vec(a-a0+1,B,lambda2)
		#for (b in 0:B) {
		#	for (m in 0:B) {				
		#			if(b<=m) {						
    	#				if (b==m) fac = 1 else {
    	#					fac = prod_mu2[b+2,m+1]
    	#				}						
    	#				B1 = BidBj[b+1,m+1]
    	#				B2 = 1/Bk1dBk[m+1]
    	#				v = fac * B1 / (B2 + lentz[m+1])
    	#				} else {
    	#					fac = prod_lambda2[m+1,b]
    	#					B1 = BidBj[m+1,b+1]
      	#					B2 = 1/Bk1dBk[b+1]
    	#					v = fac * B1 / (B2 + lentz[b+1]) 
  		#					}
  		#			phi[a+1,b+1,m+1] = v				
		#		}
		#	}
		tmp = .C("phi_routine", as.integer(B), as.double(as.vector(prod_mu2)), as.double(as.vector(prod_lambda2)), as.complex(Bk1dBk), as.complex(as.vector(BidBj)), as.complex(lentz), as.complex(as.vector(phi[a+1,,])))
		phi[a+1,,] = matrix(tmp[[7]], nrow=B+1, byrow=T)
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
