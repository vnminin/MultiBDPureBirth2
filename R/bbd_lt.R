### Sep 10 2014
### Lam Ho
### Laplace transform of a birth/birth-death process

bbd_lt <- function(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B) {
	# initialize phi
	#phi = bbd_phi(s,a0,b0,lambda2,mu2,x,y,A,B) 
  phi = array(phi_Cpp(2,s,a0,b0,lambda2,mu2,x,y,A,B),dim=c(A-a0+1,B+1,B+1))
	# phi[a,b,m]
	f = matrix(0,nrow=A+1-a0,ncol=B+1)
	# f[a,b]
	f[1,] = phi[1,,b0+1]
	if (a0 < A) {
		for (i in (a0+1):A) {
			br1 = lambda1[i-a0,1:(B+1)]
			g = c(gamma[i-a0,2:(B+1)],0)
			for (j in 0:B) {
				v1 = br1
				v2 = f[i-a0,]
				v3 = phi[i-a0+1,j+1,]
				v4 = g
				v5 = c(f[i-a0,2:(B+1)],0)
				f[i+1-a0,j+1] = sum((v1*v2+v4*v5)*v3)
			  #f[i+1-a0,j+1] =sum((lambda1[i-a0,1:(B+1)]*f[i,]+c(gamma[i-a0,2:(B+1)]*f[i,2:(B+1)],0))*phi[i+1,j+1,])
			}
		}
	}
	return(f)
}

bbd_phi <- function(s,a0,b0,lambda2,mu2,x,y,A,B) {
  # if ((a0<0)||(b0<0)||(a0>A)||(b0>B)) return(0)
	phi = array(0,dim=c(A+1-a0,B+1,B+1))
	# phi[a,b,m]
	for (a in a0:A) {
		yvec = s+y[a-a0+1,]
		#lentz = sapply(1:(B+1),cf_lentz_m,xvec = x[a-a0+1,],yvec)
		lentz = lentz_Cpp(B,xvec = x[a-a0+1,],yvec)
		#Bk1dBk = cf_Bk1dBk(B,xvec,yvec)
		Bk1dBk = Bk1dBk_Cpp(B,xvec = x[a-a0+1,],yvec)
		#BidBj = cf_BidBj(B,xvec,yvec,Bk1dBk)
		BidBj = matrix(BidBj_Cpp(B,xvec = x[a-a0+1,],yvec,Bk1dBk),nrow=(B+1),byrow=T)
    #prod_mu2 = prod_vec(a-a0+1,B,mu2)
    #prod_lambda2 = prod_vec(a-a0+1,B,lambda2)
    prod_mu2 = matrix(prod_vec_Cpp(a-a0+1,B,mu2),nrow=(B+1),byrow=T)
    prod_lambda2 = matrix(prod_vec_Cpp(a-a0+1,B,lambda2),nrow=(B+1),byrow=T)
		#tmp = .C("phi_routine", as.integer(B), as.double(prod_mu2), as.double(prod_lambda2), as.complex(Bk1dBk), as.complex(BidBj), as.complex(lentz), as.complex(phi[a+1,,]))
		#phi[a+1-a0,,] = matrix(tmp[[7]], nrow=B+1, byrow=T)
		phi[a+1-a0,,] = matrix(phi_routine_Cpp(B,prod_mu2,prod_lambda2,Bk1dBk,BidBj,lentz), nrow=B+1, byrow=T)
		}
	return(phi)
}

bbd_lt_invert = function(f,t,A=20) {
	cores = 3
	kmax = cores
  	ig = lapply(complex(real=A, imaginary=2*pi*(1:kmax))/(2*t),f) 	
  	#ig = mclapply(complex(real=A, imaginary=2*pi*(1:kmax))/(2*t),f,mc.cores = cores)  	
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
    			}   		
  			}
		result[i,j] = sk1*exp(A/2)
  		} 
  	return(result)
}
