### Sep 10 2014
### Lam Ho
### Laplace transform of a birth/birth-death process

bbd_lt <- function(s,a0,b0,brates1,brates2,drates2,A,B) {
	# initialize phi
	phi = bbd_phi(s,a0,b0,brates1,brates2,drates2,A,B)
	# phi[a,b,m]
	f = array(0,dim=A,B)
	for (j in 1:B) {
		f[1,j] = phi[1,j,b0]
	}
	for (i in 2:A) {
		for (j in 1:B) {
			v1 = as.vector(brates1[i-1,])
			v2 = as.vector(f[i-1,])
			v3 = as.vector(phi[i,j,])
			f[i,j] = sum(v1*v2*v3)
		}
	}
}

bbd_phi <-function(s,a0,b0,brates1,brates2,drates2,A,B) {
	
}