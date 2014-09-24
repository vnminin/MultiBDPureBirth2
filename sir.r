### Sep 23 2014
### Lam Ho
### Transition probabilities for SIR (Maroufy et al 2012)

sir <- function(t,n,a,f,mu){
	p = matrix(0,nrow=n+1,ncol=n+a+1)
	for (l in 0:a) {
		if (l<a) tmp = (a-l)*log(mu)+sum(log((l+1):a)) else tmp = 0
		spr = 0
		for (k in l:a){
			pr = 1
			for (kp in l:a){				
				if (!(kp==k)) pr = pr/(mu*(kp-k)+f(n,kp)-f(n,k))
			}
			pr = pr/exp((mu*k+f(n,k))*t)
			spr = spr+pr
		}
		tmp = tmp + log(spr)		
		p[n+1,l+1] = exp(tmp)
	}
	
	for (i in 0:(n-1))
		for (l in 0:(n+a-i)) {
			if (l<a) tmp = (n+a-i-l)*log(mu)+sum(log((l+1):a)) else tmp = (n-i)*log(mu)
			count = start_count(n-i,l,n+a-i)
			spr = 0			
			while (!(is.na(count[1]))) {
				lc = rep(NA,n-i+2)
				lc[1] = l
				lc[n-i+2] = a+i
				lc[2:(n-i+1)] = count				
				for (w in 0:(n-i)) {
					for (k in 0:(lc[w+2]-lc[w+1])) {
						pr = 1
						ff <- function(i,w,k) {
							return(mu*(lc[w+1]-w+k)+f(i+w,lc[w+1]-w+k))
						}
						for (wp in 0:(n-i))
							for (kp in 0:(lc[w+2]-lc[w+1])) {
								if ((abs(w-wp)>0)&&(abs(k-kp)>0)) 
									pr = pr/(ff(i,wp,kp)-ff(i,w,k))
							}
						pr = pr/exp(ff(i,w,k)*t)
						for (j in 1:(n-i)) {
							pr = pr*(lc[j+1]-j+1)*f(i+j,lc[j+1]-j)
						}
						spr = spr+pr
					}
				}
				count = add_count(count,n-i,l,n+a-i)
			}
			print(spr)
			tmp = tmp + log(spr)
			p[i+1,l+1] = exp(tmp)
		}
	return(p)
}

start_count <- function(j,l,h) {
	count = rep(1,j)
	count[1] = max(l,2)
	if (j>=2)
		for (i in 2:j)
			while ((count[i]<i)&&(count[i]<count[i-1])) count[i] = count[i]+1
	if (count[j]>h) return(NA)		
	return(count)
}

add_count <- function(count,j,l,h) {
	if (j>=2)
		for (i in 1:(j-1)) 
			if (count[i] < count[i+1]) {
				count[i] = count[i]+1
				return(count)
			}
	if (count[j] < h) {
		count[j] = count[j]+1
		return(count)
	} else return(NA)
}