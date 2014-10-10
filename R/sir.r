### Sep 23 2014
### Lam Ho
### Transition probabilities for SIR (Maroufy et al 2012)

sir <- function(t,n,a,f,mu,i,l){
	if ((i>n)||(l>(n+a-i))) return(0)
	if (i==n) {
		if (l<a) tmp = (a-l)*log(mu)+sum(log((l+1):a)) else tmp = 0
		bigs = 0
		for (k in l:a){
			pr = 1
			for (kp in l:a){				
				if (abs(kp-k)>0) pr = pr/(mu*(kp-k)+f(n,kp)-f(n,k))
			}
			pr = pr/exp((mu*k+f(n,k))*t)
			bigs = bigs+pr
		}
		tmp = tmp + log(bigs)		
		return(exp(tmp))
	}
	
	if (l<a) tmp = (n+a-i-l)*log(mu)+sum(log((l+1):a)) 
	if (l==a) tmp = (n-i)*log(mu)
	if (l>a) tmp = (n+a-i-l)*log(mu)-sum(log((a+1):l))
	count = start_count(n-i,l,n+a-i)
	bigs = 0			
	while (!(is.na(count[1]))) {
		lw = rep(NA,n-i+2)
		lw[1] = l
		lw[n-i+2] = n+a-i
		lw[2:(n-i+1)] = count					
		ff <- function(i,w,k) {
			return(mu*(lw[w+1]-w+k)+f(i+w,lw[w+1]-w+k))
		}
		smalls = 0			
		for (w in 0:(n-i)) 
			for (k in 0:(lw[w+2]-lw[w+1])) {
				pr = 1
				pc = 1						
				for (wp in 0:(n-i))
					for (kp in 0:(lw[wp+2]-lw[wp+1])) {
						if ((abs(w-wp)+abs(k-kp))>0) {
							if (abs(ff(i,wp,kp)-ff(i,w,k))==0) 											print("nooo")
							pr = pr/(ff(i,wp,kp)-ff(i,w,k))
						}
					}
				pr = pr/exp(ff(i,w,k)*t)
				smalls = smalls + pr
			}	
		for (w in 1:(n-i)) 
			smalls = smalls*(lw[w+1]-w+1)*f(i+w,lw[w+1]-w)
		bigs = bigs + smalls
		count = add_count(count,n-i,l,n+a-i)	
	}
	tmp = tmp + log(bigs)
	return(exp(tmp))
}

start_count <- function(j,l,h) {
	count = rep(1,j)
	count[1] = max(l,2)
	if (j>=2)
		for (i in 2:j)
			while ((count[i]<=i)||(count[i]<count[i-1])) count[i] = count[i]+1
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