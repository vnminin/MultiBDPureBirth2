### Sep 10 2014
### Lam Ho
### Transition probability of a birth/birth-death process


bbd_prob <- function(t,a0,b0,lambda1,lambda2,mu2,gamma,A,B,
                     nblocks=256,tol=1e-12,computeMode=0,nThreads=4,
                     doJIT=TRUE,maxdepth=400) {
	if (doJIT) enableJIT(1)
  
  ## R-C interface
#   dyn.load("src/cf_BidBj.so")
#   dyn.load("src/prod_vec.so")
#   dyn.load("src/phi_routine.so")

  if (a0<0) stop("a0 cannot be negative.")
  if (a0>A) stop("a0 cannot be bigger than A.")
  if (B<0) stop("B cannot be negative.")
  if (t<0) stop("t cannot be negative.")
  if (t<tol) {
    res = matrix(0,nrow=A-a0+1,ncol=B+1)
    res[1,b0+1] = 1
    colnames(res) = 0:B
    rownames(res) = a0:A
    return(res)
  }
	
	grid  = expand.grid(a0:A,0:(B+maxdepth))
	l1 = matrix(mapply(lambda1,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	l2 = matrix(mapply(lambda2,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	m2 = matrix(mapply(mu2,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	g = matrix(mapply(gamma,grid[,1],grid[,2]),ncol=B+1+maxdepth)
	
	xf <- function(a,b){
		if (b==0) x = 1
			else x = - l2[a-a0+1,b]*m2[a-a0+1,b+1] 	
		return(x)
	}
	x = matrix(mapply(xf,grid[,1],grid[,2]),nrow=B+1+maxdepth,byrow=TRUE)
	
	yf <- function(a,b) {
		y = l1[a-a0+1,b+1] + l2[a-a0+1,b+1] + m2[a-a0+1,b+1] + g[a-a0+1,b+1]
		return(y)
	}
	y = matrix(mapply(yf,grid[,1],grid[,2]),nrow=B+1+maxdepth,byrow=TRUE)	
					
# 	res = bbd_lt_invert(t,f=function(s) {
#     ## R
# 		#return(bbd_lt(s,a0,b0,l1,l2,m2,g,x,y,A,B))
#     
#     ## Rcpp and R
# 		return(matrix(bbd_lt_Cpp(s,a0,b0,l1,l2,m2,g,x,y,A,B),nrow=(A-a0+1),byrow=T))
# 		})

  ## Rcpp
  res = matrix(bbd_lt_invert_Cpp(t,a0,b0,l1,l2,m2,g,x,y,A,B+1,
                                 nblocks,tol,computeMode,nThreads,maxdepth),
               nrow=(A-a0+1),byrow=T)
	  
  #if(any(is.na(res))) cat("bbd_prob(",a0,",",b0,",",t,") failed\n")
  colnames(res) = 0:B
  rownames(res) = a0:A
	
	if (doJIT) enableJIT(0)
	return(abs(res))
}

dbd_prob <-function(t,a0,b0,mu1,lambda2,mu2,gamma,a,B,
                    nblocks=256,tol=1e-12,computeMode=0,nThreads=4,
                    doJIT=TRUE,maxdepth=400) {
  ## a>=0, a<=a0, B >=a0+b0-a 
  if(a<0) stop("a cannot be negative.")
  if (a>a0) stop("a0 canot be smaller than a0.")
  if (B < a0+b0-a) stop("B is too small.")
  
	l1 <- function(u,v){
    if (v>B) return(0)
    return(mu1(a0-u,B-v))
	}
	l2 <- function(u,v){
	  if (v>B) return(0)
    return(mu2(a0-u,B-v))
	}
	m2 <- function(u,v){
	  if (v>B) return(0)
    return(lambda2(a0-u,B-v))
	}
	g <- function(u,v){
	  if (v>B) return(0)
    return(gamma(a0-u,B-v))
	}
	res = matrix(0,nrow=a0-a+1,ncol=B+1)
	res[(a0-a+1):1,(B+1):1] = bbd_prob(t,0,B-b0,l1,l2,m2,g,A=a0-a,B,
	                                   nblocks,tol,computeMode,nThreads,
                                     doJIT,maxdepth)

  colnames(res) = 0:B
  rownames(res) = a:a0
	return(res)
}
