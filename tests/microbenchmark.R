library(microbenchmark)

################
### Preparing data
a0 = 0
b0 = 0
A = 50
B = A

muL = runif(1,0,1)
muM = 0.0015
eta = runif(1,0,1)
gamma = 0.04

brates1=function(a,b){muL*a+eta*a^2}
drates2=function(a,b){0}
brates2=function(a,b){muM*b}
trans=function(a,b){gamma*a} # a -> b

maxdepth = 400
grid  = expand.grid(a0:A,0:(B+maxdepth))
l1 = matrix(mapply(brates1,grid[,1],grid[,2]),ncol=B+1+maxdepth)
l2 = matrix(mapply(brates2,grid[,1],grid[,2]),ncol=B+1+maxdepth)
m2 = matrix(mapply(drates2,grid[,1],grid[,2]),ncol=B+1+maxdepth)
g = matrix(mapply(trans,grid[,1],grid[,2]),ncol=B+1+maxdepth)

xf <- function(a,b){
  if (b==0) x = 1
  else x = - l2[a-a0+1,b]*m2[a-a0+1,b+1] 	
  return(x)
}
x = matrix(mapply(xf,grid[,1],grid[,2]),ncol=B+1+maxdepth)

yf <- function(a,b) {
  y = l1[a-a0+1,b+1] + l2[a-a0+1,b+1] + m2[a-a0+1,b+1] + g[a-a0+1,b+1]
  return(y)
}
y = matrix(mapply(yf,grid[,1],grid[,2]),ncol=B+1+maxdepth)	
################
microbenchmark(phi_Cpp(complex(real=1, imaginary=1),a0,b0,l2,m2,x,y,A,B),times = 2000)

###### Results from benchmark

### Removing code duplication in phi_Cpp slows down the function.
## no code duplication: 
# min       lq     mean   median      uq      max neval
# 4.903648 5.036701 5.610665 5.212402 6.00686 39.28769  2000
## code duplication:
# min      lq     mean   median       uq      max neval 
# 4.21585 4.24143 4.759771 4.443558 5.073127 37.68972  2000

### reciprocal(complex) and one/complex are the same
## one/complex:
# min       lq     mean  median       uq      max neval
# 4.230024 4.273974 4.777658 4.46344 5.088186 36.86331  2000
## reciprocal:
# min       lq     mean   median       uq      max neval
# 4.211236 4.241655 4.829227 4.476201 5.128305 37.92462  2000
