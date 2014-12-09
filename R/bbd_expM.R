### Dec 9 2014
### Lam Ho
### Transition probability of a birth/birth-death process
### Method: matrix exponential

bbd_expM <- function(t,a0,b0,lambda1,lambda2,mu2,gamma,A,B) {
  
  Q = matrix(0,nrow = (A-a0+1)*(B+1), ncol = (A-a0+1)*(B+1))
  P0 = rep(0,(A-a0+1)*(B+1))
  
  for (i in a0:A) {
    for (j in 0:B) {
      tmp = (i-a0)*(B+1) + j+1
      Q[tmp,tmp] = - lambda1(i,j) - lambda2(i,j) - mu2(i,j) - gamma(i,j)
      if (i<A) {
        tmp1 = (i-a0+1)*(B+1) + j+1
        Q[tmp,tmp1] =  lambda1(i,j)
      }
      if (j<B) {
        tmp1 = (i-a0)*(B+1) + j+2
        Q[tmp,tmp1] =  lambda2(i,j)
      }
      if (j>0) {
        tmp1 = (i-a0)*(B+1) + j
        Q[tmp,tmp1] = mu2(i,j)
      }
      if ((i<A)&&(j>0)) {
        tmp2 = (i-a0+1)*(B+1) + j
        Q[tmp,tmp2] = gamma(i,j)
      }
    }
  }
  tmp0 = b0+1
  P0[tmp0] = 1
  
  expM <- expAtv(t(Q),P0,t)$eAtv
  P = matrix(0,nrow=A-a0+1, ncol=B+1)
  for (i in a0:A) {
    for (j in 0:B) {
      P[i-a0+1,j+1] = expM[(i-a0)*(B+1) + j+1]
    }
  }
  colnames(P) = 0:B
  rownames(P) = a0:A
  return(P)
}

dbd_expM <- function(t,a0,b0,mu1,lambda2,mu2,gamma,a,B) {
  
  Q = matrix(0,nrow = (a0-a+1)*(B+1), ncol = (a0-a+1)*(B+1))
  P0 = rep(0,(a0-a+1)*(B+1))
  
  for (i in a:a0) {
    for (j in 0:B) {
      tmp = (i-a)*(B+1) + j+1
      Q[tmp,tmp] = - mu1(i,j) - lambda2(i,j) - mu2(i,j) - gamma(i,j)
      if (i>A) {
        tmp1 = (i-a-1)*(B+1) + j+1
        Q[tmp,tmp1] =  mu1(i,j)
      }
      if (j<B) {
        tmp1 = (i-a)*(B+1) + j+2
        Q[tmp,tmp1] =  lambda2(i,j)
      }
      if (j>0) {
        tmp1 = (i-a)*(B+1) + j
        Q[tmp,tmp1] =  mu2(i,j)
      }
      if ((i>a)&&(j<B)) {
        tmp2 = (i-a-1)*(B+1) + j+2
        Q[tmp,tmp2] =  gamma(i,j)
      }
    }
  }
  tmp0 = (a0-a)*(B+1) + b0+1
  P0[tmp0] = 1
  
  expM <- expAtv(t(Q),P0,t)$eAtv
  P = matrix(0,nrow=a0-a+1, ncol=B+1)
  for (i in a:a0) {
    for (j in 0:B) {
      P[i-a+1,j+1] = expM[(i-a)*(B+1) + j+1]
    }
  }
  colnames(P) = 0:B
  rownames(P) = a:a0
  return(P)
}