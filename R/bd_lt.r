# this file contains functions for evaluating and inverting laplace transforms of BD probabilities.

################################################################################
# lentz_vec_m evaluates the continued fraction corresponding to 
# the BD probability starting at state m using the modified 
# Lentz algorithm (Press et al).

lentz_vec_m = function(s0, m, brates, drates, maxdepth=400) {
  #cat("s0 =", s0, "\n")
  
  tiny = 1e-30
  eps = 1e-8

  j = m  

  k = length(s0)

  fj1 = rep(tiny,k) 
  Cj = rep(0,k) 
  Cj1 = rep(tiny,k) 
  Dj = rep(0,k) 
  Dj1 = rep(0,k) 
  jdiff = rep(2,k) 
  jbound = rep(1,k)

  while(any(jbound > eps)) {

    #cat("j =", j, "\n")

    if(j==0) {aj = 1}
    else { aj = -1*brates(j-1)*drates(j) }
    bj = s0 + brates(j) + drates(j)
    Dj = bj + aj*Dj1 
    #cat(" D[",j,"] = ", Dj, "\n", sep="")
    Dj[Dj == 0] = tiny 
    Cj = bj + aj/Cj1
    Cj[Cj == 0] = tiny
    Dj = 1/Dj 
    jdiff = Cj*Dj # vec element-wise mult
    fj = fj1*jdiff # vec element-wise mult

    truncerr = abs(fj-fj1)
    truncerr[truncerr==0] = tiny
    jbound = ( abs(1/Dj)/abs(Im(1/Dj)) ) * truncerr
    #cat("j =", j, "\n")
    #cat(" Dj =", Dj, "\n")
    #cat(" jbound =", jbound, "\n")
    #cat(" truncerr =", truncerr, "\n")
    #cat(" fj =", fj, "\n")
    #cat(" fj1 =", fj1, "\n")
    jbound[Im(Dj)==0] = abs(jdiff-1)
    

    j = j + 1
    fj1 = fj 
    Dj1 = Dj
    Cj1 = Cj

    if((j-m) > maxdepth) {
      cat("\n\n***************************\nlentz: maxdepth =",maxdepth," reached.\n")
      cat("m =", m, ", s0 =", s0, "\n")
      cat("fj =", fj, "\n")
      return(NA)
    }
  }
  return(fj)
}

###############################################
# always returns Ba/Bb = small/large

BidBj = function(s0, i, j, brates, drates) {
  #cat("BidBj(s0,",i,",",j,")\n")

  if(i==j) return(rep(1,length(s0)))

  a = min(i,j)
  b = max(i,j)
  if(b==(a+1)) return(Bk1dBk(s0,a+1,brates,drates)) 

  ans = array(0,dim=c(b+1,length(s0)))
  ans[1,] = rep(1,length(s0)) # Ba/Ba
  ans[2,] = 1/Bk1dBk(s0,a+1,brates,drates) # Ba+1/Ba

  idx = 3
  k   = a + 2
  while(k<=b) {
    ak = -1*brates(k-2)*drates(k-1)
    bk = s0 + brates(k-1) + drates(k-1)
    ans[idx,] = bk*ans[idx-1,] + ak*ans[idx-2,]

    idx = idx + 1
    k   = k + 1
  }

  return(1/ans[idx-1,])

}

#################################################
# returns B_{k-1}/B_k

Bk1dBk = function(s0, k, brates, drates) {
  #cat("Bk1Bk(",k,")\n")

  if(k==0) {cat("k=0 not allowed!\n"); return(NA)}
  
  tiny = 1e-30
  j = 0  
  Dj = rep(0,length(s0))
  Dj1 = rep(0,length(s0))

  while(j<k) {

    if(j==0) {aj = 1}
    else { aj = -1*brates(j-1)*drates(j) }
    bj = s0 + brates(j) + drates(j)

    Dj = bj + aj*Dj1 
    Dj[Dj == 0] = tiny 
    Dj1 = 1/Dj 
    j = j + 1

  }
  return(Dj1)
}



###############################################
# wallisB3_lentz returns matrix of partial denominators 

wallisB3_lentz = function(j1,j2,j3,s0,brates,drates) {
  n = length(s0)
  if((j1==0)&&(j2==0)&&(j3==0)) {return(rbind(rep(1,n),rep(1,n),rep(1,n)))}
  jmax = max(j1,j2,j3)
  B = array(0,dim=c(jmax+2,n))
  B[1,] = rep(1,n) # starting values
  B[2,] = s0+brates(0)

  for(i in 3:(jmax+2)) {
    B[i,] = (s0+brates(i-2)+drates(i-2))*B[i-1,] - brates(i-3)*drates(i-2)*B[i-2,]
  }

  return(rbind(B[j1+1,],B[j2+1,],B[j3+1,]))
}

####################################################
# lt_eval evaluates the laplace transform of the transition probability.

lt_eval_mod = function(m,n,s,brates,drates) {
  #cat("lt_eval(s0 =", s, ")\n")
  #num_lt_evals <<- num_lt_evals + 1

  #cat("  lt_eval(",m,",",n,",",s,")\n")
  if((m<0)||(n<0)) return(0)
  x = s 
  v = 0
  if(n<=m) {
    if(n==m) { fac = 1 } 
    else { fac = prod(drates((n+1):m)) }
    #wB = wallisB3_lentz(n,m,m+1,x,brates,drates)
    #v = fac * wB[1,] / (wB[3,] + wB[2,] * lentz_vec_m(x,m+1,brates,drates)) 
    B1 = BidBj(x,n,m,brates,drates)
    B2 = 1/Bk1dBk(x,m+1,brates,drates)
    #cat("B1 =", B1[1], "\n")
    #cat("B2 =", B2[1], "\n")
    v = fac * B1 / (B2 + lentz_vec_m(x,m+1,brates,drates)) 
    #cat("Bn/Bm:\n") 
    #cat("  wallis =",wB[1,]/wB[2,],"\n")
    #cat("  new    =",B1,"\n")
    #cat("Bm+1/Bm:\n") 
    #cat("  wallis =",wB[3,]/wB[2,],"\n")
    #cat("  new    =",B2,"\n")

  } else { # n > m
    fac = prod(brates(m:(n-1)))
    B1 = BidBj(x,m,n,brates,drates)
    B2 = 1/Bk1dBk(x,n+1,brates,drates)

    #cat("B1 =", B1, "\n")
    #cat("B2 =", B2, "\n")
    #cat("fac =", fac, "\n")

    
    v = fac * B1 / (B2 + lentz_vec_m(x,n+1,brates,drates)) 

  }
  return(v)
}


####################################################
# lt_eval evaluates the laplace transform of the transition probability.

#lt_eval = function(m,n,s,brates,drates) {
  #cat("lt_eval(s0 =", s, ")\n")
  #num_lt_evals <<- num_lt_evals + 1

  #cat("  lt_eval(",m,",",n,",",s,")\n")
  #if((m<0)||(n<0)) return(0)
  #x = s 
  #v = 0
  #if(n<=m) {
    #if(n==m) { fac = 1 } 
    #else { fac = prod(drates((n+1):m)) }
    #wB = wallisB3_lentz(n,m,m+1,x,brates,drates)
    #v = fac * wB[1,] / (wB[3,] + wB[2,] * lentz_vec_m(x,m+1,brates,drates)) 
    #B1 = BidBj(x,n,m,brates,drates)
    #B2 = 1/Bk1dBk(x,m+1,brates,drates)
    #v = fac * B1 / (B2 + lentz_vec_m(x,m+1,brates,drates)) 

  #} else { # n > m
    #fac = prod(brates(m:(n-1)))
    #wB = wallisB3_lentz(m,n,n+1,x,brates,drates)
    ##numerator = fac*wB[1,]
    ##denominator = wB[3,] + wB[2,] * lentz_vec_m(x,n+1,brates,drates) 
    ##cat("num =", numerator, ", denom =", denominator,"\n")
    ##cat("lentz =",lentz_vec_m(x,n+1,brates,drates),"\n")
    #cat("wB3 =",wB[3,],"\n")
    #cat("wB2 =",wB[2,],"\n")
    #v = fac * wB[1,] / (wB[3,] + wB[2,] * lentz_vec_m(x,n+1,brates,drates) )
    #v = numerator / denominator
    #B1 = BidBj(x,m,n,brates,drates)
    #B2 = 1/Bk1dBk(x,n+1,brates,drates)

    #cat("Bm/Bn:\n") 
    #cat("  wallis =",wB[1,]/wB[2,],"\n")
    #cat("  new    =",B1,"\n")

    #cat("Bn+1/Bn:\n") 
    #cat("  wallis =",wB[3,]/wB[2,],"\n")
    #cat("  new    =",B2,"\n")

    #v = fac * B1 / (B2 + lentz_vec_m(x,n+1,brates,drates)) 

  #}
  #return(v)
#}

###########################################################
# inverts the laplace transform f using the method popularized by abate and whitt

lt_invert = function(f,t,A=20) {
  #cat("lt_invert(f,",t,"\n")
  kmax = 5
  ig = f(complex(real=A, imaginary=2*pi*(1:kmax))/(2*t))
  #ig = sapply(complex(real=A, imaginary=2*pi*(1:kmax))/(2*t),f)

      #vals = complex(real=A,im=2*pi*(1:kmax))/(2*t)
  #cat("ig =", ig, "\n")
  psum = Re(f(A/(2*t))) / (2*t)
  s = rep(0,40)
  tol = 1e-12
  levin_init(20,tol) # initialize the levin acceleration method.
  term = 10e30
  sdiff = 10e30
  k = 1
  while((abs(sdiff) > 1e-16)||(abs(term) > 1e-3)) {
 
    #cat("k =", k, "\n")

    term = (-1)^k * Re(ig[k]) / t
    #cat("k =",k,", ig[k] =", ig[k], "term =",term,"\n")

    psum = psum + term
    omega = k*term

    # omega = k*term = 0 is causing next_approx to fail
    # when term =  0.

    #cat("psum =",psum,", omega =",omega,"\n")
    s[k] = next_approx(psum,omega)


    #cat("term =", term, "\n")
    #cat("psum =", psum, "\n")
    #cat("k =", k, "term =", term, "p_approx =", s[k], "\n")

    if(is.na(s[k])) return(NA)

    if(k>1) {
      sdiff = s[k] - s[k-1]
    } else {sdiff = 10e30}
    k = k+1

    if(k>kmax) {
      #vals = complex(real=A,im=2*pi*((kmax+1):(kmax+20)))
      #cat("************\n\n\nvals =", vals, "\n")
      ig[(kmax+1):(kmax+20)] = f(complex(real=A,imaginary=2*pi*((kmax+1):(kmax+20)))/(2*t))
      #ig[(kmax+1):(kmax+20)] = sapply(complex(real=A,imaginary=2*pi*((kmax+1):(kmax+20)))/(2*t),f)
      kmax = kmax + 20
    }
  }
  #cat("k =", k, "\n")
  return(s[k-1]*exp(A/2))
}


