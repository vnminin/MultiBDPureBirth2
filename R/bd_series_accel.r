
levin.env <- new.env()

levin.env$levin_n <- 1
levin.env$levin_ncv <- 0
levin.env$levin_nmax <- 20
levin.env$levin_numer <- rep(NA,levin.env$levin_nmax)
levin.env$levin_denom <- rep(NA,levin.env$levin_nmax)
levin.env$levin_cnvgd <- 0
levin.env$levin_small <- 1e-8
levin.env$levin_big <- Inf
levin.env$levin_eps <- 1e-8
levin.env$levin_lastval <- 0
levin.env$levin_lasteps <- NA

levin_init = function(nmax, epss) {
  levin_n <<- 1
  levin_ncv <<- 0
  levin_lastval <<- 0
  levin_eps = epss
  levin_numer <<- rep(NA,nmax)
  levin_denom <<- rep(NA,nmax)
}

# s is the nth partial sum of the sequence
# omega is the nth remainder estimate
next_approx = function(s, omega, beta=1) {
  #cat("next(",s,",",omega,",",beta,")\n")
  if((s==0)&&(omega==0)) {return(0)}
  term = 1/(beta+levin_n-1)
  levin_denom[levin_n] <<- term/omega
  levin_numer[levin_n] <<- s*levin_denom[levin_n]
  if(levin_n > 1) {
    ratio = (beta + levin_n - 1)*term
    for(j in 1:levin_n ) {
      fact = (levin_n-1-j+beta)*term
      levin_numer[levin_n-j] <<- levin_numer[levin_n-j+1] - fact*levin_numer[levin_n-j]
      levin_denom[levin_n-j] <<- levin_denom[levin_n-j+1] - fact*levin_denom[levin_n-j]
      term = term * ratio
    }
  }
  levin_n <<- levin_n+1
  if(abs(levin_denom[1]) < levin_small) {
    val = levin_lastval
  } else {
    val = levin_numer[1]/levin_denom[1]
  }
  levin_lasteps <<- abs(val - levin_lastval)

  if(is.na(levin_lasteps)) {
    cat("next_approx(",s,",",omega,")\n")
    cat("levin_lasteps is NA\n")
    cat("val =",val,"\n")
    cat("levin_lastval =",levin_lastval,"\n")
    cat("levin_ncv =",levin_ncv,"\n")
    cat("levin_n =",levin_n,"\n")
    cat("levin_numer =",levin_numer,"\n")
    cat("levin_denom =",levin_denom,"\n")
    return(NA)
  }


  if(levin_lasteps <= levin_eps) {levin_ncv <<- levin_ncv + 1}
  if(levin_ncv >= 5) {
    levin_cnvgd <<- 1
  }
  levin_lastval <<- val
  return(val)
}

environment(levin_init) <- levin.env
environment(next_approx) <- levin.env


