


eu = function(a,b,t,brates,drates,w=NULL) {
  if(is.null(w)) { w = function(z) {return(1)} }
  tol = 1e-8
  kmin = min(a,b)
  kmax = max(a,b)
  fconv = function(s) {
    psum = 0
    for(k in kmin:kmax) {
      term = w(k)*brates(k)*lt_eval_mod(a,k,s,brates,drates) * 
                                   lt_eval_mod(k+1,b,s,brates,drates)
      psum = psum + term
      #cat("k =", k, ", term[1] =", term[1], "\n")
    }
    if(kmin>0) {
      for(k in (kmin-1):0) {
        term = w(k)*brates(k)*lt_eval_mod(a,k,s,brates,drates) * 
                            lt_eval_mod(k+1,b,s,brates,drates)
        psum = psum + term
        #cat("k =", k, ", term[1] =", term[1], "\n")
        if(sum(abs(term)) < tol) break
      }
    }
    k = kmax+1
    term = 1e30
    while(sum(abs(term)) > tol) {
      term = w(k)*brates(k)*lt_eval_mod(a,k,s,brates,drates) * 
                            lt_eval_mod(k+1,b,s,brates,drates)

      psum = psum + term
      #cat("k =", k, ", term[1] =", term[1], "\n")
      k = k+1
    }
    return(psum)
  }
  val = lt_invert(fconv,t) / bd_prob(a,b,t,brates,drates) 
  return(val)
}


#######################################################

ed = function(a,b,t,brates,drates, w=NULL) {
  if(is.null(w)) {w = function(z) {return(1)}}
  tol = 1e-8
  kmin = min(a,b)
  kmax = max(a,b)
  fconv = function(s) {
    psum = 0
    for(k in kmin:kmax) {
      #cat("k =", k,"\n")
      if(k==0) next
      psum = psum + w(k)*drates(k)*lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k-1,b,s,brates,drates)
    }
    if(kmin>0) {
      for(k in (kmin-1):0) {
        #cat("k =", k,"\n")
        if(drates(k) < 1e-8) { term = 0}
        else {
          term = w(k)*drates(k)*lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k-1,b,s,brates,drates)
        }
        psum = psum + term
        if(sum(abs(term)) < tol) break
      }
    }
    k = kmax+1
    term = 1e30
    while(sum(abs(term)) > tol) {
      #cat("k =", k,"\n")
      term = w(k)*drates(k)*lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k-1,b,s,brates,drates)
      psum = psum + term
      k = k+1
    }
    #cat("-------------------------\n")
    return(psum)
  }
  return(lt_invert(fconv,t) / bd_prob(a,b,t,brates,drates) )
}


########################################################


et = function(a,b,t,brates,drates, w=NULL) {
  if(is.null(w)) {w = function(z){return(1)}}
  tol = 1e-8
  kmin = min(a,b)
  kmax = max(a,b)
  fconv = function(s) {
    psum = 0
    for(k in kmin:kmax) {
      #cat("k =", k,"\n")
      psum = psum + w(k)*lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k,b,s,brates,drates)
    }
    if(kmin>0) {
      for(k in (kmin-1):0) {
        #cat("k =", k,"\n")
        term = w(k)*lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k,b,s,brates,drates)
        psum = psum + term
        if(sum(abs(term)) < tol) break
      }
    }
    k = kmax+1
    term = 1e30
    while(sum(abs(term)) > tol) {
      #cat("k =", k,"\n")
      term = w(k)*lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k,b,s,brates,drates)
      psum = psum + term
      k = k+1
    }
    #cat("-------------------------\n")
    return(psum)
  }
  return(lt_invert(fconv,t) / bd_prob(a,b,t,brates,drates) )
}

###################################################


euk = function(a,b,t,k,brates,drates) {
  tol = 1e-8
  fconv = function(s) {
      lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k+1,b,s,brates,drates)
  }
  return(brates(k) * lt_invert(fconv,t) / bd_prob(a,b,t,brates,drates) )
}

edk = function(a,b,t,k,brates,drates) {
  tol = 1e-8
  fconv = function(s) {
      lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k-1,b,s,brates,drates)
  }
  return(drates(k) * lt_invert(fconv,t) / bd_prob(a,b,t,brates,drates) )
}

etk = function(a,b,t,k,brates,drates) {
  tol = 1e-8
  fconv = function(s) {
      lt_eval_mod(a,k,s,brates,drates) * lt_eval_mod(k,b,s,brates,drates)
  }
  return(lt_invert(fconv,t) / bd_prob(a,b,t,brates,drates) )
}

