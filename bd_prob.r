
bd_prob = function(m,n,t,brates,drates) {
  fhat = function(s) {return(lt_eval_mod(m,n,s,brates,drates))}
  res = sapply(t,lt_invert,f=fhat)
  if(any(is.na(res))) cat("bd_prob(",m,",",n,",",t,",) failed\n")
  return(res)
}

