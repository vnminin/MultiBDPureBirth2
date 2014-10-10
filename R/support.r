#prod_vec0 <- function(a,B,mat) {
#	v = mat[a,]
#	res = matrix(NA,nrow=(B+1),ncol=B+1)
#	for (i in 1:(B+1))
#		for (j in i:(B+1)) {
#			if (j==i) res[i,j] = v[j] 
#			else {
#				res[i,j] = res[i,j-1]*v[j]
#			}
#		}
#	return(res)
#}

prod_vec <- function(a,B,mat) {
	v = mat[a,]
	res = matrix(0,nrow=B+1,ncol=B+1)
	tmp = .C("prod_vec",as.integer(B),as.double(v),as.double(as.vector(res)))
	res = matrix(tmp[[3]], nrow=B+1, byrow=T)
	return(res)
}