#include <R.h>

void prod_vec(int *B,double *v, double *res) {
	int i,j;
	for (i=0; i<=B[0]; i++) {
		for (j=i; j<=B[0]; j++) {
			if (j==i) {
				res[i*(B[0]+1)+j] = v[j];	
			} 
			else res[i*(B[0]+1)+j] =res[i*(B[0]+1)+j-1]*v[j];  
		}
	}
}