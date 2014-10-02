#include <R.h>
#include <complex.h>

void cf_BidBj(int *B,double *xvec, double complex *yvec, double complex *Bk1dBk, double complex *res, double complex *ans) {
	int i,j;
	for (i=0; i<=B[0]; i++) {
		for (j=i; j<=B[0]; j++) {
			if (j==i) {
				ans[j] = 1;	
			} else if (j==(i+1)) {
				ans[j] = 1/Bk1dBk[j-1];	
			}
				else ans[j] = yvec[j-1]*ans[j-1] + xvec[j-1]*ans[j-2];
			res[i*(B[0]+1) + j] = 1/ans[j];
		}
	}
}