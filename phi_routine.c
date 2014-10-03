#include <R.h>
#include <complex.h>

void phi_routine(int *B, double *prod_mu2, double *prod_lambda2, double complex *Bk1dBk, double complex *BidBj, double complex *lentz, double complex *phi){
	int i,j;
	double complex fac,B1,B2,v;
	for (i=0; i<=B[0]; i++) {
		for (j=0; j<=B[0]; j++) {
			if (i<=j) {
				if (i==j) {
					fac = 1;
				} else {
					fac = prod_mu2[i+1+j*(B[0]+1)];
				}
				B1 = BidBj[i+j*(B[0]+1)];
				B2 = 1/Bk1dBk[j];
				v = fac*B1/(B2+lentz[j]);
			} else {
				fac = prod_lambda2[j+(i-1)*(B[0]+1)];
				B1 = BidBj[j+i*(B[0]+1)];
				B2 = 1/Bk1dBk[i];
				v = fac*B1/(B2+lentz[i]);
			}
			phi[i*(B[0]+1)+j] = v;		
		}
	}
}