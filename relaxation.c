#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include "relative_error.h"
#define PI acos(-1)


// arguments: (1)phi matirx, (2)rho matrix, (3)size of the matrix, (4)spatial interval,
//             (5)convergence criteria, (6)omega for SOR
void relaxation( double *phi_guess, double *rho, int n, double h, double conv_error, float omega ) {
	double error = 1.0; // give an arbitary initial convergence error
	int itera = 0;
	double *phi_old = (double *)malloc( n*n*sizeof(double) );
	while( itera< 1000 ) {
		itera += 1;
		memcpy( phi_old, phi_guess, n*n); // copy old potential
		for( int i=1; i<(n-1); i++ )
		for( int j=1; j<(n-1); j++ ) {
			phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
			    			             + phi_guess[ind(i-1, j, n)]
						             + phi_guess[ind(i, j+1, n)]
						             + phi_guess[ind(i, j-1, n)]
						             - phi_guess[ind(i, j, n)]*4
						             - rho[ind(i, j, n)] * pow(h,2) );
		}
		error = relative_error( phi_guess, phi_old, n);
		if(itera%50==1) printf("error:%e\n",error);
	}
	printf("error:%e\n",error);
	free(phi_old);
	printf("Finish relaxation.\n");
}
