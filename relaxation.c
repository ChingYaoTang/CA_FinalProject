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
	int itera = 0;
	double *error = (double *)malloc( sizeof(double) );
	*error = 1;
	double *phi_old = (double *)malloc( n*n*sizeof(double) );
//	while( itera< 1000 ) {
	while( *error>conv_error ) {
		itera += 1;
		*error = 0;
//	       	copy old potential
//		memcpy( phi_old, phi_guess, n*n); // IT IS NOT WORKING!!
		for ( int i=0; i<n; i++ )
        	for ( int j=0; j<n; j++ ) {
        	   phi_old[ind(i, j, n)] = phi_guess[ind(i, j, n)];
        	}

		for( int i=1; i<(n-1); i++ )
		for( int j=1; j<(n-1); j++ ) {
			phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
			    			             + phi_guess[ind(i-1, j, n)]
						             + phi_guess[ind(i, j+1, n)]
						             + phi_guess[ind(i, j-1, n)]
						             - phi_guess[ind(i, j, n)]*4
						             - rho[ind(i, j, n)] * pow(h,2) );
			*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
		}
//		relative_error( phi_guess, phi_old, n, error);
//		if(itera%100==1) {
//			printf("error in while = %g\n", *error);
//			print( phi_guess, n );
//			print( phi_old , n );
//		}
	}
	printf( "loop: %d, final error = %e\n", itera, *error);
	free( phi_old );
	free( error );
	printf( "Finish relaxation.\n");
}
