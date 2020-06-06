#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include "relative_error.h"
#define PI acos(-1)
extern const float L;

// arguments: (1)phi matirx, (2)rho matrix, (3)size of the matrix,(4)convergence criteria, 
// 	      (5)omega for SOR (should be 1 for smoothing => GS), (6)1="normal", 2="even odd"
void relaxation( double *phi_guess, double *rho, int n, double conv_error, float omega, int method ) {
	double h = L/(n-1);
	int itera = 0;
	double *error = (double *)malloc( sizeof(double) );
	*error = 1;
	double *phi_old = (double *)malloc( n*n*sizeof(double) );
	if( method==1 ) {
		while( *error>conv_error ) {
			itera += 1;
			*error = 0;
//		       	copy old potential
			memcpy( phi_old, phi_guess, n*n*sizeof(double) );
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
		}
	} else if( method==2 ) {
//	while( itera< 1000 ) 
		while( *error>conv_error ) {
			itera += 1;
			*error = 0;
//	       		copy old potential
			memcpy( phi_old, phi_guess, n*n*sizeof(double) );
//      	        update odd part
			for( int i=1; i<(n-1); i++ )
			for( int j=( i%2 ) ? 1: 2; j<(n-1); j+=2 ) {
				phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
				    			             + phi_guess[ind(i-1, j, n)]
							             + phi_guess[ind(i, j+1, n)]
							             + phi_guess[ind(i, j-1, n)]
							             - phi_guess[ind(i, j, n)]*4
						        	     - rho[ind(i, j, n)] * pow(h,2) );
				*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
			}
//			update even part
			for( int i=1; i<(n-1); i++ )
			for( int j=( i%2 ) ? 2: 1; j<(n-1); j+=2 ) {
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
	}
	printf( "Total iteration = %d, final conv error = %e\n", itera, *error);
	free( phi_old );
	free( error );
	printf( "Finish relaxation with n = %d.\n", n);
}
