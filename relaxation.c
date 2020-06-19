#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include "relative_error.h"
#define PI acos(-1)
extern const float L;
extern bool sor_method; 


// arguments: (1)phi matirx, (2)rho matrix, (3)size of the matrix,(4)convergence criterion, 
// 	      (5)updating method: 1="normal", 0="even odd", (6)omega for SOR (should be 1 for smoothing => GS), 
// 	      (7)which equation are we dealing with: 0 for Poisson eq., 1 for residual eq. 


void relaxation( double *phi_guess, double *rho, int n, double *conv_criterion, float omega, bool w ) {
//	Determine the physical grid size
	double h = L/(n-1);
//	Two end criteria for relaxation
	double *itera = (double *)malloc( sizeof(double) );
	*itera = 0;
	double *error = (double *)malloc( sizeof(double) );
	*error = 1;
//	Store the primitive input to make the comparison with the up-to-date result
	double *phi_old = (double *)malloc( n*n*sizeof(double) );

//	Set the end criterion
	double *condition1;
	double *condition2;
	if( *conv_criterion<1.0 ) {
		condition1 = error;
		condition2 = conv_criterion;
	} else {
		condition1 = conv_criterion;
		condition2 = itera;
	}

	if( sor_method==1 ) {
		while( *condition1 > *condition2 ) {
			*itera += 1;
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
						        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
				*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
			}
		}
	} else if( sor_method==0 ) {
 	     	int is, js;
		while( *condition1 > *condition2 ) {
			*itera += 1;
			*error = 0;
//	       		copy old potential
			memcpy( phi_old, phi_guess, n*n*sizeof(double) );
//      	        update by odd-even
			is = 1;
			for( int oe=1; oe<=2; oe++ ) {
				js = is;
#ifdef OPENMP
#pragma omp parallel for
#endif
				for( int i=1; i<(n-1); i++ ) {
//#ifdef OPENMP
//#pragma omp parallel for
//#endif
					for( int j=js; j<(n-1); j+=2 ) {
						phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
				    				             	+ phi_guess[ind(i-1, j, n)]
								        	     + phi_guess[ind(i, j+1, n)]
									             + phi_guess[ind(i, j-1, n)]
									             - phi_guess[ind(i, j, n)]*4
								        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
						*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
					}
					js = 3-js;
				}
				is = 3-is;
			}
//		relative_error( phi_guess, phi_old, n, error);
//		if(itera%100==1) {
//			printf("error in while = %g\n", *error);
//			print( phi_guess, n );
//			print( phi_old , n );
//		}
		}
	}
	if( *conv_criterion>1.0 ) {
		printf( "[N = %3d               ] Finish relaxation. Total iteration = %g, final conv error = %e\n", n, *itera, *error);
	} else {
		printf("Exact solver by relaxation terminated. Total iteration = %g, final conv error = %e\n", *itera, *error);
	}

	free( phi_old );
	free( error );
	free( itera );
}
