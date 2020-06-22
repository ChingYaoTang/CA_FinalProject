#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "time.h"
#include "basic.h"
#include "relative_error.h"
#define PI acos(-1)
#include <omp.h>
//#define OPENMP
extern const float L;
extern bool sor_method; 


// arguments: (1)phi matirx, (2)rho matrix, (3)size of the matrix,(4)convergence criterion, 
// 	      (5)updating method: 1="normal", 0="even odd", (6)omega for SOR (should be 1 for smoothing => GS), 
// 	      (7)which equation are we dealing with: 0 for Poisson eq., 1 for residual eq. 


void relaxation( double *phi_guess, double *rho, int n, double *conv_criterion, float omega, bool w ) {
	double tr;
	tr	= omp_get_wtime();
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
		while( *condition1 > *condition2 ) {
			*itera += 1;
			*error = 0;
//	       		copy old potential
			memcpy( phi_old, phi_guess, n*n*sizeof(double) );
//      	        update odd part
#ifdef OPENMP
#pragma omp parallel for
#endif
			for( int i=1; i<(n-1); i++ )
			for( int j=1; j<(n/2+1); j++ ) {
				j = (j+i%2)*2;
				phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
				    			             + phi_guess[ind(i-1, j, n)]
							             + phi_guess[ind(i, j+1, n)]
							             + phi_guess[ind(i, j-1, n)]
							             - phi_guess[ind(i, j, n)]*4
						        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
//				*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
			}
#ifdef OPENMP
#pragma omp barrier
#pragma omp parallel for
#endif
//			update even part
			for( int i=1; i<(n-1); i++ )
			for( int j=1; j<(n/2+1); j++ ) {
				j = (j+(i+1)%2)*2;
				phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
				    			             + phi_guess[ind(i-1, j, n)]
							             + phi_guess[ind(i, j+1, n)]
							             + phi_guess[ind(i, j-1, n)]
							             - phi_guess[ind(i, j, n)]*4
						        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
//				*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
			}
			
			relative_error( phi_guess, phi_old, n, error );
		}
	}
	tr = omp_get_wtime()-tr;
	if( *conv_criterion>1.0 ) {
		printf( "[N = %3d               ] Finish relaxation. Total iteration = %g, final conv error = %e (%.3f sec)\n", n, *itera, *error, tr);///(double)CLOCKS_PER_SEC);
	} else {
		printf("Exact solver by relaxation terminated. Total iteration = %g, final conv error = %e\n", *itera, *error);
	}

	free( phi_old );
	free( error );
	free( itera );
}
