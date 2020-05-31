#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#define PI acos(-1)
extern const float L;
extern float dx, bc;

// arguements: (1)phi matirx, (2)rho matrix, (3)size of the matrix, (4)spatial interval,
//             (5)convergence criteria, (6)relaxation method
int relaxation( double *phi, double *rho, int n, double delta, double conv_error, int method ) {
	double error = 1.0; // give an arbitary initial convergence error
	double *phi_old = (double *)malloc( n*n*sizeof(double) );
	while( error>conv_error ) {
		memcpy( phi_old, phi, n*n); // copy old potential
		for( int i=1; i<(n-1); i++ )
		for( int j=1; j<(n-1); j++ ) {
			phi[ind(i, j, n)] += 1/4 * ( phi[ind(i+1, j, n)]
			    			       + phi[ind(i-1, j, n)]
						       + phi[ind(i, j+1, n)]
						       + phi[ind(i, j-1, n)]
						       - phi[ind(i, j, n)]*4
						       - rho[ind(i, j, n)] * pow(delta,2) );
		}
		error = relative_error( phi, phi_old, n);
	}
	return 0;
}
