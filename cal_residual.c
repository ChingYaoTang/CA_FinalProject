// input matrices and their size, calculate residual by exact L
// res = L*phi_guess "-" rho.
// For original Poisson equation, this is just the case;
// but for successive residual equation, we must impose an extra minus sign on rho.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const float L;

// arguments: (1)phi_guess matrix, (2)rho matrix, (3)residual matrix, (4)matrix size, (5)which equation are we dealing with
void cal_residual( double *phi_guess, double *rho, double *residual, int n, bool w ) {
	double h = L/(n-1);

//      calculate interior points
//      0 for original Poisson equation, 1 for residual equation
	for( int i=1; i<n-1; i++ )
	for( int j=1; j<n-1; j++ ) {
		residual[ind(i, j, n)] = 1/pow(h,2) * ( phi_guess[ind(i+1, j, n)]
	  				 	      + phi_guess[ind(i-1, j, n)]
					 	      + phi_guess[ind(i, j+1, n)]
						      + phi_guess[ind(i, j-1, n)]
						      - phi_guess[ind(i, j, n)]*4 )
				         - rho[ind(i, j, n)]*pow(-1,w);
	}
//	impose homogeneous boundary condition
	for( int i=0; i<n; i++ ) {
		residual[ind(i, 0, n)] = residual[ind(i, n-1, n)] = residual[ind(0, i, n)] = residual[ind(n-1, i, n)] = 0.0;
	}

	printf("[N = %3d               ] Finish residual calculation.\n", n);
=======
	printf("Finish residual calculation with N = %d.\n", n);

}
