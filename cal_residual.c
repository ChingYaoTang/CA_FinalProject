/// input matrices and their size, calculate residual by exact L

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const float L;

// arguments: (1)phi_guess matrix, (2)rho matrix, (3)residual matrix, (4)matrix size
void cal_residual( double *phi_guess, double *rho, double *residual, int n ) {
	double h = L/(n-1);
//      calculate interior points
	for( int i=1; i<n-1; i++ )
	for( int j=1; j<n-1; j++ ) {
		residual[ind(i, j, n)] = 1/pow(h,2) * ( phi_guess[ind(i+1, j, n)]
						      + phi_guess[ind(i-1, j, n)]
					 	      + phi_guess[ind(i, j+1, n)]
						      + phi_guess[ind(i, j-1, n)]
						      - phi_guess[ind(i, j, n)]*4 )
				         - rho[ind(i, j, n)];
	}
//	impose homogeneous boundary condition
	for( int i=0; i<n; i++ ) {
		residual[ind(i, 0, n)] = residual[ind(i, n-1, n)] = residual[ind(0, i, n)] = residual[ind(n-1, i, n)] = 0.0;
	}
	printf("Finish residual calculation with n = %d.\n", n);
//	print(residual,n);
}
