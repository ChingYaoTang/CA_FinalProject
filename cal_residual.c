/// input spatial interval, matrices and their size, calculate residual by exact L

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const float L;
extern float dx, bc;
extern const int N;
extern double *analytic,*potential,*density;

// arguments: (1)phi_guess matrix, (2)rho matrix, (3)residual matrix, (4)matrix size, (5)spatial interval
void cal_residual( double *phi_guess, double *rho, double *residual, int n, double h ) {
//      start from 0 because residual matrix has no constrain on the boundary
	for( int i=0; i<n; i++ )
	for( int j=0; j<n; j++ ) {
		residual[ind(i, j, n)] = 1/pow(h,2) * ( phi_guess[ind(i+1, j, n)]
						      + phi_guess[ind(i-1, j, n)]
					 	      + phi_guess[ind(i, j+1, n)]
						      + phi_guess[ind(i, j-1, n)]
						      - phi_guess[ind(i, j, n)]*4 )
				         - rho[ind(i, j, n)];
	}
//	printf("residual\n");
//	print(residual,n);
}
