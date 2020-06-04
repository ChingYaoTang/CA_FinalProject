#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"


// arguments: (1)fine matrix, (2)matrix size of fine matrix, (3)coarse matrix
void restriction( double *matrix_f, int n_f, double *matrix_c ) {
//	double 	*result;
	int	n_c = (n_f+1)/2;
//	result = (double *)malloc(N_*N_*sizeof(double));
	for( int i=1;i<n_/2-1;i++ )
	for( int j=1;j<n_/2-1;j++ ){
		matrix_c[ind(i, j, n_)] = matrix_f[ind(2*i, 2*j, n)]/4
				        + ( matrix_f[ind(2*i+1, 2*j, n)]
				          + matrix_f[ind(2*i-1, 2*j, n)]
			                  + matrix_f[ind(2*i, 2*j+1, n)]
				          + matrix_f[ind(2*i, 2*j-1, n)] )/8
				        + ( matrix_f[ind(2*i+1, 2*j+1, n)]
		 		          + matrix_f[ind(2*i-1, 2*j-1, n)]
			 	          + matrix_f[ind(2*i+1, 2*j-1, n)]
				          + matrix_f[ind(2*i-1, 2*j+1, n)])/16;
	}
	printf("test of restriction\n");
	print(result,N_);
//	return result;
}
