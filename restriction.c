#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"


// arguments: (1)fine matrix, (2)matrix size of fine matrix, (3)coarse matrix
void restriction( double *matrix_f, int n_f, double *matrix_c ) {
	int	n_c = (n_f+1)/2;
	for( int i=1;i<n_c/2-1;i++ )
	for( int j=1;j<n_c/2-1;j++ ){
		matrix_c[ind(i, j, n_c)] = matrix_f[ind(2*i, 2*j, n_f)]/4
				        + ( matrix_f[ind(2*i+1, 2*j, n_f)]
				          + matrix_f[ind(2*i-1, 2*j, n_f)]
			                  + matrix_f[ind(2*i, 2*j+1, n_f)]
				          + matrix_f[ind(2*i, 2*j-1, n_f)] )/8
				        + ( matrix_f[ind(2*i+1, 2*j+1, n_f)]
		 		          + matrix_f[ind(2*i-1, 2*j-1, n_f)]
			 	          + matrix_f[ind(2*i+1, 2*j-1, n_f)]
				          + matrix_f[ind(2*i-1, 2*j+1, n_f)])/16;
	}
}
