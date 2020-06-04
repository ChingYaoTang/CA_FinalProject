#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const float L;
extern float dx, bc;
extern const int N;


// arguments: (1)coarse matrix, (2)matrix size of fine matrix, (3)fine matrix
void prolongation( double *matrix_c, int n_c, double matrix_f) {
//	double 	*result;
	int	n_f = n_c*2-1;
//	result = (double *)malloc(n_f*n_f*sizeof(double));
	for( int i=1; i<n_f-1; i++ )
	for( int j=1; j<n_f-1; j++ ) {
		matrix_f[ind(i, j, n_f)] = matrix_c[ind(i/2, j/2, n_c)]
			        	 + ( matrix_c[ind(i/2+1, j/2, n_c)]
				  	   + matrix_c[ind(i/2-1, j/2, n_c)]
				   	   + matrix_c[ind(i/2, j/2+1, n_c)]
				   	   + matrix_c[ind(i/2, j/2-1, n_c)] )/2
		                	 + ( matrix_c[ind(i/2+1, j/2+1, n_c)]
			 	   	   + matrix_c[ind(i/2-1, j/2-1, n_c)]
					   + matrix_c[ind(i/2+1, j/2-1, n_c)]
				   	   + matrix_c[ind(i/2-1, j/2+1, n_c)] )/4;
	}
	printf("test of prolongation\n");
//	print(result,n_);
//	return result;
}
