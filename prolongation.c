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


double* prolongation(double *matrix){
	double 	*result;
	int	N_ = N*2-1;
	result = (double *)malloc(N_*N_*sizeof(double));
	for( int i=1;i<N_-1;i++ )
	for( int j=1;j<N_-1;j++ ){
		result[ind(i, j, N)] = matrix[ind(i/2, j/2, N)]
			         + ( matrix[ind(i/2+1, j/2, N)]
				   + matrix[ind(i/2-1, j/2, N)]
				   + matrix[ind(i/2, j/2+1, N)]
				   + matrix[ind(i/2, j/2-1, N)] )/2
		                 + ( matrix[ind(i/2+1, j/2+1, N)]
				   + matrix[ind(i/2-1, j/2-1, N)]
				   + matrix[ind(i/2+1, j/2-1, N)]
				   + matrix[ind(i/2-1, j/2+1, N)] )/4;
	}
	printf("test of prolongation\n");
	print(result,N_);
	return result;
}
