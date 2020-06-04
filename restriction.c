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


double* restriction(double *matrix){
	double 	*result;
	int	N_ = (N+1)/2;
	result = (double *)malloc(N_*N_*sizeof(double));
	for( int i=1;i<N_/2-1;i++ )
	for( int j=1;j<N_/2-1;j++ ){
		result[ind(i, j, N)] = matrix[ind(2*i, 2*j, N)]/4
				     + ( matrix[ind(2*i+1, 2*j, N)]
				       + matrix[ind(2*i-1, 2*j, N)]
			               + matrix[ind(2*i, 2*j+1, N)]
				       + matrix[ind(2*i, 2*j-1, N)] )/8
				     + ( matrix[ind(2*i+1, 2*j+1, N)]
		 		       + matrix[ind(2*i-1, 2*j-1, N)]
			 	       + matrix[ind(2*i+1, 2*j-1, N)]
				       + matrix[ind(2*i-1, 2*j+1, N)])/16;
	}
	printf("test of restriction\n");
	print(result,N_);
	return result;
}
