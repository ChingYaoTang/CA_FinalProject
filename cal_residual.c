/// input matrix and the size of the matrix, calculate residual


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const float L;
extern float dx, bc;
extern const int N;
extern double *analytic,*potential,*density,*residual;


int cal_residual(double *matrix,int n){
	residual = (double *)malloc(n*n*sizeof(double));
	for( int i=1;i<n-1;i++ )
	for( int j=1;j<n-1;j++ ){
		residual[ind(i,j)] = matrix[ind(i+1,j)]+matrix[ind(i-1,j)]+matrix[ind(i,j+1)]+matrix[ind(i,j-1)]-4*matrix[ind(i,j)]-dx*dx*density[ind(i,j)];
	}
	printf("residual\n");
	print(residual,n);
	return 0;
}
