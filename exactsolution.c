#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#define PI acos(-1)
extern const float L;
extern float dx, bc;
extern const int N;
extern double *analytic,*potential,*density;

int exactsolution(int n){
	double error=10.0;
	double *potential_bk;
	potential_bk = (double *)malloc(N*N*sizeof(double));
	memcpy(potential_bk,potential,sizeof(potential));
	while(error>1e-5){
		error=0.0;
		for( int i=1;i<(n-1);i++ )
		for( int j=1;j<(n-1);j++ ){
			potential[ind(i,j)] += potential[ind(i+1,j)]+potential[ind(i-1,j)]+potential[ind(i,j+1)]+potential[ind(i,j-1)]-4*potential[ind(i,j)]-dx*dx*density[ind(i,j)];
			error += fabs(potential_bk[ind(i,j)]-potential[ind(i,j)])/potential[ind(i,j)]/(N-2)/(N-2);
		}
	}
	return 0;
}
