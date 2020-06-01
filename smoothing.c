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

int smoothing(int n){
	for( int k=0;k<n;k++ ){
		for( int i=1;i<(N-1);i++ ){
			for( int j=1;j<(N-1);j++ ){
				potential[ind(i, j, n)] += potential[ind(i+1, j, n)]
							+ potential[ind(i-1, j, n)]
							+ potential[ind(i, j+1, n)]
							+ potential[ind(i, j-1, n)]
							- 4*potential[ind(i, j, n)]
							- dx*dx*density[ind(i, j, n)];
			}
		}
	}
	return 0;
}
