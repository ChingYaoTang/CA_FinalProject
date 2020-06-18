#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"
#include "time.h"

// Input the fine gird matrix and calculate the corresponding coarse grid matrix by full weighting operator

// arguments: (1)fine matrix, (2)matrix size of fine matrix, (3)coarse matrix
void restriction( double *matrix_f, int n_f, double *matrix_c ) {
	clock_t tr;
	tr = clock();
	int n_c = (n_f+1)/2;
	int i_c, j_c, i_f, j_f;
//	Interior points
#ifdef OPENMP
#pragma omp parallel for collapse( 2 )
#endif

	for( i_c=1, i_f=2; i_c<n_c; i_c++, i_f+=2 )
	for( j_c=1, j_f=2; j_c<n_c; j_c++, j_f+=2 ) {
		matrix_c[ind(i_c, j_c, n_c)] = matrix_f[ind(i_f, j_f, n_f)]/4
				               + ( matrix_f[ind(i_f+1, j_f, n_f)]
				                 + matrix_f[ind(i_f-1, j_f, n_f)]
			                         + matrix_f[ind(i_f, j_f+1, n_f)]
				                 + matrix_f[ind(i_f, j_f-1, n_f)] )/8
				               + ( matrix_f[ind(i_f+1, j_f+1, n_f)]
		 		                 + matrix_f[ind(i_f-1, j_f-1, n_f)]
			 	                 + matrix_f[ind(i_f+1, j_f-1, n_f)]
				                 + matrix_f[ind(i_f-1, j_f+1, n_f)] )/16;
	}
//	Boundary points
//	Up & down boundaries
#ifdef OPENMP
#pragma omp parallel for collapse( 2 )
#endif
	for( j_c=0, j_f=0; j_c<n_c; j_c++, j_f+=2  ) {
		matrix_c[ind(0, j_c, n_c)]     = matrix_f[ind(0, j_f, n_f)];
		matrix_c[ind(n_c-1, j_c, n_c)] = matrix_f[ind(n_f-1, j_f, n_f)];
	}
//	Left & right boundaries
#ifdef OPENMP
#pragma omp parallel for collapse( 2 )
#endif
	for( i_c=0, i_f=0; i_c<n_c; i_c++, i_f+=2  ) {
		matrix_c[ind(i_c, 0, n_c)]     = matrix_f[ind(i_f, 0, n_f)];
		matrix_c[ind(i_c, n_c-1, n_c)] = matrix_f[ind(i_f, n_f-1, n_f)];
	}

	tr = clock()-tr;
	printf("[N_f = %3d -> N_c = %3d] Finish restriction.(Using %.3f sec) \n", n_f, n_c,tr/(double)CLOCKS_PER_SEC);
}
