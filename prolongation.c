#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"
#include "time.h"

// Input the coarse gird matrix and calculate the corresponding fine grid matrix by bilinear operator

// arguments: (1)coarse matrix, (2)matrix size of fine matrix, (3)fine matrix
void prolongation( double *matrix_c, int n_c, double *matrix_f) {
	int n_f = 2*n_c-1;
	int i_c, j_c, i_f, j_f;
	clock_t t;
	t = clock();
//	Copy the points with factor 1 to the fine matrix
#ifdef OPENMP
#pragma omp parallel for collapse( 2 )
#endif
	for( i_f=0, i_c=0; i_f<n_f; i_f+=2, i_c++ )
	for( j_f=0, j_c=0; j_f<n_f; j_f+=2, j_c++ ) {
		matrix_f[ind(i_f, j_f, n_f)] = matrix_c[ind(i_c, j_c, n_c)];
	}
#ifdef OPENMP
#pragma omp parallel for collapse( 2 )
#endif
//	Compute the rest points
//	Compute even row
	for( i_f=0; i_f<n_f; i_f+=2 )
	for( j_f=1; j_f<n_f-1; j_f+=2 ) {
		matrix_f[ind(i_f, j_f, n_f)] = ( matrix_f[ind(i_f, j_f+1, n_f)] + matrix_f[ind(i_f, j_f-1, n_f)] )/2;
	}
#ifdef OPENMP
#pragma omp parallel for collapse( 2 )
#endif

//	Compute odd row
	for( i_f=1; i_f<n_f-1; i_f+=2 )
	for( j_f=0; j_f<n_f; j_f++ ) {
		matrix_f[ind(i_f, j_f, n_f)] = ( matrix_f[ind(i_f+1, j_f, n_f)] + matrix_f[ind(i_f-1, j_f, n_f)] )/2;
	}
	t = clock()-t;

	printf("[N_c = %3d -> N_f = %3d] Finish prolongation.(Use %.3f sec)\n", n_c, n_f,t/(double)CLOCKS_PER_SEC);
}
