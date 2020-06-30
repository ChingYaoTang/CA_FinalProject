#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include <omp.h>

extern const int   NThread;
// Input the fine gird matrix and calculate the corresponding coarse grid matrix by full weighting operator

// arguments: (1)fine matrix, (2)matrix size of fine matrix, (3)coarse matrix

void restriction( double *matrix_f, int n_f, double *matrix_c ) {
#	ifdef DEBUG
	double tr;
	tr = omp_get_wtime();
#	endif
	int n_c = (n_f+1)/2;
	int i_c, j_c, i_f, j_f;
	
//	Interior points
#	ifdef OPENMP
	omp_set_num_threads( NThread );
#	pragma omp parallel for collapse( 2 ) private( i_f, j_f )
#	endif
	for( i_c=1; i_c<n_c; i_c++ )
	for( j_c=1; j_c<n_c; j_c++ ) {
		i_f = 2*i_c;
		j_f = 2*j_c;
#		ifdef	FULL_WEIGHTING
		matrix_c[ind(i_c, j_c, n_c)] = matrix_f[ind(i_f, j_f, n_f)]/4
				               + ( matrix_f[ind(i_f+1, j_f, n_f)]
				                 + matrix_f[ind(i_f-1, j_f, n_f)]
			                         + matrix_f[ind(i_f, j_f+1, n_f)]
				                 + matrix_f[ind(i_f, j_f-1, n_f)] )/8
				               + ( matrix_f[ind(i_f+1, j_f+1, n_f)]
		 		                 + matrix_f[ind(i_f-1, j_f-1, n_f)]
			 	                 + matrix_f[ind(i_f+1, j_f-1, n_f)]
				                 + matrix_f[ind(i_f-1, j_f+1, n_f)] )/16;
#		endif
#		ifdef	HALF_WEIGHTING
		matrix_c[ind(i_c, j_c, n_c)] = matrix_f[ind(i_f, j_f, n_f)]/2
				               + ( matrix_f[ind(i_f+1, j_f, n_f)]
				                 + matrix_f[ind(i_f-1, j_f, n_f)]
			                         + matrix_f[ind(i_f, j_f+1, n_f)]
				                 + matrix_f[ind(i_f, j_f-1, n_f)] )/8;
#		endif
	}

//	Boundary points
#	ifdef OPENMP
#	pragma omp parallel for private( i_f )
#	endif
	for( i_c=0; i_c<n_c; i_c++  ) {
		i_f = i_c*2;
//	Up & down boundaries
		matrix_c[ind(i_c, 0, n_c)]     = matrix_f[ind(i_f, 0, n_f)];
		matrix_c[ind(i_c, n_c-1, n_c)] = matrix_f[ind(i_f, n_f-1, n_f)];
		
//	Left & right boundaries
		matrix_c[ind(0, i_c, n_c)]     = matrix_f[ind(0, i_f, n_f)];
		matrix_c[ind(n_c-1, i_c, n_c)] = matrix_f[ind(n_f-1, i_f, n_f)];
//		printf("Using openmp with number of threads = %d\n", omp_get_num_threads());
	}

#	ifdef DEBUG
	tr = omp_get_wtime()-tr;

#	ifdef	FULL_WEIGHTING
	printf("[N_f = %4d -> N_c = %4d] Finish full-weighting restriction.(Duration = %.3f sec) \n", n_f, n_c, tr);
#	endif
#	ifdef	HALF_WEIGHTING
	printf("[N_f = %4d -> N_c = %4d] Finish half-weighting restriction.(Duration = %.3f sec) \n", n_f, n_c, tr);
#	endif

#	endif

}


