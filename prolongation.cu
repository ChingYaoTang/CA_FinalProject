#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include <omp.h>

// Input the coarse gird matrix and calculate the corresponding fine grid matrix by bilinear operator

	// arguments: (1)coarse matrix, (2)matrix size of fine matrix, (3)fine matrix

__global__
void prolongation_gpu( double (*matrix_c), int n_c, double (*matrix_f) ){
	int n_f = 2*n_c-1;
	const int i_f = blockIdx.x;
       	const int j_f = threadIdx.x;
	if( i_f<n_f && j_f<n_f ){	
		matrix_f[i_f*n_f+j_f] = 1.0;	
		/*if( i_f%2==0 && j_f%2==0 ){
			int i_c = i_f/2;
	                int j_c = j_f/2;
        	        matrix_f[i_f*n_f+j_f] = matrix_c[i_c*n_c+j_c];
		}
		__syncthreads();
		if( i_f%2==0 && j_f%2==1 ){
		matrix_f[i_f*n_f+j_f] = ( matrix_f[i_f*n_f+(j_f+1)] + matrix_f[i_f*n_f+(j_f-1)] )/2;
		}
		__syncthreads();
		if( i_f%2==1 ){
		matrix_f[i_f*n_f+j_f] = ( matrix_f[(i_f+1)*n_f+j_f] + matrix_f[(i_f-1)*n_f+j_f] )/2;
		}*/
	}
}


void prolongation( double *matrix_c, int n_c, double *matrix_f) {
#	ifdef OPENMP
	double t;
	t = omp_get_wtime();
#	endif	

	int n_f = 2*n_c-1;
	int i_c, j_c, i_f, j_f;
#	ifdef GPU
	double (*d_matrix_f),(*d_matrix_c);
        cudaMalloc( &d_matrix_f, n_f*n_f*sizeof(double));
        cudaMalloc( &d_matrix_c, n_c*n_c*sizeof(double));
        cudaMemcpy( d_matrix_c, matrix_c, n_c*n_c*sizeof(double), cudaMemcpyHostToDevice );
        prolongation_gpu  <<< GRID_SIZE, BLOCK_SIZE >>> ( d_matrix_c, n_c, d_matrix_f );
 	cudaMemcpy( matrix_f, d_matrix_f, n_f*n_f*sizeof(double), cudaMemcpyDeviceToHost );
        cudaFree(d_matrix_f);
        cudaFree(d_matrix_c);
        printf("Using gpu prolongate.\n");
	print(matrix_f,n_f);
#	endif



/*
//	Copy the points with factor 1 to the fine matrix
#	ifdef OPENMP
#	pragma omp parallel for collapse(2) private( i_c, j_c ) 
#endif
	for( i_f=0; i_f<n_f; i_f+=2 ) 
	for( j_f=0; j_f<n_f; j_f+=2 ) {
		i_c = i_f/2;
		j_c = j_f/2;
		matrix_f[ind(i_f, j_f, n_f)] = matrix_c[ind(i_c, j_c, n_c)];
	}
//	Compute the rest points
//	Compute even row
#ifdef OPENMP
#	pragma omp parallel for collapse( 2 )
#endif
	for( i_f=0; i_f<n_f; i_f+=2 )
	for( j_f=1; j_f<n_f-1; j_f+=2 ) {
		matrix_f[ind(i_f, j_f, n_f)] = ( matrix_f[ind(i_f, j_f+1, n_f)] + matrix_f[ind(i_f, j_f-1, n_f)] )/2;
	}

//	Compute odd row
#ifdef OPENMP
#	pragma omp parallel for collapse( 2 )
#endif
	for( i_f=1; i_f<n_f-1; i_f+=2 )
	for( j_f=0; j_f<n_f; j_f++ ) {
		matrix_f[ind(i_f, j_f, n_f)] = ( matrix_f[ind(i_f+1, j_f, n_f)] + matrix_f[ind(i_f-1, j_f, n_f)] )/2;
	}
*/	
#	ifdef DEBUG
	t = omp_get_wtime()-t;
	printf("[N_c = %4d -> N_f = %4d] Finish prolongation.(Duration = %.3f sec)\n", n_c, n_f, t);
#	endif

}
