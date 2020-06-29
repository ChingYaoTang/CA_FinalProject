#include <cstdio>
#include <cstdlib>
#include <cstring>

#define N		11
#define GRID_SIZE	N
#define BLOCK_SIZE  	N


int ind( int i, int j, int NGrid ) {
        return i * NGrid + j;
}

void print( double *matrix, int n) {
    for( int i=0; i<n; i++ ) {
        for( int j=0; j<n; j++ ) {
            int index = i*n + j;
                printf("%.3f\t", matrix[index]);
        }
    printf("\n");
    }
}

__global__
void restriction_gpu( double (*matrix_f), int n_f, double (*matrix_c) ){
//      const int p = blockDim.x*blockIdx.x + threadIdx.x;
        int n_c = (n_f+1)/2;
        //int p   = blockDim.x*blockIdx.x + threadIdx.x;
        int i_c = blockIdx.x;
        int j_c = threadIdx.x;
        int i_f = 2*i_c;
        int j_f = 2*j_c;
        if( i_c<n_c && j_c<n_c){
        matrix_c[i_c*n_c+j_c] = matrix_f[i_c*n_c+j_c]/4
                              + ( matrix_f[(i_f+1)*n_f+j_f]
                                + matrix_f[(i_f-1)*n_f+j_f]
                                + matrix_f[i_f*n_f+(j_f+1)]
                                + matrix_f[i_f*n_f+(j_f-1)] )/8
                                + ( matrix_f[(i_f+1)*n_f+(j_f+1)]
                                  + matrix_f[(i_f-1)*n_f+(j_f-1)]
                                  + matrix_f[(i_f+1)*n_f+(j_f-1)]
                                  + matrix_f[(i_f-1)*n_f+(j_f+1)] )/16;
        }
}




int main( void  ){
	printf( "test restriction\n" );
        double *phi_corr_h_ = (double *)malloc( N * N * sizeof(double) );
        for( int i=0; i<N; i++) {
                for( int j=0; j<N; j++) {
                        if( i==0 || j==0 || i==N-1 || j==N-1) phi_corr_h_[ind( i, j, N )] = 0.0;
                        else phi_corr_h_[ind( i, j, N )] = 1.0;
                }
        }
        printf( "phi_corr_h\n" );
        print( phi_corr_h_, N );
        double *phi_corr_2h_ = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
        

	int n_f = N;
	int n_c = (n_f+1)/2;
//	restriction( phi_corr_h_, N, phi_corr_2h_ );
       	double (*d_matrix_f),(*d_matrix_c);
        cudaMalloc( &d_matrix_f, n_f*n_f*sizeof(double));
        cudaMalloc( &d_matrix_c, n_c*n_c*sizeof(double));
        cudaMemcpy( d_matrix_f, phi_corr_h_, n_f*n_f*sizeof(double), cudaMemcpyHostToDevice );
        restriction_gpu  <<< GRID_SIZE, BLOCK_SIZE >>> ( d_matrix_f, n_f, d_matrix_c );
        cudaMemcpy( phi_corr_2h_, d_matrix_c, n_c*n_c*sizeof(double), cudaMemcpyDeviceToHost );
        cudaFree(d_matrix_f);
        cudaFree(d_matrix_c);
        printf("Using gpu restrict.\n");
 
	
	printf( "phi_corr_2h after restriction\n" );
        print( phi_corr_2h_, (N+1)/2 );
        free(phi_corr_h_);
        free(phi_corr_2h_);
	
	return 0;
}
