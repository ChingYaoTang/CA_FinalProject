#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"

extern const float L;

// input matrices and their size, calculate residual by exact L
// res = L*phi_guess "-" rho.
// For original Poisson equation, this is just the case;
// but for successive residual equation, we must impose an extra minus sign on rho.


// arguments: (1)phi_guess matrix, (2)rho matrix, (3)residual matrix, (4)matrix size, (5)which equation are we dealing with
__global__
void residual_gpu( double *phi_guess, double *rho, double *residual, int n, bool w , double h) {
	const int i = blockIdx.x+1;
        const int j = threadIdx.x+1;
	residual[i*n+j] = 1/pow(h,2) * ( phi_guess[(i+1)*n+j]
                                       + phi_guess[(i-1)*n+j]
                                       + phi_guess[i*n+(j+1)]
                                       + phi_guess[i*n+(j-1)]
                                       - phi_guess[i*n+j]*4 )
                                       - rho[i*n+j]*pow(-1,w);
}

__global__
void zero_gpu( double *residual, int n){
	const int i = blockIdx.x;
        const int j = threadIdx.x;
        residual[i*n+j] = 0.0;
}



void cal_residual( double *phi_guess, double *rho, double *residual, int n, bool w ) {
	double h = L/(n-1);
/*
#ifdef GPU
	double *d_residual,*d_phi_guess,*d_rho;
	cudaMalloc( &d_residual, n*n*sizeof(double));
	cudaMalloc( &d_phi_guess, n*n*sizeof(double));
	cudaMalloc( &d_rho,  n*n*sizeof(double));
	cudaMemcpy( d_phi_guess, phi_guess, n*n*sizeof(double), cudaMemcpyHostToDevice );	
	cudaMemcpy( d_rho, rho, n*n*sizeof(double), cudaMemcpyHostToDevice );
	zero_gpu     <<<n,n>>> ( d_residual, n);
	residual_gpu <<<n-2,n-2>>> ( d_phi_guess, d_rho, d_residual, n, w ,h);
	cudaMemcpy( residual, d_residual, n*n*sizeof(double), cudaMemcpyDeviceToHost );
	cudaFree(d_residual);
	cudaFree(d_phi_guess);
	cudaFree(d_rho);
#else*/
//      calculate interior points
//      0 for original Poisson equation, 1 for residual equation
#	ifdef OPENMP
#	pragma omp parallel for collapse(2)
#	endif
	for( int i=1; i<n-1; i++ )
	for( int j=1; j<n-1; j++ ) {
		residual[ind(i, j, n)] = 1/pow(h,2) * ( phi_guess[ind(i+1, j, n)]
	  				 	      + phi_guess[ind(i-1, j, n)]
					 	      + phi_guess[ind(i, j+1, n)]
						      + phi_guess[ind(i, j-1, n)]
						      - phi_guess[ind(i, j, n)]*4 )
				         - rho[ind(i, j, n)]*pow(-1,w);
	}

//	impose homogeneous boundary condition
#	ifdef OPENMP
#	pragma omp parallel for
#	endif
	for( int i=0; i<n; i++ ) {
		residual[ind(i, 0, n)] = residual[ind(i, n-1, n)] = residual[ind(0, i, n)] = residual[ind(n-1, i, n)] = 0.0;
	}
//#endif
#	ifdef DEBUG
	printf("[N = %4d                ] Finish residual calculation.\n", n);
#	endif

}
