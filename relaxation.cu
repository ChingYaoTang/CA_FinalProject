#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include "relative_error.h"
#include <omp.h>

extern const float L;
extern const bool sor_method; 


// arguments: (1)phi matirx, (2)rho matrix, (3)size of the matrix,(4)convergence criterion, 
// 	      (5)updating method: 1="normal", 0="even odd", (6)omega for SOR (should be 1 for smoothing => GS), 
// 	      (7)which equation are we dealing with: 0 for Poisson eq., 1 for residual eq. 

__global__
void relaxation_gpu_odd( double (*phi_guess), double (*rho), int n, double omega, bool w, double h ){
	int job = n/BLOCK_SIZE+1;

	for( int a=0;a<job;a++ )
	for( int b=0;b<job;b++ ){
	const int i = blockIdx.x*job+a+1;
	const int j = threadIdx.x*job+b+1;
		if( i<n-1 && j<n-1 ){
//	Compute odd cells
		if( (i%2+j%2)%2==0 ){
		 double r = omega/4 * ( phi_guess[(i+1)*n+j]+ phi_guess[(i-1)*n+j]
				 		 +phi_guess[i*n+(j+1)]+ phi_guess[i*n+(j-1)]-phi_guess[i*n+j]*4\\
						 -rho[i*n+j]*pow(h,2)*pow(-1,w));
		 //error[i*n+j]     = fabs(r/(phi_guess[i*n+j]*pow(n,2)));
		 phi_guess[i*n+j] += r;
		}
		}
	}
}
__global__
void relaxation_gpu_even( double (*phi_guess), double (*rho), int n, double omega, bool w, double h ){
        int job = n/BLOCK_SIZE+1;

        for( int a=0;a<job;a++ )
        for( int b=0;b<job;b++ ){
        const int i = blockIdx.x*job+a+1;
        const int j = threadIdx.x*job+b+1;
                if( i<n-1 && j<n-1 ){
//      Compute odd cells
                if( (i%2+j%2)%2==1 ){
                 double r = omega/4 * ( phi_guess[(i+1)*n+j]+ phi_guess[(i-1)*n+j]
                                                 +phi_guess[i*n+(j+1)]+ phi_guess[i*n+(j-1)]-phi_guess[i*n+j]*4\\
                                                 -rho[i*n+j]*pow(h,2)*pow(-1,w));
                 //error[i*n+j]     = fabs(r/(phi_guess[i*n+j]*pow(n,2)));
                 phi_guess[i*n+j] += r;
                }
                }
        }
}

__global__
void compute_error( double (*error), int n, double (*result)){
	int i = threadIdx.x;
	for( int j=1;j<n;j++ )	error[i*n]+=error[i*n+j];
	__syncthreads();
	if( threadIdx.x == 0 ){
		*result = 0.0;
		for( int k=0;k<n;k++ ) *result += error[k*n];
	}
}

__global__
void relative_error_gpu( double *expe, double *theo, int n, double *error ) {
	int i = blockIdx.x;
	int j = threadIdx.x;
	error[i*n+j] = fabs( ( expe[i*n+j] - theo[i*n+j] ) / theo[i*n+j] );
}




void relaxation( double *phi_guess, double *rho, int n, double *conv_criterion, float omega, bool w ) {
#ifdef DEBUG
	double tr;
	tr	= omp_get_wtime();
#endif
//	Determine the physical grid size
	double h = L/(n-1);

//	Two end criteria for relaxation
	double *itera = (double *)malloc( sizeof(double) );
	*itera = 0;
	double *error = (double *)malloc( sizeof(double) );
	*error = 1;

//	Store the primitive input to make the comparison with the up-to-date result
	double *phi_old = (double *)malloc( n*n*sizeof(double) );

//	Set the end criterion
	double *condition1;
	double *condition2;
	if( *conv_criterion<1.0 ) {
		condition1 = error;
		condition2 = conv_criterion;
	} else {
		condition1 = conv_criterion;
		condition2 = itera;
	}

//	Relaxation
	if( sor_method==1 ) {
		while( *condition1 > *condition2 ) {
			*itera += 1;
			*error = 0;
//		       	copy old potential
			memcpy( phi_old, phi_guess, n*n*sizeof(double) );
			for( int i=1; i<(n-1); i++ )
			for( int j=1; j<(n-1); j++ ) {
				phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
			    				             + phi_guess[ind(i-1, j, n)]
							             + phi_guess[ind(i, j+1, n)]
							             + phi_guess[ind(i, j-1, n)]
						        	     - phi_guess[ind(i, j, n)]*4
						        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
				*error += fabs( ( phi_guess[ind(i, j, n)] - phi_old[ind(i, j, n)] ) / phi_old[ind(i, j, n)] );
			}
		}
	} else if( sor_method==0 ) {
#ifdef PARALLEL_GPU
		double (*d_phi), (*d_rho);//, (*d_error), (*d_result), (*d_phi_old);
#endif
		while( *condition1 > *condition2 ) {
			*itera += 1;
			*error = 0;
//	       		copy old potential
			memcpy( phi_old, phi_guess, n*n*sizeof(double) );
#ifdef PARALLEL_GPU
			double *error_tot;
			error_tot = (double *)malloc(n*n*sizeof(double));
			cudaMalloc( &d_phi, n*n*sizeof(double));
			//cudaMalloc( &d_phi_old, n*n*sizeof(double));
			cudaMalloc( &d_rho, n*n*sizeof(double));
			//cudaMalloc( &d_error, n*n*sizeof(double));
			//cudaMalloc( &d_result, sizeof(double));
			//cudaMemcpy( d_phi_old, phi_guess, n*n*sizeof(double), cudaMemcpyHostToDevice );
			cudaMemcpy( d_phi, phi_guess, n*n*sizeof(double), cudaMemcpyHostToDevice );
			cudaMemcpy( d_rho, rho, n*n*sizeof(double), cudaMemcpyHostToDevice );
			relaxation_gpu_odd <<< BLOCK_SIZE,GRID_SIZE >>> ( d_phi, d_rho, n, omega, w, h);// d_error);
			relaxation_gpu_even <<< BLOCK_SIZE,GRID_SIZE >>> ( d_phi, d_rho, n, omega, w, h);// d_error);
		//	relative_error_gpu <<<n,n>>> (d_phi,d_phi_old,n,d_error);
		//	compute_error	<<<1,n>>> ( d_error, n, d_result);
		//	cudaMemcpy( error, d_result, sizeof(double), cudaMemcpyDeviceToHost );
			cudaMemcpy( phi_guess, d_phi, n*n*sizeof(double), cudaMemcpyDeviceToHost );
		//	cudaMemcpy( error_tot, d_error, n*n*sizeof(double), cudaMemcpyDeviceToHost );
		/*	for( int i=1;i<(n-1);i++ )
			for( int j=1;j<(n-1);j++ ){
				*error+=error_tot[i*n+j]/pow(n,2);//(phi_guess[i*n+j]-phi_old[i*n+j])/phi_old[i*n+j]/(n*n);
			}*/
			cudaFree(d_phi);
	                //cudaFree(d_error);
        	        cudaFree(d_rho);
		//	printf("Using GPU.\n");
#endif

#ifdef WO_OMP
		//	printf("Not Using GPU.\n");
			
//			update odd part
			for( int i=1; i<(n-1); i++ )
 			for( int j=( i%2 + (i+1)%2*2 ); j<(n-1); j+=2 ) {
 				phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
 				    			             + phi_guess[ind(i-1, j, n)]
 							             + phi_guess[ind(i, j+1, n)]
 							             + phi_guess[ind(i, j-1, n)]
 							             - phi_guess[ind(i, j, n)]*4
 						        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
 			}
//			update even part
 			for( int i=1; i<(n-1); i++ )
 			for( int j=( (i+1)%2 + i%2*2 ); j<(n-1); j+=2 ) {
 				phi_guess[ind(i, j, n)] += omega/4 * ( phi_guess[ind(i+1, j, n)]
 				    			             + phi_guess[ind(i-1, j, n)]
 							             + phi_guess[ind(i, j+1, n)]
 							             + phi_guess[ind(i, j-1, n)]
 							             - phi_guess[ind(i, j, n)]*4
 						        	     - rho[ind(i, j, n)] * pow(h,2) * pow(-1,w) );
 			}
#endif
			relative_error( phi_guess, phi_old, n, error );
//#endif
		}//end of while
#ifdef PARALLE_GPU
//		cudaMemcpy( phi_guess, d_phi_old, n*n*sizeof(double), cudaMemcpyDeviceToHost );
//		cudaFree(d_phi_old);
 //       	cudaFree(d_error);
//		cudaFree(d_rho);
#endif
	}
#ifdef DEBUG	

	tr = omp_get_wtime()-tr;
	
#endif
	if( *conv_criterion>1.0 ) {
		printf( "[N = %4d                ] Finish relaxation. Total iteration = %g, final conv error = %e \n", n, *itera, *error);
	} else {
		printf("Exact solver by relaxation terminated. Total iteration = %g, final conv error = %e\n", *itera, *error);
#ifdef DEBUG
		printf("Duration of exact solver = %.3f sec. \n", tr);
#endif
	}

	free( phi_old );
	free( error );
	free( itera );
}
