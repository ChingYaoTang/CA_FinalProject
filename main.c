//#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "init_sin.h"			
//usage:( analytic, potential, density, a, b, c)
#include "restriction.h"
#include "prolongation.h"
#include "basic.h"
#include "cal_residual.h"
#include "relaxation.h"			
//usage:( potential, density, size of matrix, spatial interval,convergence criteria, omega of SOR)
#include "relative_error.h"
//usage:double value of error = ((1)experimental value, (2)theoretical value, (3)matrix size)
#include "exact_im.h"
#include "up_down.h"

//	Set the basic parameters
const float  L                = 1; 	    	// Boxsize in the solver
const int    N                = 65;              // Number of the resolution
const double dx               = L/(N-1);	// Spatial interval 
const int    cycle_num        = 1;		// Number of cylces
int          cycle_type       = 3;		// 1:two grid, 2:V cycle, 3:W cycle, 4:SOR
int          final_level      = 4;		// Final level of V cycle or W cycle
bool         sor_method       = 0;		// 0:even-odd, 1:normal
float        omega            = 1.5;
cal_fn       exact_solver     = relaxation;	// function name of the exact solver


//main function
int main( int argc, char *argv[] ) {
//	test_prol_rest(N);	

	double *conv_loop        = (double *)malloc( sizeof(double) );		// Criterion for the smoothing
	*conv_loop               = 10;
	double *conv_precision   = (double *)malloc( sizeof(double) );		// Criterion for exact relaxation solver
	*conv_precision          = 1e-10;

	double *analytic         = (double *)malloc( N * N * sizeof(double) );	// analytic potential matrix
	double *potential        = (double *)malloc( N * N * sizeof(double) );	// potential matrix of initial guess
	double *density          = (double *)malloc( N * N * sizeof(double) );	// density matrix
	double *error_rel        = (double *)malloc( sizeof(double) );		// Relative error with analytic solution
	
//	Initialize the Poisson solver problem
	const double bc         = 1.0;        // Boundary condition
	const double kx         = PI/L;
	const double ky         = PI/L;
	init_sin( analytic, potential, density, kx, ky, bc );
	//print( analytic, N );

//	which cycle do u wanna use?
// 	two grid
	if (cycle_type==1) {
		for( int times=0; times<cycle_num; times++ ) {
			printf("\nTwo-grid No.%d\n", times);	
			//	Pre-smoothing
			relaxation( potential, density, N, conv_loop, 1, 0 );
			relative_error( potential, analytic, N, error_rel );

			//	Calculate the residual in finest grid
			double *residual_h = (double *)malloc( N * N * sizeof(double) );
			cal_residual( potential, density, residual_h, N, 0 );
			//print( residual, N );

			//	Restrict the residual from h to 2h
			double *residual_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
			restriction( residual_h, N, residual_2h );
			//print( residual_2h, (N+1)/2 );

			//	Solve exact solution of phi_corr_2h
			double *phi_corr_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
			exact_solver( phi_corr_2h, residual_2h, (N+1)/2, conv_precision, omega, 1 );
			//print( phi_corr_2h, (N+1)/2 );

			//	Prolongate the phi_corr_2h to phi_corr_h
			double *phi_corr_h = (double *)malloc( N * N * sizeof(double) );
			prolongation( phi_corr_2h, (N+1)/2, phi_corr_h );	

			//	Update potential
			add_correction( potential, phi_corr_h, N );

			//	Post-smoothing
			relaxation( potential, density, N, conv_loop, 1, 0 );

			//	Compute error
			//print( potential, N );
			relative_error( potential, analytic, N, error_rel );
			
			free(residual_h);
			free(residual_2h);
			free(phi_corr_h);
			free(phi_corr_2h);
		}
	}

//	V cycle 
	else if (cycle_type==2) {
		int tot_length = 0;
		int *nn        = (int *)malloc( final_level * sizeof(int));
		int *level_ind = (int *)malloc( final_level * sizeof(int));
		int m          = N;
		for( int l=0; l<final_level; l++ ) {
			nn[l]        =  m;		// dim for each level
			level_ind[l] =  tot_length;	// index for each level
			tot_length   += pow(m,2);	// total length
			m            = (m+1)/2;
		}
		
		double *phi      = (double *)malloc( tot_length * sizeof(double) );
		double *rho      = (double *)malloc( tot_length * sizeof(double) );
		memcpy( (phi + level_ind[0] ), potential, pow(nn[0],2) * sizeof(double) );
		memcpy( (rho + level_ind[0]), density, pow(nn[0],2) * sizeof(double) );

		for( int times=0; times<cycle_num; times++ ) {
			printf("====================================================================================================\n                                             V cycle No.%d\n====================================================================================================\n", times+1);

			down(phi,rho, 0, final_level, nn, level_ind, conv_loop, conv_precision);
			up(phi,rho, final_level, 0, nn, level_ind, conv_loop);
			
			//	Compute error
			//	print( potential, N );
			relative_error( (phi + level_ind[0]), analytic, nn[0], error_rel );
			
		}
		free(nn);
		free(level_ind);
		free(phi);
		free(rho);
	}

//	W cycle
	else if (cycle_type==3) {
		int tot_length = 0;
		int *nn        = (int *)malloc( final_level * sizeof(int));
		int *level_ind = (int *)malloc( final_level * sizeof(int));
		int m          = N;
		for( int l=0; l<final_level; l++ ) {
			nn[l]        =  m;		// dim for each level
			level_ind[l] =  tot_length;	// index for each level
			tot_length   += pow(m,2);	// total length
			m            = (m+1)/2;
		}
		
		double *phi      = (double *)malloc( tot_length * sizeof(double) );
		double *rho      = (double *)malloc( tot_length * sizeof(double) );
		memcpy( (phi + level_ind[0] ), potential, pow(nn[0],2) * sizeof(double) );
		memcpy( (rho + level_ind[0]), density, pow(nn[0],2) * sizeof(double) );

		for( int times=0; times<cycle_num; times++ ) {
			printf("====================================================================================================\n                                             W cycle No.%d\n====================================================================================================\n", times+1);

			down(phi, rho, 0, final_level, nn, level_ind, conv_loop, conv_precision);
			
			for (int i=final_level-3; abs(i)<final_level-2; i--){
				up(phi, rho, final_level, abs(i)+1, nn, level_ind, conv_loop);
				down(phi, rho, abs(i)+1, final_level, nn, level_ind, conv_loop, conv_precision);
			}
			
			up(phi, rho, final_level, 0, nn, level_ind, conv_loop);
			
			//	Compute error
			//	print( potential, N );
			relative_error( (phi + level_ind[0]), analytic, nn[0], error_rel );
			
		}
		free(nn);
		free(level_ind);
		free(phi);
		free(rho);
	}

//	SOR
	else if (cycle_type==4) {
		printf("\nSOR:\nOmega = %f \nConv_precision = %e \n", omega, *conv_precision);	
		relaxation( potential, density, N, conv_precision, omega, 0 );
		relative_error( potential, analytic, N, error_rel );
	}

//	V cycle prototype
	else if (cycle_type==5) {
	}

	free(conv_loop);
	free(conv_precision);
	free(analytic);
	free(potential);
	free(density);
	free(error_rel);
	return EXIT_SUCCESS;
}

