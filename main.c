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


//	Set the basic parameters
const float  L          = 1;                  // Boxsize in the solver
const int    N          = 9;                  // Number of the resolution
const double dx         = L/(N-1);            // Spatial interval 
int     error_criteria  = 10;                 // Iteration for the smoothing


int main( int argc, char *argv[] ) {
//	test_prol_rest(N);	

	double *analytic, *potential, *density;
	analytic  = (double *)malloc( N * N * sizeof(double) );
	potential = (double *)malloc( N * N * sizeof(double) );
	density   = (double *)malloc( N * N * sizeof(double) );
	
//	Initialize the Poisson solver problem
	const double bc         = 1.0;        // Boundary condition
	const double kx         = PI/L;
	const double ky         = PI/L;
	init_sin( analytic, potential, density, kx, ky, bc );
    //print( analytic, N );

//	which cycle do u wanna use?
//	1:two grid	2:V cycle	3:W cycle
	int cycle =1;
	if (cycle=1){
		for(int times=0;times<2;times++){
			double *error_ = (double *)malloc(sizeof(double));
			
			//	Pre-smoothing up to certain error_conv
			relaxation( potential, density, N, error_criteria, 1, 1 );
			relative_error( potential, analytic, N, error_ );

			//	Calculate the residual in finest grid
			double *residual_h = (double *)malloc( N * N * sizeof(double) );
			cal_residual( potential, density, residual_h, N );
			//print( residual, N );

			//	Restrict the residual from h to 2h
			double *residual_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
			restriction( residual_h, N, residual_2h );
			//print( residual_2h, (N+1)/2 );

			//	Solve exact solution of phi_corr_2h
			double *phi_corr_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
			exact_im( residual_2h, (N+1)/2, phi_corr_2h );

			//	Prolongate the phi_corr_2h to phi_corr_h
			double *phi_corr_h = (double *)malloc( N * N * sizeof(double) );
			prolongation( phi_corr_2h, (N+1)/2, phi_corr_h );	

			//	Update potential
			add_correction( potential, phi_corr_h, N );

			//	Post-smoothing
			relaxation( potential, density, N, error_criteria, 1, 1 );

			//	Compute error
			//print( potential, N );
			relative_error( potential, analytic, N, error_ );
			
			free(error_ );
			free(residual_h);
			free(residual_2h);
			free(phi_corr_h);
			free(phi_corr_2h);
		}
	}
	//a test for writing residual and phi_corr together
	if (cycle=2){
		for(int times=0;times<2;times++){
			double *error_ = (double *)malloc(sizeof(double));
			
			//	Pre-smoothing up to certain error_conv
			relaxation( potential, density, N, error_criteria, 1, 1 );
			relative_error( potential, analytic, N, error_ );

			//total length of residual = total length of phi_corr
			int tot_length = N*N + (N+1)/2 * (N+1)/2;
			double *residual_h = (double *)malloc( tot_length * sizeof(double) );

			//	Calculate the residual in finest grid
			cal_residual( potential, density, residual_h[0], N );
			//print( residual, N );

			//	Restrict the residual from h to 2h
			double *residual_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
			restriction( residual_h, N, residual_2h );
			//print( residual_2h, (N+1)/2 );

			//	Solve exact solution of phi_corr_2h
			double *phi_corr_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
			exact_im( residual_2h, (N+1)/2, phi_corr_2h );

			//	Prolongate the phi_corr_2h to phi_corr_h
			double *phi_corr_h = (double *)malloc( N * N * sizeof(double) );
			prolongation( phi_corr_2h, (N+1)/2, phi_corr_h );	

			//	Update potential
			add_correction( potential, phi_corr_h, N );

			//	Post-smoothing
			relaxation( potential, density, N, error_criteria, 1, 1 );

			//	Compute error
			//print( potential, N );
			relative_error( potential, analytic, N, error_ );
			
			free(error_ );
			free(residual_h);
			free(residual_2h);
			free(phi_corr_h);
			free(phi_corr_2h);
		}
	}
	if (cycle=3){
		for(int times=0;times<1;times++){
			double *error_ = (double *)malloc(sizeof(double));
			
			double *residual_h = (double *)malloc( N * N * sizeof(double) );
		}

	}
	if (cycle=4){
		for(int times=0;times<1;times++){
			double *error_ = (double *)malloc(sizeof(double));
		}
	}

	free(analytic);
	free(potential);
	free(density);
	
	return EXIT_SUCCESS;
}

