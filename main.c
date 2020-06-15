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
const float  L                = 1; 	    	// Boxsize in the solver
const int    N                = 17;              // Number of the resolution
const double dx               = L/(N-1);	// Spatial interval 
const int    cycle_num        = 2;		// Number of cylces
int          cycle_type       = 3;		// 1:two grid, 2:W cycle, 3:V cycle, 4:SOR
int          final_level      = 3;		// Final level of V cycle or W cycle
bool         sor_method       = 0;		// 0:even-odd, 1:normal

int main( int argc, char *argv[] ) {
//	test_prol_rest(N);	

	double *conv_loop        = (double *)malloc( sizeof(double) );		// Criterion for the smoothing
	*conv_loop               = 10;
	double *conv_precision   = (double *)malloc( sizeof(double) );		// Criterion for exact relaxation solver
	*conv_precision          = 1e-14;

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
	if (cycle_type==1) {
		for( int times=0; times<cycle_num; times++ ) {
			printf("\nTwo-grid No.%d\n", times);	
			//	Pre-smoothing
			relaxation( potential, density, N, conv_loop, 1, 1, 0 );
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
			exact_im( residual_2h, (N+1)/2, phi_corr_2h );
			//print( phi_corr_2h, (N+1)/2 );

			//	Prolongate the phi_corr_2h to phi_corr_h
			double *phi_corr_h = (double *)malloc( N * N * sizeof(double) );
			prolongation( phi_corr_2h, (N+1)/2, phi_corr_h );	

			//	Update potential
			add_correction( potential, phi_corr_h, N );

			//	Post-smoothing
			relaxation( potential, density, N, conv_loop, 1, 1, 0 );

			//	Compute error
			//print( potential, N );
			relative_error( potential, analytic, N, error_rel );
			
			free(residual_h);
			free(residual_2h);
			free(phi_corr_h);
			free(phi_corr_2h);
		}
	}
/*
	//a test for writing residual and phi_corr together
	else if (cycle_type==2) {
		for( int times=0; times<cycle_num; times++ ) {
			
			//	Pre-smoothing up to certain error_criterion
			relaxation( potential, density, N, error_criterion, 1, 1, 0 );
			relative_error( potential, analytic, N, error_rel );

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
			relaxation( potential, density, N, error_criterion, 1, 1, 0 );

			//	Compute error
			//print( potential, N );
			relative_error( potential, analytic, N, error_rel );
			
			free(residual_h);
			free(residual_2h);
			free(phi_corr_h);
			free(phi_corr_2h);
		}
	}
*/
//	V cycle
	else if (cycle_type==3) {
		// total length of residual 
		// = total length of phi_old 
		// = total length of rho 
		// = total length of phi_corr - pow(nn[final_level-1],2)
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

		double *phi_old  = (double *)malloc( tot_length * sizeof(double) );
		double *phi_corr = (double *)malloc( (tot_length - pow(nn[final_level-1],2)) * sizeof(double) );
		double *residual = (double *)malloc( tot_length * sizeof(double) );
		double *rho      = (double *)malloc( tot_length * sizeof(double) );
		memcpy( (phi_old + level_ind[0] ), potential, pow(nn[0],2) * sizeof(double) );
		memcpy( (rho + level_ind[0] ), density, pow(nn[0],2) * sizeof(double) );
		
		for( int times=0; times<cycle_num; times++ ) {
			printf("====================================================================================================\nV cycle No.%d\n", times);

			//	down-sample until the level before final one
			for( int l=0; l<final_level-1; l++ ) {
				printf("----------------------------------------------------------------------------------------------------\nLevel:%d\n", l);
				bool w=1;
				if( l==0 ) w=0 ;
				//	Pre-smoothing
				relaxation( (phi_old + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, 1, w );
				if( l==0 ) relative_error( (phi_old + level_ind[l]), analytic, nn[l], error_rel );

				//	Calculate the residual
				cal_residual( (phi_old + level_ind[l]), (rho + level_ind[l]), (residual + level_ind[l]), nn[l], w );

				//	Restrict the residual, which is rho in next level
				restriction( (residual + level_ind[l]), nn[l], (rho + level_ind[l+1]) );
				printf("Down-sample to next level.\n");
			}
			printf("----------------------------------------------------------------------------------------------------\nReach the final level.\n");
			printf("Level:%d (Coarsest level)\n", final_level-1);
			//	Solve exact solution in coarsest level
			exact_im( (rho + level_ind[final_level-1]), nn[final_level-1], (phi_old + level_ind[final_level-1]) );
			//print((phi_old + level_ind[final_level-1]),nn[final_level-1]);
			
			//	Prolongate the phi_old, which is phi_corr in previous level
			prolongation( (phi_old + level_ind[final_level-1]), nn[final_level-1], (phi_corr + level_ind[final_level-2]) );	
			printf("Up-sample to previous level.\n");

			//	up-sample until the level before finest one
			for( int l=final_level-2; l>0; l-- ) {	
				printf("----------------------------------------------------------------------------------------------------\nLevel:%d\n", l);
				//	Update phi_old
				add_correction( (phi_old + level_ind[l]), (phi_corr + level_ind[l]), nn[l] );
	
				//	Post-smoothing phi_old
				relaxation( (phi_old + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, 1, 1 );

				//	Prolongate the phi_old, which is phi_corr in previous level
				prolongation( (phi_old + level_ind[l]), nn[l], (phi_corr + level_ind[l-1]) );	
				printf("Up-sample to previous level.\n");
			}

			printf("----------------------------------------------------------------------------------------------------\nLevel:%d (Finest level) \n", 0);
			//	Update phi_old in finest level
			add_correction( (phi_old + level_ind[0]), (phi_corr + level_ind[0]), nn[0] );

			//	Post-smoothing phi_old in finest level
			relaxation( (phi_old + level_ind[0]), (rho + level_ind[0]), nn[0], conv_loop, 1, 1, 0 );
			//	Compute error
			//print( potential, N );
			relative_error( (phi_old + level_ind[0]), analytic, nn[0], error_rel );
			
		}
		free(nn);
		free(level_ind);
		free(phi_old);
		free(phi_corr);
		free(residual);
		free(rho);

	}
//	SOR
	else if (cycle_type==4) {
		float omega = 1.5;
		printf("\nSOR:\nOmega = %f \nConv_precision = %e \nMethod = Even-Odd \n", omega, *conv_precision);	
		relaxation( potential, density, N, conv_precision, sor_method, omega, 0 );
		relative_error( potential, analytic, N, error_rel );
	}

	free(conv_loop);
	free(conv_precision);
	free(analytic);
	free(potential);
	free(density);
	free(error_rel);
	return EXIT_SUCCESS;
}

