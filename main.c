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


//	Set the basic parameters
const float L          = PI;                // Boxsize in the solver
const int   N          = 9;                   // Number of the resolution
float	    dx	       = L/(N-1);             // Spatial interval 
float	    bc	       = 0.0;                 // Boundary condition
double      error_conv = 1e-3;                 // Convergence error for the smoothing


int main( int argc, char *argv[] ) {

	double *analytic, *potential, *density, *residual;
	analytic  = (double *)malloc( N * N * sizeof(double) );
	potential = (double *)malloc( N * N * sizeof(double) );
	density   = (double *)malloc( N * N * sizeof(double) );
	residual  = (double *)malloc( N * N * sizeof(double) );
//	Initialize the Poisson solver problem
	init_sin( analytic, potential, density, 1.0, 1.0, 0.0 );
     	print(potential,N);	
//      Pre-smoothing up to certain error_conv
	relaxation( potential, density, N, dx, error_conv, 1.0 );
//	smoothing(potential,density,10);
	print(potential,N);
/*//      Calculate the residual in finest grid
	cal_residual( potential, density, residual, N, dx );      
//      Restrict the residual from h to 2h
        double *residual_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
	restriction( residual, N, residual_2h );
//      Solver exact solution
	double *phi_corr_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
	exact_im( residual_2h, (N+1)/2, phi_corr_2h );
//      Prolongation the phi_corr_2h to phi_corr_h
	double *phi_corr_h = (double *)malloc( N * N * sizeof(double) );
	prolongation( phi_corr_2h, (N+1)/2, phi_corr_h );	
//	Update potential
//	potential += phi_corr_h;

//      Post-smoothing
	relaxation( potential, density, N, dx, error_conv, 1.0 );

//	exactsolution((N+1)/2);	
*/
	free(analytic);
	free(potential);
	free(density);

	return EXIT_SUCCESS;
}

