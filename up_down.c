#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include "relaxation.h"
#include "cal_residual.h"
#include "restriction.h"
#include "prolongation.h"
#include "exact_im.h"
#include "basic.h"
#include <omp.h>

extern cal_fn exact_solver;

//	down-sample with an exact solver
void down( double *phi, double *rho, int pulse_level, int final_level, int *nn, int *level_ind, double *conv_loop, double *conv_precision ) {
	//	distinguish the finest grid from others
	bool w;

	for( int l=pulse_level; l<final_level-1; l++ ) {
		printf("----------------------------------------------------------------------------------------------------\n                                                Level:%d\n", l);
		w=1;
		if( l==0 ) w=0 ;

		//	Pre-smooth the phi_old
		relaxation( (phi + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, w );

		//	Calculate the residual
		double *residual = (double *)malloc( pow(nn[l],2) * sizeof(double) );
		cal_residual( (phi + level_ind[l]), (rho + level_ind[l]), (residual), nn[l], w );

		//	Restrict the residual, which is rho in next level
		restriction( (residual), nn[l], (rho + level_ind[l+1]) );
		free(residual);

		printf("Down-sample to next level.\n");
	}
	
	printf("----------------------------------------------------------------------------------------------------\n                                                Level:%d\n", final_level-1);
	printf("Reach the final level N= %d (Coarsest level).\n", nn[final_level-1]);
	
	double t_exact;
	t_exact = omp_get_wtime();
	
	//	Solve exact solution in coarsest level
	exact_solver( (phi + level_ind[final_level-1]), (rho + level_ind[final_level-1]), nn[final_level-1], conv_precision, 1, 1 );
	
	t_exact = omp_get_wtime()-t_exact;
	printf("Duration of exact solver = %.3f sec. \nUp-sample to previous level.\n", t_exact);

}

//	up-sample without an exact solver
void up( double *phi, double *rho, int final_level, int pulse_level, int *nn, int *level_ind, double *conv_loop ) {
	//	distinguish the finest grid from others
	bool w;
	
	for( int l=final_level-2; l>=pulse_level; l-- ) {	
		printf("----------------------------------------------------------------------------------------------------\n                                                Level:%d\n", l);
		w=1;	
		if( l==0 ) {
			w=0;
			printf("Back to the original level (Finest level).\n");
		}

		//	Prolongate the phi_old from previous level, which is phi_corr in this level
		double *phi_corr = (double *)malloc( pow(nn[l],2) * sizeof(double) );
		prolongation( (phi + level_ind[l+1]), nn[l+1], (phi_corr) );	

		//	Update the phi_old of this level
		add_correction( (phi + level_ind[l]), (phi_corr), nn[l] );
		free(phi_corr);
	
		//	Post-smooth the phi_old
		relaxation( (phi + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, w );
		if ( l!=pulse_level ) printf("Up-sample to previous level.\n");
	}
}
