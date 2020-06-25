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

extern const float omega_sor;
extern cal_fn exact_solver;

//	down-sample with an exact solver
void down( double *phi, double *rho, int pulse_level, int final_level, int *nn, int *level_ind, double *conv_loop, double *conv_precision ) {
	//	distinguish the finest grid from others
	bool w;

	for( int l=pulse_level; l<final_level-1; l++ ) {
		printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d\n", l, l+1);
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
	exact_solver( (phi + level_ind[final_level-1]), (rho + level_ind[final_level-1]), nn[final_level-1], conv_precision, omega_sor, 1 );
	
	t_exact = omp_get_wtime()-t_exact;
	printf("Duration of exact solver = %.3f sec. \nUp-sample to previous level.\n", t_exact);

}

//	up-sample without an exact solver
void up( double *phi, double *rho, int final_level, int pulse_level, int *nn, int *level_ind, double *conv_loop ) {
	//	distinguish the finest grid from others
	bool w;
	
	for( int l=final_level-2; l>=pulse_level; l-- ) {	
		printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d\n", l+1, l);
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
		

		if ( l!=pulse_level || l==0 ) {
			//	Post-smooth the phi_old
			relaxation( (phi + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, w );
			printf("Up-sample to previous level.\n");
		}  
		

	}
}


//	Recursive W cycle
void W_cycle( double *phi, double *rho, int l, int final_level, int *nn, int *level_ind, double *conv_loop, double *conv_precision ) {	
	bool w;
	w=1;
	if( l==0 ) w=0 ;

	printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d                                gamma=1\n", l, l+1);
	//	Pre-smooth
	relaxation( (phi + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, w );	

	//	Calculate the residual
	double *residual = (double *)malloc( pow(nn[l],2) * sizeof(double) );
	cal_residual( (phi + level_ind[l]), (rho + level_ind[l]), (residual), nn[l], w );

	//	Restrict the residual, which is rho in next level
	restriction( (residual), nn[l], (rho + level_ind[l+1]) );

	//	Fill zeros in phi_old in next level
	fill_zero( (phi + level_ind[l+1]), nn[l+1] );


	if( (l+1) != (final_level-1) ) {
		W_cycle( phi, rho, l+1, final_level, nn, level_ind, conv_loop, conv_precision );
	}
	else {
		printf("----------------------------------------------------------------------------------------------------\n                                              Final Level:%d\n", l+1);
		exact_solver( (phi + level_ind[final_level-1]), (rho + level_ind[final_level-1]), nn[final_level-1], conv_precision, omega_sor, 1 );
	}
	

	printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d\n", l+1, l);
	//	Prolongate the phi_old from previous level, which is phi_corr in this level
	double *phi_corr = (double *)malloc( pow(nn[l],2) * sizeof(double) );
	prolongation( (phi + level_ind[l+1]), nn[l+1], (phi_corr) );	

	//	Update the phi_old of this level
	add_correction( (phi + level_ind[l]), (phi_corr), nn[l] );

	//	Re-smooth
	relaxation( (phi + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, w );	
	
	printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d                                gamma=2\n", l, l+1);
	//	Calculate the residual
	cal_residual( (phi + level_ind[l]), (rho + level_ind[l]), (residual), nn[l], w );

	//	Restrict the residual, which is rho in next level
	restriction( (residual), nn[l], (rho + level_ind[l+1]) );
	free(residual);


	if( (l+1) != (final_level-1) ) {
		W_cycle( phi, rho, l+1, final_level, nn, level_ind, conv_loop, conv_precision );
	} else { 
		printf("----------------------------------------------------------------------------------------------------\n                                              Final Level:%d\n", l+1);
		exact_solver( (phi + level_ind[final_level-1]), (rho + level_ind[final_level-1]), nn[final_level-1], conv_precision, omega_sor, 1 );
	}
	
	printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d\n", l+1, l);
	//	Prolongate the phi_old from previous level, which is phi_corr in this level
	prolongation( (phi + level_ind[l+1]), nn[l+1], (phi_corr) );	

	//	Update the phi_old of this level
	add_correction( (phi + level_ind[l]), (phi_corr), nn[l] );
	free(phi_corr);
	
	//	Post-smooth
	relaxation( (phi + level_ind[l]), (rho + level_ind[l]), nn[l], conv_loop, 1, w );

}



//	Down 1 step from level l
void down_1step( double *phi, double *rhs, int l, int *nn, int *level_ind, double *conv_loop, bool w ) {
	printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d\n", l, l+1);

	//	Pre-smooth the phi_old
	relaxation( (phi + level_ind[l]), (rhs + level_ind[l]), nn[l], conv_loop, 1, w );

	//	Calculate the residual
	double *residual = (double *)malloc( pow(nn[l],2) * sizeof(double) );
	cal_residual( (phi + level_ind[l]), (rhs + level_ind[l]), (residual), nn[l], w );

	//	Restrict the residual, which is rhs on next coarser level
	restriction( (residual), nn[l], (rhs + level_ind[l+1]) );
	free(residual);
	
	//	Fill zeros in phi_old on next coarser level
	fill_zero( (phi + level_ind[l+1]), nn[l+1] );

	

}
//	Up 1 step from level l
void up_1step( double *phi, double *rhs, int l, int *nn, int *level_ind, double *conv_loop, bool w ) {
	printf("----------------------------------------------------------------------------------------------------\n                                           Level:%d -> Level:%d\n", l, l-1);

	//	Prolongate the phi_old from coarser level l, which is phi_corr on finer level l-1
	double *phi_corr = (double *)malloc( pow(nn[l-1],2) * sizeof(double) );
	prolongation( (phi + level_ind[l]), nn[l], (phi_corr) );	

	//	Update the phi_old of finer level l-1
	add_correction( (phi + level_ind[l-1]), (phi_corr), nn[l-1] );
	free(phi_corr);
	

	//	Post-smooth the phi_old
	relaxation( (phi + level_ind[l-1]), (rhs + level_ind[l-1]), nn[l-1], conv_loop, 1, w );
	
}
