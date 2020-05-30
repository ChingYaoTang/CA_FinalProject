//#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "init_sin.h"
#include "smoothing.h"
#include "restriction.h"
#include "prolongation.h"
#include "basic.h"

//	Set the basic parameters
const float L   = 2*PI;                 // boxsize in the solver
const int   N   = 5;                   // Number of the resolution
float	    dx	= L/(N-1);
float	    bc	= 0.0;
double 	*analytic,*potential,*density;



int main( int argc, char *argv[] ){

	init_sin(1.0,1.0,0.0);		//Initialize the Poisson solver problem
	smoothing(3);			//Applying smoothing for several times
	printf("after:\n");
	print(analytic,N);
	restriction(analytic);
	prolongation(analytic);

	free(analytic);
	free(potential);
	free(density);

	return 0;
}

