//#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "init_sin.h"
#include "smoothing.h"

//	Set the basic parameters
const float L   = 2*PI;                 // boxsize in the solver
const int   N   = 4;                   // Number of the resolution
float	    dx	= L/(N-1);
float	    bc	= 0.0;
double 	*analytic,*potential,*density;


int print(double *matrix){
	for( int i=0;i<N;i++ ){
                for( int j=0;j<N;j++ ){
                        int index = i*N+j;
                        printf("%.3f ",matrix[index]);
                }
                printf("\n");
        }
}


int main( int argc, char *argv[] ){

	init_sin(1.0,1.0,0.0);		//Initialize the Poisson solver problem
	smoothing(3);			//Applying smoothing for several times
	printf("after:\n");
	print(potential);



	free(analytic);
	free(potential);
	free(density);

	return 0;
}

