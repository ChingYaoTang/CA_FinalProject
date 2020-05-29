//#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "add.h"
#include "init_sin.h"

//	Set the basic parameters
const float L   = 2*PI;                 // boxsize in the solver
const int   N   = 4;                   // Number of the resolution
float	    dx	= L/(N-1);
float	    bc	= 0.0;
double 	*analytic,*potential,*density;


int  main( int argc, char *argv[] ){

	init_sin(1.0,1.0,0.0);



	for( int i=0;i<N;i++ ){
		for( int j=0;j<N;j++ ){
			int index = i*N+j;
			printf("%.3f ",analytic[index]);
		}
		printf("\n");
	}


	free(analytic);
	free(potential);
	free(density);


	int a = 3;
	int b = 4;
	int c;
	c = add(a,b);
	printf("%d\n",c);

	return 0;
}

