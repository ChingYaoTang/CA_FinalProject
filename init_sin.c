#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"
extern const float L;
extern float dx, bc;
extern const int N;

void init_sin( double *analytic, double *potential, double *density, float a,float b, float c ) {
	for( int i=0; i<N; i++) {
		for( int j=0; j<N; j++) {
			analytic[ind( i, j, N )] = sin(a*(float)i*dx)*sin(b*(float)j*dx)+c;
			density[ind( i, j, N )]  = -(a*a+b*b)*(analytic[ind( i, j, N )]-c);
			if( i==0 || j==0 || i==(N-1) || j==(N-1)) potential[ind( i, j, N)] = bc;
                        else potential[ind( i,j, N )] = 0;
		}
	}
}
