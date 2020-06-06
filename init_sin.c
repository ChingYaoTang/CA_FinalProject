#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"
extern const float L;
extern const double dx;
extern const int N;

void init_sin( double *analytic, double *potential, double *density, const double kx, const double ky, double bc) {
	for( int i=0; i<N; i++) {
		for( int j=0; j<N; j++) {
			analytic[ind( i, j, N )] = sin( kx*i*dx ) * sin( ky*j*dx ) + bc;
			density[ind( i, j, N )]  = -( pow(kx,2) + pow(ky,2) ) * ( analytic[ind( i, j, N )] - bc );
			if( i==0 || j==0 || i==(N-1) || j==(N-1)) potential[ind( i, j, N )] = bc;
			else potential[ind( i, j, N )] = 0.0;
		//	potential[ind( i, j, N )] = 0.3 * sin( kx*i*dx ) * sin( ky*j*dx ) + bc;
		}
	}
}
