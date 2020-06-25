#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"
extern const float L;
extern const double dx;
extern const int N;

void init_sin( double *analytic, double *potential, double *density, const double kx, const double ky, double bc ) {
	for( int i=0; i<N; i++) {
		for( int j=0; j<N; j++) {
			analytic[ind( i, j, N )] = sin( kx*i*dx ) * sin( ky*j*dx ) + bc;
			density[ind( i, j, N )]  = -( pow(kx,2) + pow(ky,2) ) * ( analytic[ind( i, j, N )] - bc );
			if( i==0 || j==0 || i==(N-1) || j==(N-1)) potential[ind( i, j, N )] = bc;
			else potential[ind( i, j, N )] = bc;
		}
	}
	printf("Using sin test problem, N=%d\n",N);
}

void init_sin_rho( double *density, const double kx, const double ky, double bc, int n ) {
	double h = L/(n-1);
	for( int i=0; i<n; i++) 
	for( int j=0; j<n; j++) {
		density[ind( i, j, n )]  = -( pow(kx,2) + pow(ky,2) ) * ( sin( kx*i*h ) * sin( ky*j*h ) );
	}
}

//	Test problem of eq.: L*u = 2x(y-1)(y-2x+xy+2)exp(x-y) with homogeneous BC
void init_test2_anal( double *analytic, int n ) {
	double h = L/(n-1);
	for( int i=0; i<n; i++) 
	for( int j=0; j<n; j++) {
		analytic[ind( i, j, n )]  = ( i*h )*( 1 - i*h )*( j*h )*( 1 - j*h )*exp( i*h -j*h );
	}
}

void init_test2_rho( double *density, int n ) {
	double h = L/(n-1);
	for( int i=0; i<n; i++) 
	for( int j=0; j<n; j++) {
		density[ind( i, j, n )]  = 2*( i*h )*( j*h - 1 )*( j*h - 2*i*h + i*h*j*h + 2 )*exp( i*h - j*h );
	}
}


