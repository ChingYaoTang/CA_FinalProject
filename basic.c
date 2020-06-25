#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "prolongation.h"
#include "restriction.h"

//	1D index of the matrix element
int ind( int i, int j, int NGrid ) {
	return i * NGrid + j;
}


//	print out the matrix 
void print( double *matrix, int n) {
    for( int i=0; i<n; i++ ) {
        for( int j=0; j<n; j++ ) {
            int index = i*n + j;
                printf("%.3e\t", matrix[index]);
        }
    printf("\n");
    }
}


//	add the correction to phi_odd
void add_correction( double *phi_old, double *phi_corr, int n  ) {
	int i, j;
	for( i=0; i<n; i++)
	for( j=0; j<n; j++) {
		phi_old[ind(i, j, n)] += phi_corr[ind(i, j, n)];
	}
	printf("[N = %3d               ] Finish correction addition.\n", n);
}

//	fill zero
void fill_zero( double *phi_guess, int n  ) {
	int i, j;
	for( i=0; i<n; i++)
	for( j=0; j<n; j++) {
		phi_guess[ind(i, j, n)] = 0;
	}
	printf("[N = %3d               ] Fill zero.\n", n);
}

void test_prol_rest( const int N ) {
	printf( "test restriction\n" );
	double *phi_corr_h_ = (double *)malloc( N * N * sizeof(double) );
	for( int i=0; i<N; i++) {
		for( int j=0; j<N; j++) {
			if( i==0 || j==0 || i==N-1 || j==N-1) phi_corr_h_[ind( i, j, N )] = 0.0;
			else phi_corr_h_[ind( i, j, N )] = 1.0;
		}
	}
	printf( "phi_corr_h\n" );
	print( phi_corr_h_, N );
	double *phi_corr_2h_ = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
	restriction( phi_corr_h_, N, phi_corr_2h_ );
	printf( "phi_corr_2h after restriction\n" );
	print( phi_corr_2h_, (N+1)/2 );
	free(phi_corr_h_);
	free(phi_corr_2h_);

	printf( "test prolongation\n");
	double *phi_corr_2h = (double *)malloc( (N+1)/2 * (N+1)/2 * sizeof(double) );
	for( int i=0; i<(N+1)/2; i++) {
		for( int j=0; j<(N+1)/2; j++) {
			if( i==0 || j==0 || i==(N+1)/2-1 || j==(N+1)/2-1) phi_corr_2h[ind( i, j, (N+1)/2 )] = 0.0;
			else phi_corr_2h[ind( i, j, (N+1)/2 )] = 1.0;
		}
	}
//      Prolongate the phi_corr_2h to phi_corr_h
	double *phi_corr_h = (double *)malloc( N * N * sizeof(double) );
	printf( "phi_corr_2h\n" );
	print( phi_corr_2h, (N+1)/2 );
	prolongation( phi_corr_2h, (N+1)/2, phi_corr_h );
	printf( "phi_corr_h after prolongation\n" );
	print( phi_corr_h, N );
	free(phi_corr_h);
	free(phi_corr_2h);
}
