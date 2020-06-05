#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)

int ind( int i, int j, int NGrid ) {
	return i * NGrid + j;
}


void print( double *matrix, int n) {
        for( int i=0; i<n; i++ ) {
                for( int j=0; j<n; j++ ) {
                        int index = i*n + j;
                        printf("%.3f\t", matrix[index]);
                }
                printf("\n");
        }
}

