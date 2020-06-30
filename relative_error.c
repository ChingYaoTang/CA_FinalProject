#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "basic.h"
#include <omp.h>

extern const int N;
extern const int   NThread;

// arguments: (1)experimental value, (2)theoretical value, (3)matrix size

void relative_error( double *expe, double *theo, int n, double *error ) {
	double sum = 0;
#	ifdef OPENMP
	omp_set_num_threads( NThread );
#	pragma omp parallel for collapse(2) reduction( +:sum )	
#	endif
	for( int i=1; i<n-1; i++ )
	for( int j=1; j<n-1; j++ ) {
		sum += fabs( ( expe[ind(i, j, n)] - theo[ind(i, j, n)] ) / theo[ind(i, j, n)] );
	}

	*error = sum/pow(n,2);
//	if( n==N ) printf("Relative error = %g\n", *error);
}
