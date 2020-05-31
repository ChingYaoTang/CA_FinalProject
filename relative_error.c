#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

// arguments: (1)experimental value, (2)theoretical value, (3)matrix size
double relative_error( double *expe, double *theo, int n ) {
	double error = 0.;
	for( int i=1; i<n-1; i++ )
	for( int j=1; j<n-1; j++ ) {
		error += fabs( ( expe[ind(i, j, n)] - theo[ind(i, j, n)] ) / theo[ind(i, j, n)] )
	}
	error /= pow(n,2);
	return error;
}
