//input residual and size of the matrix, output phi^tilda^corr_gridsize. By method of fast Fourier transform.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const int N;
extern const double dx;
void exact_fft( double *phi_corr, double *residual, int n, double *conv_criterion, float omega, bool w ) {
    //compute Laplacian(phi_corr)=-residual

}