//input residual and size of the matrix, output phi^tilda^corr_gridsize. By method of fast Fourier transform.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"
#include <fftw3.h>

void exact_fft( double *phi_corr, double *residual, int n, double *conv_criterion, float omega, bool w ) {
    //compute Laplacian(phi_corr)=-residual
  
    /// Grid size
    int Nh=(n/2+1);

    /// Declare FFTW components.
    fftw_complex *residualk;
    residualk = (fftw_complex*) fftw_malloc( n * n * sizeof(fftw_complex));
    fftw_plan fwrd = fftw_plan_dft_r2r_2d(n,n,residual,residualk,FFTW_ESTIMATE);
    fftw_plan bwrd = fftw_plan_dft_c2r_2d(n,n,phi_corr,residualk,FFTW_ESTIMATE);
        
    fftw_execute(fwrd);

    double k1,k2;
    for (int i=0;i<n;i++){
        if (2*i<n)
            k1 = i;
        else
            k1 = n-i;                        
        for (int j=0;j<Nh;j++){
            k2 = j;
            double fac = pow(k1,2)+pow(k2,2);
            if (fabs(fac) < 1e-14){
                residualk[j+Nh*(i)][0] = 0.0;
                residualk[j+Nh*(i)][1] = 0.0;
            }
            else{
                residualk[j+Nh*(i)][0] /= fac;
                residualk[j+Nh*(i)][1] /= fac;
            }               
        }
    }
        
    fftw_execute(bwrd);

    fftw_destroy_plan(fwrd); fftw_destroy_plan(bwrd);
    fftw_free(residualk);
}