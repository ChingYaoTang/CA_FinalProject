#ifndef CAL_RESIDUAL_H
#define CAL_RESIDUAL_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density,*residual;
void cal_residual( double *phi_guess, double *rho, double *residual, int n, double h );

#endif
