#ifndef RELAXATION_H
#define RELAXATION_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
int relaxation( double *phi, double *rho, int n, double delta, double conv_error, int method );

#endif
