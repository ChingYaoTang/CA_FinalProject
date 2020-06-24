#ifndef INIT_H
#define INIT_H

extern const float L;
extern const double dx;
extern const int N;

void init_sin( double *analytic, double *potential, double *density ,const double kx, const double ky, double bc );
void init_sin_rho( double *density ,const double kx, const double ky, double bc, int n );

#endif
