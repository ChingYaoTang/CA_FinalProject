#ifndef INIT_SIN_H
#define INIT_SIN_H

void init_sin( double *analytic, double *potential, double *density ,const double kx, const double ky, double bc );
void init_sin_d( double *density ,const double kx, const double ky, double bc, int n );

#endif
