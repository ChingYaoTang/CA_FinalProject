#ifndef INIT_SIN_H
#define INIT_SIN_H

extern const float L;
extern const double dx, kx, ky, bc;
extern const int N;

void init_sin( double *analytic, double *potential, double *density );

#endif
