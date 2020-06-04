#ifndef RELAXATION_H
#define RELAXATION_H

extern float dx, bc;
extern const int N;
extern const float L;
void relaxation( double *phi_guess, double *rho, int n, double h, double conv_error, float omega );

#endif
