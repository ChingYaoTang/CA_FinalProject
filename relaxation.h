#ifndef RELAXATION_H
#define RELAXATION_H

extern const float L;
extern const bool sor_method; 
void relaxation( double *phi_guess, double *rho, int n, double *conv_criterion, float omega, bool w );

#endif
