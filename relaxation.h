#ifndef RELAXATION_H
#define RELAXATION_H

extern double final_conv_rate;
extern const int N;
extern const float L;
extern const bool sor_method; 
void relaxation( double *phi_guess, double *rho, int n, double *conv_criterion, float omega, bool w );

#endif
