#ifndef EXACT_IM_H
#define EXACT_IM_H

double comatrix(double *matrix, int p, int q, int m);
double det(double *matrix, int m);
double inversematrix(double *matrix, int n);
void exact_im( double *phi_corr, double *residual, int n, double *conv_criterion, float omega, bool w );

#endif
