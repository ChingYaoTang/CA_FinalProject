#ifndef RESTRICTION_H
#define RESTRICTION_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
double* restriction(double *matrix);

#endif