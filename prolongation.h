#ifndef PROLONGATION_H
#define PROLONGATION_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
double* prolongation(double *matrix);

#endif
