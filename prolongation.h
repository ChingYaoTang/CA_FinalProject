#ifndef PROLONGATION_H
#define PROLONGATION_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
void prolongation( double *matrix_c, int n_c, double *matrix_f );

#endif
