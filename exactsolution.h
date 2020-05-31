#ifndef EXACTSOLUTION_H
#define EXACTSOLUTION_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
int exactsolution(int n);

#endif
