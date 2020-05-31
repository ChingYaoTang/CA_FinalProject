#ifndef CAL_RESIDUAL_H
#define CAL_RESIDUAL_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density,*residual;
int cal_residual(double *matrix, int n);

#endif
