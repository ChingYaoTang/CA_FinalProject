#ifndef SMOOTHING_H
#define SMOOTHING_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
int smoothing(int n);

#endif
