#ifndef INIT_SIN_H
#define INIT_SIN_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic;
void init_sin( double *analytic, double *potential, double *density, float a,float b, float c );

#endif
