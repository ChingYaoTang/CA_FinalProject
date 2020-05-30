#ifndef BASIC_H
#define BASIC_H

extern float dx, bc;
extern const int N;
extern const float L;
extern double *analytic,*potential,*density;
int ind(int i,int j);
int print(double *matric,int n);

#endif
