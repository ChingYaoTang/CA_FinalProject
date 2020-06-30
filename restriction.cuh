

#ifndef RESTRICTION_H
#define RESTRICTION_H

void restriction( double *matrix_f, int n_f, double *matrix_c );
__global__
void restriction_gpu( double (*matrix_f), int n_f, double (*matrix_c) );



#endif
