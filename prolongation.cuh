#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#ifndef PROLONGATION_H
#define PROLONGATION_H

void prolongation( double *matrix_c, int n_c, double *matrix_f );
__global__
void prolongation_gpu( double (*matrix_c), int n_c, double (*matrix_f) );

#endif
