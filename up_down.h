#ifndef UP_DOWN_H
#define UP_DOWN_H

void up( double *phi, double *rho, int final_level, int pulse_level, int *nn, int *level_ind, double *conv_loop );
void down( double *phi, double *rho, int pulse_level, int final_level, int *nn, int *level_ind, double *conv_loop, double *conv_precision );

#endif
