#ifndef UP_DOWN_H
#define UP_DOWN_H

extern const float omega_sor;
extern cal_fn exact_solver;

void up( double *phi, double *rho, int final_level, int pulse_level, int *nn, int *level_ind, double *conv_loop );
void down( double *phi, double *rho, int pulse_level, int final_level, int *nn, int *level_ind, double *conv_loop, double *conv_precision );
void W_cycle( double *phi, double *rho, int l, int final_level, int *nn, int *level_ind, double *conv_loop, double *conv_precision );	
void up_1step( double *phi, double *rho, int l, int *nn, int *level_ind, double *conv_loop, bool w );
void down_1step( double *phi, double *rho, int l, int *nn, int *level_ind, double *conv_loop, bool w );

#endif
