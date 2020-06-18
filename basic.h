#ifndef BASIC_H
#define BASIC_H

int ind( int i, int j, int n);
void print( double *matric, int n);
void add_correction( double *phi_old, double *phi_corr, int n );
void test_prol_rest( const int N );

typedef void (*cal_fn)( double*, double*, int, double*, float, bool );
void choose_solver( cal_fn solver_name, double *, double *, int, double *, float, bool );

#endif
