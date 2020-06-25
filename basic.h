#ifndef BASIC_H
#define BASIC_H

#define OPENMP
#define DEBUG
#define PARALLEL_FOR //Set the parallel method of relaxation.c
#define FULL_WEIGHTING
#define TEST2


int ind( int i, int j, int n);
void print( double *matric, int n);
void add_correction( double *phi_old, double *phi_corr, int n );
void fill_zero( double *phi_guess, int n  );
void test_prol_rest( const int N );

//	define the calculation function type for exact solvers
typedef void (*cal_fn)( double*, double*, int, double*, float, bool );

#endif
