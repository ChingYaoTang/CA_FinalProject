#ifndef BASIC_H
#define BASIC_H

//#define OPENMP
//#define DEBUG
//#define PARALLEL_GPU  //Set the parallel method of relaxation.c
#define WO_OMP
#define FULL_WEIGHTING// Set method of restriction
#define TEST2         //Set the initial problem

#define GPU

int ind( int i, int j, int n);
void print( double *matric, int n);
void add_correction( double *phi_old, double *phi_corr, int n );
void fill_zero( double *phi_guess, int n  );
void test_prol_rest( const int N ,double *phi_corr_h); 
//	define the calculation function type for exact solvers
typedef void (*cal_fn)( double*, double*, int, double*, float, bool );



#endif
