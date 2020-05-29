#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
extern const float L;
extern float dx, bc;
extern const int N;
extern double *analytic,*potential,*density;

int init_sin(float a,float b,float c){
	analytic = (double *)malloc( N * N * sizeof(double) );
	potential = (double *)malloc( N * N * sizeof(double) );
	density = (double *)malloc( N * N * sizeof(double) );
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			int index = i*N + j;
			analytic[index] = sin(a*float(i)*dx)*sin(b*float(j)*dx)+c;
			density[index]   = -(a*a+b*b)*(analytic[index]-c);
			if( i==0 || j==0 || i==(N-1) || j==(N-1)){potential[index]=bc;}
                        else{potential[index]=0;}
		}
	}
	return 0;
}
