#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
extern const float L;
extern float dx, bc;
extern const int N;
extern double *analytic,*potential,*density;

int ind( int i,int j, int NGrid ) {
	return i * NGrid + j;
}


int print(double *matrix,int n){
        for( int i=0;i<n;i++ ){
                for( int j=0;j<n;j++ ){
                        int index = i*N+j;
                        printf("%.3f\t",matrix[index]);
                }
                printf("\n");
        }
	return 0;
}

