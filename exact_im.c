//input residual and size of the matrix, output phi^tilda^corr_gridsize. By method of inverse matrix.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"

extern const float L;
extern float dx, bc;
extern const int N;
extern double *residual;

int exact_im( double *matrix,int n ){
    //turn into Au=b where u and b are tempararily vectors
    //construct A
    double *A;
    const int n_ext = pow(n-2,2); // abbr for n_extension
    A = (double *)malloc(n_ext*n_ext*sizeof(double));
    for(int i=0;i<n_ext-1;i++)
	for(int j=0;j<n_ext-1;j++){
        if (i==j) A[ind(i, j, n_ext)]=4;
        else if (i==j+1||i==j-1||i==j+3||i==j-3) A[ind(i, j, n_ext)]=-1;
        else A[ind(i, j, n_ext)]=0;
    }
    //construct b
    double *b;
    b = (double *)malloc(n_ext*sizeof(double));
    for(int i=1;i<n-1;i++)
	for(int j=1;j<n-1;j++){
        b[ind(i-1 ,j-1, n)]=pow((N/n)*dx,2)*matrix[ind(j, i, n)];
        if(i==1 && j==1) b[ind(i-1 ,j-1, n)] += matrix[ind(0, 1, n)]+matrix[ind(1, 0, n)];
        else if(i==2 && j==1) b[ind(i-1 ,j-1, n)] += matrix[ind(2, 0, n)];
        else if(i==3 && j==1) b[ind(i-1 ,j-1, n)] += matrix[ind(4, 1, n)]+matrix[ind(3, 0, n)];
        else if(i==1 && j==2) b[ind(i-1 ,j-1, n)] += matrix[ind(0, 2, n)];
        else if(i==3 && j==2) b[ind(i-1 ,j-1, n)] += matrix[ind(4, 2, n)];
        else if(i==1 && j==3) b[ind(i-1 ,j-1, n)] += matrix[ind(0, 3, n)]+matrix[ind(1, 4, n)];
        else if(i==2 && j==3) b[ind(i-1 ,j-1, n)] += matrix[ind(2, 4, n)];
        else if(i==3 && j==3) b[ind(i-1 ,j-1, n)] += matrix[ind(4, 3, n)]+matrix[ind(3, 4, n)];
    }
    //calculate A^(-1)


    
    //calculate u=A^(-1)b


    //turn u into a matrix
    

}