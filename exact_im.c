//input residual and size of the matrix, output phi^tilda^corr_gridsize. By method of inverse matrix.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#define PI acos(-1)
#include "basic.h"


//find comatrix
//delete row and column the new matrix is a (m-1)*(m-1) matrix
void comatrix(double *matrix, int p, int q, int m, double *comx){
    for (int i=0;i<m;i++){
        for (int j=0;j<m;j++){    
            if(i<p && j<q) comx[ind(i,j,m-1)]=matrix[ind(i,j,m)];
            else if(i<p && j>q) comx[ind(i,j-1,m-1)]=matrix[ind(i,j,m)];
            else if(i>p && j<q) comx[ind(i-1,j,m-1)]=matrix[ind(i,j,m)];
            else if(i>p && j>q) comx[ind(i-1,j-1,m-1)]=matrix[ind(i,j,m)];
        }
    }   
}

//find determinant
//by method of expanding out the first row recursively until the comatrix is 1 dimenison
double det(double *matrix, int m){
    int D=0;
    int sign=1;
    if (m==1) return matrix[ind(0,0,m)];
    double* comx;
    comx = (double*)malloc((m-1) * (m-1) * sizeof(double));
    for (int f=0;f<m;f++){
        comatrix(matrix, 0, f, m, comx);
        D += sign*matrix[ind(0,f,m)]*det(comx,m-1);
        sign =-sign;
    }
    free(comx);
    return D;
}

//find inverse matrix
//if C is the matrix of cofactors (ie every element is the determinant of its comatrix)
//then the inverse matrix = (-1)^(i+j)*(transpose of C)/(determinant of matrix)
void inversematrix(double *matrix, int n, double *invmx){
    int sign=1;
    double* comx;
    comx = (double*)malloc((n-1) * (n-1) * sizeof(double));
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            sign=pow(-1, i+j);
            comatrix(matrix, j, i, n, comx);
            invmx[ind(i,j,n)]=sign*det(comx,n-1)/det(matrix,n);         
        }
    }
    free(comx);
}

//main function
extern const int N;
extern const double dx;
void exact_im( double *phi_corr, double *residual, int n, double *conv_criterion, float omega, bool w ) {
    //convert Laplacian(phi_corr)=-residual to Au=b 
    //where u and b are vectors of dimension = (n-2)^2 
    const int n_ext = pow(n - 2, 2);
    
    //construct A
    double* A;
    A = (double*)malloc(n_ext * n_ext * sizeof(double));
    for (int i = 0; i < n_ext; i++) {
        for (int j = 0; j < n_ext; j++) {
            if (i == j) A[ind(i, j, n_ext)] = 4;
            else if (i == j + 1 || i == j - 1 || i == j + 3 || i == j - 3) A[ind(i, j, n_ext)] = -1;
            else A[ind(i, j, n_ext)] = 0;
        }
    }
   
    //construct b
    double* b;
    b = (double*)malloc((n-2) * (n-2) * sizeof(double));
    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < n - 2; j++) {
            b[ind(i, j, n - 2)] = pow((N / n) * dx, 2) * residual[ind(j + 1, i + 1, n)];
        }
    }
    b[ind(0, 0, n_ext)] += residual[ind(0, 1, n)] + residual[ind(1, 0, n)];
    b[ind(0, 1, n_ext)] += residual[ind(2, 0, n)];
    b[ind(0, 2, n_ext)] += residual[ind(4, 1, n)] + residual[ind(3, 0, n)];
    b[ind(0, 3, n_ext)] += residual[ind(0, 2, n)];
    b[ind(0, n_ext-4, n_ext)] += residual[ind(4, 2, n)];
    b[ind(0, n_ext-3, n_ext)] += residual[ind(0, 3, n)] + residual[ind(1, 4, n)];
    b[ind(0, n_ext-2, n_ext)] += residual[ind(2, 4, n)];
    b[ind(0, n_ext-1, n_ext)] += residual[ind(4, 3, n)] + residual[ind(3, 4, n)];

    //calculate A^(-1)
    double* invmx_A;
    invmx_A = (double*)malloc(n_ext * n_ext * sizeof(double));
    inversematrix(A, n_ext, invmx_A);

    //calculate u=A^(-1)b
    double *phi_corr_nb; //nb stands for "no boundary"
    phi_corr_nb = (double *)malloc(n_ext*sizeof(double));
    for(int i=0;i<n_ext;i++){
	    for(int j=0;j<n_ext;j++){
            phi_corr_nb[ind(0,i,n_ext)]+=invmx_A[ind(i,j,n_ext)]*b[ind(0,j,n_ext)];
        }
    }
    //phi_corr_nb is now a matrix phi_corr_nb[ind(i,j,n-2)], without boundary
    
    //paste the boundary
    for(int i=0;i<n;i++){
	    for(int j=0;j<n;j++){
            if( i==0 || j==0 || i==(n-1) || j==(n-1)) 
                phi_corr[ind(i,j,n)]=0;
            else 
                phi_corr[ind(i,j,n)]=phi_corr_nb[ind(i-1,j-1,n-2)];
        }
    }
    free(A);
    free(b);
    free(invmx_A);
    free(phi_corr_nb);
    printf("Exact solver by inverse matrix terminated\n");
}
