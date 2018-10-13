#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    double sum=0;
    for(i=0;i<size*size;i++) {
        sum+=(bref[i]-b[i]);
    }
    if (sum>2)
        return 0;
    else
        return 1;
}
//Doolittle Algorithm
double *my_dgesv(double *A, int n, double *B, int nrhs) {
    double *L = (double*)calloc(n * n, sizeof(double));
    double *U = (double*)calloc(n * n, sizeof(double));
    double *Y = (double*)calloc(n*n, sizeof(double));
    double *X = (double*)calloc(n*n, sizeof(double));
    if ((U == NULL)||(L == NULL)||(X == NULL)||(Y == NULL))
        exit(EXIT_FAILURE);

    int i,j,k;
    for (j = 0; j < n; j++){
        for (i = 0; i <n; i++) {
            if(i<=j){
                U[i*n+j]=A[i*n+j];
                for(k=0;k<=(i-1);k++){
                    U[i*n+j]-=L[i*n+k]*U[k*n+j];
                }
                if (i==j)
                    L[i*n+j]=1;
                else
                    L[i*n+j]=0;
            }
            else{
                L[i*n+j]=A[i*n+j];
                for(k=0; k<=j-1; k++)
                    L[i*n+j]-=L[i*n+k]*U[k*n+j];
                L[i*n+j]/=U[j*n+j];
                U[i*n+j]=0;
            }
        }
    }
 //Conseguidas LU solucionamos o sistema
    for (k=0;k<n;k++){
        for(i=0; i<nrhs; i++){
            Y[i*nrhs+k]=B[i*nrhs+k];
            for(j=0; j<i; j++)
                Y[i*n+k]-=L[i*n+j]*Y[j*nrhs+k];
        }
    }

    for (k=0;k<nrhs;k++){
        for(i=n-1; i>=0; i--){
            X[i * nrhs + k]= Y[i * nrhs + k];
            for(j=i+1; j<n; j++)
                X[i*n+k]-=U[i*n+j]*X[j*n+k];
            X[i*n+k]/=U[i*n+i];
        }
     }
    
    return X;

}

    void main(int argc, char *argv[])
    {
    	
       
        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);
        
        /*
        matrices para o debugeo
        double m2[] = {4, 12,  10,  10,
             12, 6,  7,  4,
             2, 8, 10, 11,
             7, 10, 6, 12};

        double B[]={166,112,161,173,
        100,12,18,13,
        16,11,11,73,
        66,12,61,17};
        */
        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        double *res =my_dgesv(a,n,b,nrhs);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        if (check_result(bref,res,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
    }
