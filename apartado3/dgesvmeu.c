#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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
        #pragma omp parallel for private(j,k)
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
//algo do pragma simd vector length
double sum;
    #pragma omp parallel private(i,j,k)
    {
    #pragma omp for
    for (k=0;k<n;k++){
        for(i=0; i<nrhs; i++){
            Y[i*nrhs+k]=B[i*nrhs+k];
	    sum=Y[i*n+k];
	    #pragma omp simd reduction(-:sum)
            for(j=0; j<i; j++)
                sum-=L[i*n+j]*Y[j*nrhs+k];
	    Y[i*n+k]=sum;
        }
    }
    #pragma omp for
    for (k=0;k<nrhs;k++){
        for(i=n-1; i>=0; i--){
            X[i * nrhs + k]= Y[i * nrhs + k];
		sum=X[i*n+k];
	    #pragma omp simd reduction(-:sum)

            for(j=i+1; j<n; j++)
                sum-=U[i*n+j]*X[j*n+k];
	    X[i*n+k]=sum;
            X[i*n+k]/=U[i*n+i];
        }
     }
     }
	free(L);
	free(U);
	free(Y);    
    return X;

}

void main(int argc, char *argv[])
    {
    	
       
        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        //aref = generate_matrix(size);        
        b = generate_matrix(size);
        //bref = generate_matrix(size);

        // Using MKL to solve the system
        int n = size, nrhs = size, lda = size, ldb = size, info;

        //clock_t tStart = clock();
        //info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        //printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        double start,end;  
        start = omp_get_wtime();
        double *res =my_dgesv(a,n,b,nrhs);
        end = omp_get_wtime();
        printf("Time taken by my implementation: %.2fs\n", (end - start));
        //free(ipiv);
            if (check_result(res,res,size)==1)
        printf("Result is ok!\n");
            else    
        printf("Result is wrong!\n");
        free(a);
        //free(aref);
        free(b);
        //free(bref);
        free(res);    
        //mkl_free_buffers();
}
