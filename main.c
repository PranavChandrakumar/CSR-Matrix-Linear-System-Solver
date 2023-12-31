#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"



int main(int argc, char *argv[])
{
    // <Handle the inputs here>
    const char *filename = argv[1];
    if (argc != 2){
        printf("Fatal error. Invalid number of arguments. Please input only one .mtx file\n");
        exit(0);
    }

    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    // Initializing all the vector b (in Ax=b)
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    // Set all elements of b to 1
    for (int i = 0; i < A.num_rows; ++i)
    {
        b[i] = 1.0;
    }
    
    double x[A.num_rows];
    for (int i = 0; i < A.num_rows; i++){
        x[i] = 1.0;
    }
    double y[A.num_rows];

    char outputFileName[20] = "output.png";
    
    save_sparsity_pattern(&A,outputFileName);

    clock_t start = clock();

    Jacobi(&A,x,y,10000); //10000 iterations of Jacobi max

    clock_t end = clock();
    double runtime = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Program runtime: %.16f\n", runtime);
    
    
    free(A.row_ptr);
    free(A.col_ind);    
    free(A.csr_data);
}