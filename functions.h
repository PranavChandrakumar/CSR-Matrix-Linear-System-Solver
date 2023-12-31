#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
    double *csr_data;   // Array of non-zero values
    int *col_ind;       // Array of column indices
    int *row_ptr;       // Array of row pointers
    int num_non_zeros;  // Number of non-zero elements
    int num_rows;       // Number of rows in matrix
    int num_cols;       // Number of columns in matrix
} CSRMatrix;


void ReadMMtoCSR(const char *filename, CSRMatrix *matrix);
void spmv_csr(const CSRMatrix *A, const double *x, double *y);

// ###########################################################

int NumDiagonals(const CSRMatrix *A);

void MatrixMultiply(const CSRMatrix *A, const int *all_cols, const double *all_values, const double *x, double *result, int num_diagonals);

void Jacobi(const CSRMatrix *A, double *b, double *x, int num_iterations);
void ComputeResidual(const CSRMatrix *A, double *b, double *Ax);
double ComputeNorm(const CSRMatrix *A, double *b, double *residual, double *Ax);

int CompareElements(const void *a, const void *b);

void save_sparsity_pattern(const CSRMatrix *A, const char *filename);
#endif