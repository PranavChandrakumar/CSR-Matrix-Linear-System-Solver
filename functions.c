#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <png.h>
#include "functions.h"

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{
    FILE *source;
    source = fopen(filename, "r");

    if (source == NULL)
    {
        printf("Error in opening file\n");
        exit(0);
    }

    char test_line[256];
    int ignore = 0;
    int line_num = 0;
    int ignore_first_line = 0;
    int *all_rows;
    int *all_cols;
    double *all_values;

    while (fgets(test_line, sizeof(test_line), source) != NULL)
    {
        for (int i = 0; i < strlen(test_line); i++) // Checks if the line is a valid line
        {
            if (!isdigit(test_line[i]) && (test_line[i] != '.') && (test_line[i] != ' ') && (test_line[i] != '\n') && (test_line[i] != '-'))
            {
                ignore++; // Adds one to a "checker value" which indicates the current line is a line to be ignored
                break;
            }
        }

        if (ignore == 1) // Condition to check if the line is a line to be ignored
        {
            ignore--;
        }
        else if (ignore == 0 && line_num == 0) // Condition to check if the line is the first line of "important" values (#Rows, #Columns, #Non-zero entries)
        {
            sscanf(test_line, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

            all_rows = (int *)malloc(matrix->num_non_zeros * sizeof(int));
            all_cols = (int *)malloc(matrix->num_non_zeros * sizeof(int));
            all_values = (double *)malloc(matrix->num_non_zeros * sizeof(double));

            if (all_rows == NULL || all_cols == NULL || all_values == NULL)
            {
                printf("Memory allocation failed\n");
                return;
            }
            //printf("%d %d %d\n", matrix->num_rows, matrix->num_cols, matrix->num_non_zeros);
            line_num++;
        }

        if (ignore == 0 && line_num >= 1 && ignore_first_line != 0) // Executes for all lines that contain the non-zero values and their indices
        {
            sscanf(test_line, "%d %d %lf", &all_rows[line_num - 1], &all_cols[line_num - 1], &all_values[line_num - 1]);
            // printf("%d %d %.16f\n", all_rows[line_num - 1], all_cols[line_num - 1], all_values[line_num - 1]);
            line_num++;
        }
        else if (ignore == 0 && line_num >= 1 && ignore_first_line == 0) // Don't assign first line values to actual non-zero values arrays
        {
            ignore_first_line++;
        }
    }

    matrix->row_ptr = (int *)calloc(matrix->num_rows + 1, sizeof(int)); // Allocate +1 to accommodate for num_non_zeros
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));

    if (matrix->row_ptr == NULL || matrix->col_ind == NULL || matrix->csr_data == NULL)
    {
        printf("Memory allocation failed\n");
        return;
    }

    matrix->row_ptr[matrix->num_rows + 1] = matrix->num_non_zeros; // Sets the last element to number of nonzeros

    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        matrix->row_ptr[all_rows[i]]++;
    }

    for (int i = 0; i <= matrix->num_rows; i++)
    {
        matrix->row_ptr[i] += matrix->row_ptr[i - 1];
    }

    int indexes[matrix->num_non_zeros];
    int index_val = 0;
    for (int i = 0; i <= matrix->num_rows; i++)
    {
        for (int j = 0; j < matrix->num_non_zeros; j++)
        {
            if (all_rows[j] == i)
            {
                indexes[index_val] = j;
                index_val++;
            }
        }
    }

    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        matrix->col_ind[i] = all_cols[indexes[i]] - 1;
        matrix->csr_data[i] = all_values[indexes[i]];
    }

    // Print Statements
    // Uncomment these print statements for the same output as the first part of the assignment requires.
    /*
    printf("Number of non-zero entries: %d\n", matrix->num_non_zeros);
    printf("Row Pointers: ");
    for (int i = 0; i <= matrix->num_rows; i++)
    {
        printf("%d ", matrix->row_ptr[i]);
    }
    printf("\n");

    printf("Column Indexes: ");
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        printf("%d ", matrix->col_ind[i]);
    }
    printf("\n");

    printf("Values: ");
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        printf("%.4f ", matrix->csr_data[i]);
    }
    printf("\n\n");
    */
    

    // Print statements required by the sample output for later components of the assignment
    printf("Matrix Name: %s\n", filename);
    printf("Number of non-zero entries: %d\n", matrix->num_non_zeros);
    printf("The dimensions of the matrix: %d by %d\n", matrix->num_rows, matrix->num_cols);


    free(all_rows);
    free(all_cols);
    free(all_values);

    fclose(source);
}

void spmv_csr(const CSRMatrix *A, const double *x, double *y) // Caluclating A*x and writing the result to y. IMPORTANT: This function does not consider if a matrix is symmetrical or not. Therefore, for matricies like b1_ss.mtx, this function is called and only prints out the result of matrix multiplication with x[i] = 1.0.
{
    for (int i = 0; i < A->num_rows; i++)
    { // Initialize all elements of y to 0
        y[i] = 0.0;
    }

    int row_indexes[A->num_non_zeros];
    int row_value = 0;
    int row_index_index = 0;

    for (int i = 1; i < A->num_cols + 1; i++)
    { // Creating a new array which contains all the row indices for all non zero values
        for (int j = 0; j < (A->row_ptr[i] - A->row_ptr[i - 1]); j++)
        {
            row_indexes[row_index_index] = row_value;
            row_index_index++;
        }
        row_value++;
    }

    for (int i = 0; i < A->num_rows; i++)
    { // Vector x HAS to equal the number of columns of A, otherwise we can't do a matrix multiplication
        for (int j = 0; j < A->num_non_zeros; j++)
        {
            if (i == row_indexes[j])
            {
                y[i] += x[A->col_ind[j]] * A->csr_data[j];
            }
        }
    }
    
    printf("Result: \n");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("%f ", y[i]);
    }
    printf("\n");
}

int NumDiagonals(const CSRMatrix *A)
{
    int row_indexes[A->num_non_zeros];
    int row_value = 0;
    int row_index_index = 0;

    for (int i = 1; i < A->num_cols + 1; i++)
    { // Creating a new array which contains all the row indices for all non zero values
        for (int j = 0; j < (A->row_ptr[i] - A->row_ptr[i - 1]); j++)
        {
            row_indexes[row_index_index] = row_value;
            row_index_index++;
        }
        row_value++;
    }

    int num_diagonal_entries = 0;
    for (int i = 0; i < A->num_non_zeros; i++)
    {
        if (row_indexes[i] == A->col_ind[i])
        {
            num_diagonal_entries++;
        }
    }
    return num_diagonal_entries;
}

typedef struct { //Structure which holds the all the information about a matrix element, which is useful for soring algorithm in Symmetrify, as element-wise relationships between row, column and value are retained
    int row;
    int col;
    double value;
} MatrixElement;

int CompareElements(const void *a, const void *b) { //Function which is used to compare elements in qsort of Symmetrify
    return ((MatrixElement*)a)->row - ((MatrixElement*)b)->row;
}

void Symmetrify(const CSRMatrix *A, int *symmetrified_rows, int *symmetrified_col, double *symmetrified_values, const int num_diagonal_entries)
{
    int *row_indexes = (int *)calloc(A->num_non_zeros, sizeof(int));
    int row_value = 0;
    int row_index_index = 0;

    for (int i = 1; i < A->num_cols + 1; i++)
    { // Creating a new array which contains all the row indices for all non zero values
        for (int j = 0; j < (A->row_ptr[i] - A->row_ptr[i - 1]); j++)
        {
            row_indexes[row_index_index] = row_value;
            row_index_index++;
        }
        row_value++;
    }

    int pos = 0;//This is the index, relative to the end of the list prior to reflection, that all reflected values will have

    for (int i = 0; i < A->num_non_zeros; i++) //
    {
        symmetrified_rows[i] = row_indexes[i];
        symmetrified_col[i] = A->col_ind[i];
        symmetrified_values[i] = A->csr_data[i];
        if (row_indexes[i] != A->col_ind[i])
        {
            symmetrified_rows[A->num_non_zeros + pos] = A->col_ind[i]; //Adds the row and column indices of the symmetrified values to the list of all indices
            symmetrified_col[A->num_non_zeros + pos] = row_indexes[i];
            symmetrified_values[A->num_non_zeros + pos] = A->csr_data[i];
            pos++;
        }
    }

    MatrixElement *elements = malloc((2 * A->num_non_zeros - num_diagonal_entries) * sizeof(MatrixElement)); //Allocating memory for the matrix elements

    for (int i = 0; i < 2 * A->num_non_zeros - num_diagonal_entries; i++) //Assigning values to the members of the struct (elements)
    {
        elements[i].row = symmetrified_rows[i];
        elements[i].col = symmetrified_col[i];
        elements[i].value = symmetrified_values[i];
    }

    qsort(elements, 2 * A->num_non_zeros - num_diagonal_entries, sizeof(MatrixElement), CompareElements); // Sort elements using qsort, it's so much faster than the bubble sort algorithm that I previously used

    for (int i = 0; i < 2 * A->num_non_zeros - num_diagonal_entries; i++)//Updating input arrays with sorted entries
    {
        symmetrified_rows[i] = elements[i].row;
        symmetrified_col[i] = elements[i].col;
        symmetrified_values[i] = elements[i].value;
    }

    memset(A->row_ptr, 0, (A->num_rows+1) * sizeof(int)); //Resetting all elemenets or row_ptr to 0

    for (int i = 0; i < 2 * A->num_non_zeros - num_diagonal_entries; i++) // Writing arranged row, column, and values list to output lists
    {
        A->row_ptr[symmetrified_rows[i]+1]++;
    }
    for (int i = 1; i < A->num_rows + 1; i++)
    {
        A->row_ptr[i] += A->row_ptr[i-1];
    }

    /*
    for (int i = 0; i < A->num_rows + 1; i++)
    {
        printf("%d ", A->row_ptr[i]);
    }
    printf("\n");

    for (int i = 0; i < 2 * A->num_non_zeros - num_diagonal_entries; i++) // Writing arranged row, column, and values list to output lists
    {
        printf("%d ", symmetrified_col[i]);
    }
    printf("\n");
    */

   free(row_indexes);
}

void MatrixMultiply(const CSRMatrix *A, const int *all_cols, const double *all_values, const double *x, double *result, int num_diagonals)
{
    // Initialize result vector to zero
    for (int i = 0; i < A->num_rows; i++)
    {
        result[i] = 0.0;
    }

    // Perform matrix-vector multiplication
    for (int i = 0; i < A->num_rows; i++)
    {
        for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
        {
            result[i] += x[all_cols[j]]*all_values[j];
        }
    }
}

void Jacobi(const CSRMatrix *A, double *b, double *x, int num_iterations)
{
    
    int num_diagonal_entries = NumDiagonals(A);
    if (num_diagonal_entries != A->num_cols)
    {
        printf("Not all diagonal entries are non-zero. Jacobi approach will fail.\n");
        for (int i = 0; i < A->num_cols; i++)
        { // Initial guess for x is all 1 for simple matrix multiplication.
            x[i] = 1.0;
        }
        spmv_csr(A, x, b);
        return;
    }
    
    int *AllRows = (int *)calloc(2 * A->num_non_zeros - num_diagonal_entries, sizeof(int));
    int *AllCols = (int *)calloc(2 * A->num_non_zeros - num_diagonal_entries, sizeof(int));
    double *AllVals = (double *)calloc(2 * A->num_non_zeros - num_diagonal_entries, sizeof(double));

    printf("Valid lower triangular matrix, symmetrifying values...\n");
    Symmetrify(A, AllRows, AllCols, AllVals, num_diagonal_entries); //Reflect all non-zeros across diagonal, and then with new matrix, multiply 
    printf("Symmetrification complete, commencing Jacobi iterations\n");
    double diagonal_entries[num_diagonal_entries]; // A list of all diagonal entries, used for algebra.
    int index_of_diagonal_entry = 0;               // Need this for program logic

    for (int i = 0; i < 2 * A->num_non_zeros - num_diagonal_entries; i++)
    {
        if (AllRows[i] == AllCols[i]) //Makes algebraic process of Jacobi method substantially easier if there's a list of all diagonal values
        {
            diagonal_entries[index_of_diagonal_entry] = AllVals[i];
            index_of_diagonal_entry++;
        }
    }
    
    /*
    for (int i = 0; i < 2 * num_diagonal_entries; i++)
    {
        printf("%.16f\n", diagonal_entries[i]);
    }
    */

    double *result = (double *)malloc(A->num_cols * sizeof(double)); // Result vector of A*x from Jacobi method
    double *residual = (double *)malloc(A->num_cols * sizeof(double));
   
    double norm_diff = 1.0; //Initialize norm_diff

    for (int i = 0; i < A->num_cols; i++)
    { // Initial guess for x is all zeros.
        x[i] = 0.0;
    }
    for (int i = 0; i < num_iterations; i++)
    {
        if (norm_diff < 1e-7) //Condition to check if norm difference is sufficient to stop iterations
        {
            break;
        }
        MatrixMultiply(A, AllCols, AllVals, x, result, num_diagonal_entries);
        for (int i = 0; i < A->num_cols; i++)
        {
            x[i] = (b[i] - result[i] + diagonal_entries[i] * x[i]) / diagonal_entries[i]; //Algebra, this is based off of the numerical process for the Jacobi method
            //printf("%f ",x[i]);
        }
        norm_diff = ComputeNorm(A, b,residual,result);
        //printf("%f\n", norm_diff);
    }

    // This prints out the components of the solution vector
    /*
    printf("Result:\n");
    for (int i = 0; i < A->num_cols; i++)
    {
        printf("%f ", x[i]);
    }
    printf("\n");
    */

    printf("Residual Norm: %.16f\n", norm_diff);

    free(AllRows);
    free(AllCols);
    free(AllVals);
    free(residual);
    free(result);
}

void ComputeResidual(const CSRMatrix *A, double *b, double *Ax)
{
    double residual[A->num_cols];
    printf("Residual: [");
    for (int i = 0; i < A->num_cols; i++)
    {
        residual[i] = b[i] - Ax[i];
        printf("%f, ", residual[i]);
    }
    printf("]\n");
}

double ComputeNorm(const CSRMatrix *A, double *b, double *residual, double *Ax)
{
    double norm = 0.0;
    for (int i = 0; i < A->num_cols; i++)
    {
        residual[i] = b[i] - Ax[i];
        norm += residual[i] * residual[i];
    }
    norm = sqrt(norm);
    return norm;
}

void save_sparsity_pattern(const CSRMatrix *A, const char *filename) //Function to create .png visualization, code is provided by the instructor.
{
    int width = A->num_cols * 10;  // Adjust dimensions for larger image
    int height = A->num_rows * 10; // Adjust dimensions for larger image
    png_bytep *row_pointers = (png_bytep *)malloc(height * sizeof(png_bytep));
    if (row_pointers == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }
    for (int i = 0; i < height; i++)
    {
        row_pointers[i] = (png_byte *)malloc(width * sizeof(png_byte));
        if (row_pointers[i] == NULL)
        {
            fprintf(stderr, "Memory allocation failed\n");
            return;
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            row_pointers[i][j] = 255; // Initialize to white (255 = white in grayscale)
            // Set a larger block for non-zero elements
            if (i % 10 == 0 && j % 10 == 0)
            {
                for (int k = 0; k < 10; k++)
                {
                    for (int l = 0; l < 10; l++)
                    {
                        if ((i + k) < height && (j + l) < width)
                        {
                            row_pointers[i + k][j + l] = 255; // Adjust the block color if needed
                        }
                    }
                }
            }
            // Check if the element at (i, j) is non-zero
            for (int k = A->row_ptr[i / 10]; k < A->row_ptr[(i / 10) + 1]; k++)
            {
                if (A->col_ind[k] == (j / 10))
                {
                    row_pointers[i][j] = 0; // Set to black
                    break;
                }
            }
        }
    }
    FILE *fp = fopen(filename, "wb");
    if (!fp)
    {
        fprintf(stderr, "Error opening file for writing\n");
        return;
    }
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr)
    {
        fprintf(stderr, "Error creating PNG write struct\n");
        fclose(fp);
        return;
    }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        fprintf(stderr, "Error creating PNG info struct\n");
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(fp);
        return;
    }
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
    png_init_io(png_ptr, fp);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    for (int i = 0; i < height; i++)
    {
        free(row_pointers[i]);
    }
    free(row_pointers);
}