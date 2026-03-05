/*
Header file for matrix class definition and implementation.

*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

// Define the matrix struct
typedef struct matrix matrix;
struct matrix {
    int rows;
    int cols;
    double* val;
};

// Define the vector struct
typedef struct vector vector;
struct vector {
    int size;
    double* val;
};

// Shortcut evaluate macros for easier element access
#define mget(mat,i,j) mat.val[(i-1)*mat.cols+(j-1)]
#define mgetp(mat,i,j) mat->val[(i-1)*mat->cols+(j-1)]
#define vget(vec, i) vec.val[(i-1)]
#define vgetp(vec, i) vec->val[(i-1)]
#define print_matrix(mat) print_matrix_full(mat, #mat);
#define print_vector(vec) print_vector_full(vec, #vec);
#define print_Scalar(z) print_Scalar_full(z, #z);

// Function declarations
matrix new_matrix(const int rows, const int cols);
void print_matrix_full(const matrix* mat, char* varname);
matrix matrix_add(const matrix* A, const matrix* B);
matrix matrix_sub(const matrix* A, const matrix* B);
matrix matrix_mult(const matrix* A, const matrix* B);
matrix matrix_dot_mult(const matrix* A, const matrix* B);
matrix matrix_transpose(const matrix* A);
matrix matrix_scalar_mult(const matrix* A, double scalar);
matrix matrix_identity(int n);

// Vector functions
vector new_vector(const int size);
void print_vector_full(const vector* vec, char* varname);
vector vector_add(const vector* x, const vector* y);
vector vector_sub(const vector* x, const vector* y);
double vector_dot_mult(const vector* x, const vector* y);
vector matrix_vector_mult(const matrix* A, const vector* x);
double vector_norm(const vector* x);

vector solve(const matrix* A, const vector* b);

// Memory management
void free_matrix(matrix* mat);
void free_vector(vector* vec);

// Eigenvalue algorithms
double power_iteration(const matrix* A, vector* v, double tol, int max_iters);
double shifted_inverse_iteration(const matrix* A, vector* v, double mu, double tol, int max_iters);
double rayleigh_quotient_iteration(const matrix* A, vector* v, double tol, int max_iters);

#endif