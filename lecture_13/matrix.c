/*
File containing the matrix functions 
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"

// Creates a new matrix and initializes its values to 0.0
matrix new_matrix(const int rows, const int cols) {
    matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    assert(rows > 0);
    assert(cols > 0);
    mat.val = (double*)malloc(sizeof(double) * rows * cols);
    for (int i = 0; i < (rows * cols); i++) {
        mat.val[i] = 0.0;
    }
    return mat;
}

// Generates an n x n identity matrix
matrix matrix_identity(int n) {
    matrix I = new_matrix(n, n);
    for (int i = 1; i <= n; i++) {
        mget(I, i, i) = 1.0;
    }
    return I;
}

// Multiplies a matrix by a scalar value
matrix matrix_scalar_mult(const matrix* A, double scalar) {
    matrix C = new_matrix(A->rows, A->cols);
    for (int i = 1; i <= A->rows; i++) {
        for (int j = 1; j <= A->cols; j++) {
            mget(C, i, j) = mgetp(A, i, j) * scalar;
        }
    }
    return C;
}

// Evaluates the transpose of a matrix
matrix matrix_transpose(const matrix* A) {
    matrix C = new_matrix(A->cols, A->rows);
    for (int i = 1; i <= A->rows; i++) {
        for (int j = 1; j <= A->cols; j++) {
            mget(C, j, i) = mgetp(A, i, j);
        }
    }
    return C;
}

// Adds two matrices together
matrix matrix_add(const matrix* A, const matrix* B) {
    const int rows = A->rows;
    const int cols = A->cols;
    assert(rows == B->rows);
    assert(cols == B->cols);
    matrix C = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            mget(C, i, j) = mgetp(A, i, j) + mgetp(B, i, j);
        }
    }
    return C;
}

// Subtracts matrix B from matrix A
matrix matrix_sub(const matrix* A, const matrix* B) {
    const int rows = A->rows;
    const int cols = A->cols;
    assert(rows == B->rows);
    assert(cols == B->cols);
    matrix C = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            mget(C, i, j) = mgetp(A, i, j) - mgetp(B, i, j);
        }
    }
    return C;
}

// Computes matrix multiplication (A * B)
matrix matrix_mult(const matrix* A, const matrix* B) {
    const int rowsA = A->rows;
    const int colsA = A->cols;
    const int rowsB = B->rows;
    const int colsB = B->cols;
    assert(colsA == rowsB);
    matrix C = new_matrix(rowsA, colsB);
    for (int i = 1; i <= rowsA; i++) {
        for (int j = 1; j <= colsB; j++) {
            for (int k = 1; k <= colsA; k++) {
                mget(C, i, j) += mgetp(A, i, k) * mgetp(B, k, j);
            }
        }
    }
    return C;
}

// Computes element-by-element product
matrix matrix_dot_mult(const matrix* A, const matrix* B) {
    const int rows = A->rows;
    const int cols = A->cols;
    assert(rows == B->rows);
    assert(cols == B->cols);
    matrix C = new_matrix(rows, cols);
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            mget(C, i, j) = mgetp(A, i, j) * mgetp(B, i, j);
        }
    }
    return C;
}

// Prints a matrix
void print_matrix_full(const matrix* mat, char* varname) {
    assert(mat->rows > 0);
    assert(mat->cols > 0);
    printf("\n %.100s =\n", &varname[1]);
    for(int i = 1; i <= mat->rows; i++) {
        printf(" | ");
        for(int j = 1; j <= mat->cols; j++) {
            printf("%10.3e", mgetp(mat, i, j));
            if (j < mat->cols) { printf(", "); }
            else { printf(" "); }
        }
        printf("|\n");
    }
    printf("\n");
}

// Creates a new vector
vector new_vector(const int size) {
    vector vec;
    vec.size = size;
    assert(size > 0);
    vec.val = (double*)malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++) {
        vec.val[i] = 0.0;
    }
    return vec;
}

// Adds two vectors
vector vector_add(const vector* x, const vector* y) {
    const int size = x->size;
    assert(size == y->size);
    vector z = new_vector(size);
    for (int i = 1; i <= size; i++) {
        vget(z, i) = vgetp(x, i) + vgetp(y, i);
    }
    return z;
}

// Subtracts vector y from vector x
vector vector_sub(const vector* x, const vector* y) {
    const int size = x->size;
    assert(size == y->size);
    vector z = new_vector(size);
    for (int i = 1; i <= size; i++) {
        vget(z, i) = vgetp(x, i) - vgetp(y, i);
    }
    return z;
}

// Calculates the dot product of two vectors
double vector_dot_mult(const vector* x, const vector* y) {
    const int size = x->size;
    assert(size == y->size);
    double z = 0.0;
    for (int i = 1; i <= size; i++) {
        z += vgetp(x, i) * vgetp(y, i);
    }
    return z;
}

// Computes the L2 norm of a vector
double vector_norm(const vector* x) {
    double sum = 0.0;
    for (int i = 1; i <= x->size; i++) {
        sum += vgetp(x, i) * vgetp(x, i);
    }
    return sqrt(sum);
}

// Multiplies a matrix by a vector
vector matrix_vector_mult(const matrix* A, const vector* x) {
    const int rows = A->rows;
    const int cols = A->cols;
    const int size = x->size;
    assert(cols == size);
    vector Ax = new_vector(rows);
    for (int i = 1; i <= rows; i++) {
        double tmp = 0.0;
        for (int j = 1; j <= size; j++) {
            tmp += mgetp(A, i, j) * vgetp(x, j);
        }
        vget(Ax, i) = tmp;
    }
    return Ax;
}

// Prints a vector
void print_vector_full(const vector* vec, char* varname) {
    assert(vec->size > 0);
    printf("\n %.100s =\n", &varname[1]);
    printf(" | ");
    for(int i = 1; i <= vec->size; i++) {
        printf("%10.3e", vgetp(vec, i));
        if (i < vec->size) { printf(", "); }
    }
    printf(" |^T\n\n");
}

// Linear solve via Gaussian elimination with partial pivoting
vector solve(const matrix* A, const vector* b) {
    const int rows = A->rows;
    const int cols = A->cols;
    const int size = b->size;
    assert(rows == cols);
    assert(rows == size);
    
    // Create copies so we don't modify the original matrices/vectors
    matrix Acopy = new_matrix(rows, cols);
    vector bcopy = new_vector(size);
    for (int i = 1; i <= rows; i++) {
        vget(bcopy, i) = vgetp(b, i);
        for (int j = 1; j <= cols; j++) {
            mget(Acopy, i, j) = mgetp(A, i, j);
        }
    }
    
    vector x = new_vector(rows);
    for (int i = 1; i <= (size - 1); i++) {
        int p = i;
        double maxA = -100.0e0;
        for (int j = i; j <= size; j++) {
            double tmp = fabs(mget(Acopy, j, i));
            if (tmp > maxA) {
                p = j;
                maxA = tmp;
            }
        }
        
        if (maxA <= 1.0e-14) {
            printf(" Cannot invert system\n");
            exit(1);
        }
        
        if (p != i) {
            for (int j = 1; j <= size; j++) {
                double tmp1 = mget(Acopy, i, j);
                mget(Acopy, i, j) = mget(Acopy, p, j);
                mget(Acopy, p, j) = tmp1;
            }
            double tmp2 = vget(bcopy, i);
            vget(bcopy, i) = vget(bcopy, p);
            vget(bcopy, p) = tmp2;
        }
        
        for (int j = (i + 1); j <= size; j++) {
            double dm = mget(Acopy, j, i) / mget(Acopy, i, i);
            for (int k = 1; k <= size; k++) {
                mget(Acopy, j, k) = mget(Acopy, j, k) - dm * mget(Acopy, i, k);
            }
            vget(bcopy, j) = vget(bcopy, j) - dm * vget(bcopy, i);
        }
    }
    
    vget(x, size) = vget(bcopy, size) / mget(Acopy, size, size);
    for (int j = 1; j <= (size - 1); j++) {
        double sum = 0.0e0;
        for (int k = (size - j + 1); k <= size; k++) {
            sum = sum + mget(Acopy, size - j, k) * vget(x, k);
        }
        vget(x, size - j) = (vget(bcopy, size - j) - sum) / mget(Acopy, size - j, size - j);
    }
    
    free_matrix(&Acopy);
    free_vector(&bcopy);
    return x;
}

// Memory freeing
void free_matrix(matrix* mat) {
    if (mat->val != NULL) {
        free(mat->val);
        mat->val = NULL;
    }
}

void free_vector(vector* vec) {
    if (vec->val != NULL) {
        free(vec->val);
        vec->val = NULL;
    }
}