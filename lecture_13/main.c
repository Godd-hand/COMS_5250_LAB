/*
Combined execution file for testing the matrix and noise functions in the context of a deblurring problem. 

With steps:
1. Generate a symmetric blur kernel matrix A.
2. Generate a true signal vector x (a step function).
3. Compute the observed blurred signal b = Ax.
4. Add normal distribution noise to b to create b_noisy.
5. Solve the linear system Ax_rec = b_noisy directly and compute the error.
6. Implement Tikhonov regularization to solve (A^T A + lambda I)x_rec_tik = A^T b_noisy and compute the error for different lambda values.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "noise.h"

int main() {
    srand(time(NULL));
    
    int n_values[] = {64, 128};
    double sigmas[] = {1e-6, 1e-4, 1e-2};
    double lambdas[] = {1e-6, 1e-4, 1e-2};

    for (int i = 0; i < 2; i++) {
        int n = n_values[i];
        printf("==========================================\n");
        printf("Testing with n = %d\n", n);
        printf("==========================================\n");

        // 1. Generate the symmetric blur kernel matrix A
        matrix A = new_matrix(n, n);
        for (int r = 1; r <= n; r++) {
            mget(A, r, r) = 2.0 / 4.0;
            if (r > 1) mget(A, r, r - 1) = 1.0 / 4.0;
            if (r < n) mget(A, r, r + 1) = 1.0 / 4.0;
        }

        // 2. Generate the true signal vector x
        vector x = new_vector(n);
        for (int r = 1; r <= n; r++) {
            if (r >= n / 4 && r <= n / 2) {
                vget(x, r) = 1.0;
            } else {
                vget(x, r) = 0.0;
            }
        }
        
        double norm_x = vector_norm(&x);

        // 3. Compute observed blurred signal b = Ax
        vector b = matrix_vector_mult(&A, &x);

        // Calculate A transpose
        matrix AT = matrix_transpose(&A);
        matrix ATA = matrix_mult(&AT, &A);
        matrix I = matrix_identity(n);

        for (int j = 0; j < 3; j++) {
            double sigma = sigmas[j];
            printf("\nSigma: %e \n", sigma);

            // Add noise to b
            vector b_noisy = new_vector(n);
            for (int r = 1; r <= n; r++) {
                vget(b_noisy, r) = vget(b, r) + generate_normal(sigma);
            }

            // Solve Ax_rec = b_noisy directly
            vector x_rec = solve(&A, &b_noisy);
            
            // Calculate Error: ||x - x_rec|| / ||x||
            vector diff = vector_sub(&x, &x_rec);
            double err = vector_norm(&diff) / norm_x;
            printf("Direct Solve Error: %e\n", err);
            
            free_vector(&x_rec);
            free_vector(&diff);

            // Tikhonov Regularization: (A^T A + lambda I)x = A^T b
            vector ATb = matrix_vector_mult(&AT, &b_noisy);

            for (int k = 0; k < 3; k++) {
                double lambda = lambdas[k];
                
                // Form lambda * I
                matrix lam_I = matrix_scalar_mult(&I, lambda);
                
                // Form system matrix M = A^T A + lambda I
                matrix M = matrix_add(&ATA, &lam_I);
                
                // Solve M * x_rec_tik = ATb
                vector x_rec_tik = solve(&M, &ATb);
                
                // Calculate Tikhonov Error
                vector diff_tik = vector_sub(&x, &x_rec_tik);
                double err_tik = vector_norm(&diff_tik) / norm_x;
                printf("Tikhonov Solve Error (lambda = %e): %e\n", lambda, err_tik);

                free_matrix(&lam_I);
                free_matrix(&M);
                free_vector(&x_rec_tik);
                free_vector(&diff_tik);
            }

            free_vector(&b_noisy);
            free_vector(&ATb);
        }

        // Clean up matrices and vectors for current n
        free_matrix(&A);
        free_matrix(&AT);
        free_matrix(&ATA);
        free_matrix(&I);
        free_vector(&x);
        free_vector(&b);
    }

    return 0;
}