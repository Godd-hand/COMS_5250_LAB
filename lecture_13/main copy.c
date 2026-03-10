#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

// Helper to reset the initial guess vector to 1.0s
void reset_guess(vector* v) {
    for (int i = 1; i <= v->size; i++) {
        vgetp(v, i) = 1.0;
    }
}

int main() {
    int n_values[] = {5, 10};
    double tol = 1e-8;
    int max_iters = 1000;

    for (int idx = 0; idx < 2; idx++) {
        int n = n_values[idx];
        printf("==========================================\n");
        printf("Evaluating Matrix A for n = %d\n", n);
        printf("==========================================\n");

        // Generate Blur Matrix A
        matrix A = new_matrix(n, n);
        for (int r = 1; r <= n; r++) {
            mget(A, r, r) = 2.0 / 4.0;
            if (r > 1) mget(A, r, r - 1) = 1.0 / 4.0;
            if (r < n) mget(A, r, r + 1) = 1.0 / 4.0;
        }
        
        vector v0 = new_vector(n);
        double lambda;

        // Power Iteration
        reset_guess(&v0);
        lambda = power_iteration(&A, &v0, tol, max_iters);
        printf("Power Iteration:\n");
        printf("  Dominant Eigenvalue: %.6f\n\n", lambda);

        // Shifted Inverse Iteration (Shift mu = 0.0 to find the smallest eigenvalue)
        reset_guess(&v0);
        double mu = 0.0;
        lambda = shifted_inverse_iteration(&A, &v0, mu, tol, max_iters);
        printf("Shifted Inverse Iteration (mu = %.1f):\n", mu);
        printf("  Eigenvalue closest to shift: %.6f\n\n", lambda);

        // Rayleigh Quotient Iteration
        reset_guess(&v0);
        lambda = rayleigh_quotient_iteration(&A, &v0, tol, max_iters);
        printf("Rayleigh Quotient Iteration:\n");
        printf("  Eigenvalue found: %.6f\n\n", lambda);

        free_matrix(&A);
        free_vector(&v0);
    }

    return 0;
}