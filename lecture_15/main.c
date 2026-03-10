#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "trimatrix.h"

void Hessenberg(const matrix* A, trimatrix* T);
void QRA(trimatrix* T);

int main() {
    // Array holding the sizes to be tested 
    int sizes[] = {6, 12};

    // Loop through N = 6, 12
    for (int idx = 0; idx < 2; idx++) {
        int N = sizes[idx];
        printf("================================\n");
        printf("QR Algorithm for N = %d\n", N);
        printf("================================\n");

        // 1. Generate matrix A
        matrix A = new_matrix(N, N);
        for (int r = 1; r <= N; r++) {
            // Set the center diagonal
            mget(A, r, r) = 2.0 / 4.0;
            
            // Set the sub-diagonal
            if (r > 1) mget(A, r, r - 1) = 1.0 / 4.0;
            
            // Set the super-diagonal
            if (r < N) mget(A, r, r + 1) = 1.0 / 4.0;
        }

        printf("Original Matrix A:");
        print_matrix(&A);

        // Hessenberg reduction to tridiagonal form
        trimatrix T = new_trimatrix(N);
        Hessenberg(&A, &T);

        printf("Reduction to Tridiagonal Form:");
        print_trimatrix(&T);

        // QR Algorithm to find eigenvalues of T
        
        QRA(&T);

        printf("After QR Algorithm:");
        print_trimatrix(&T);

        // Clean up memory before the next loop
        delete_matrix(&A);
        delete_trimatrix(&T);
    }

    return 0;
}