/* Lab assignment 19: Signal Normalization
 * Normalize a given signal of length N in 2-norm (energy of the signal) and
 * max norm (peak of the signal).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// function to re-initialize the vector with random values
void init_vector(double* x, int N) 
{
    for (int i = 0; i < N; i++) 
    {
        x[i] = ((double)rand() / RAND_MAX) * 100.0 - 50.0;
    }
}

int main(int argc, char* argv[]) 
{
    int N = 100000000; 
    int thread_counts[] = {1, 2, 4, 8, 16};
    
    printf("===================================================\n");
    printf(" Vector Normalization (N = %d)\n", N);
    printf("===================================================\n\n");

    // Allocate memory 
    double* x = (double*)malloc(N * sizeof(double));
    if (x == NULL) 
    {
        printf("Memory allocation failed!\n");
        return 1;
    }

    srand(time(NULL));

    // Loop through the different thread counts
    for (int t = 0; t < 5; t++) 
    {
        int num_threads = thread_counts[t];
        printf(" Running with %d Thread(s) \n", num_threads);

        double time1, time2;

        // 1. Fine-Grained 2-Norm
        
        init_vector(x, N);
        double norm2_fine = 0.0;
        
        time1 = omp_get_wtime();
        
        #pragma omp parallel for num_threads(num_threads) reduction(+:norm2_fine)
        for (int i = 0; i < N; i++) 
        {
            norm2_fine += x[i] * x[i];
        }
        norm2_fine = sqrt(norm2_fine);
        
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < N; i++) 
        {
            x[i] = x[i] / norm2_fine;
        }
        
        time2 = omp_get_wtime();
        printf(" Fine-Grained 2-Norm Time:   %.5f sec\n", time2 - time1);

        // 2. Coarse-Grained 2-Norm

        init_vector(x, N);
        double norm2_coarse = 0.0;
        
        time1 = omp_get_wtime();
        
        #pragma omp parallel num_threads(num_threads)
        {
            int my_rank = omp_get_thread_num();
            int N_per_thread = N / num_threads;
            int istart = my_rank * N_per_thread;
            // Handle uneven division for the last thread
            int iend = (my_rank == num_threads - 1) ? N : (my_rank + 1) * N_per_thread;

            double norm_thread = 0.0;
            for (int i = istart; i < iend; i++) 
            {
                norm_thread += x[i] * x[i];
            }

            #pragma omp critical
            norm2_coarse += norm_thread;

            #pragma omp barrier // Wait for all threads to add to the global sum

            // Only one thread needed for square root
            #pragma omp single
            {
                norm2_coarse = sqrt(norm2_coarse);
            } 

            // Normalize the chunk
            for (int i = istart; i < iend; i++) 
            {
                x[i] = x[i] / norm2_coarse;
            }
        }
        
        time2 = omp_get_wtime();
        printf(" Coarse-Grained 2-Norm Time: %.5f sec\n", time2 - time1);

        // 3. Fine-Grained Max-Norm
        
        init_vector(x, N);
        double norm_max_fine = 0.0;
        
        time1 = omp_get_wtime();
        
        #pragma omp parallel for num_threads(num_threads) reduction(max:norm_max_fine)
        for (int i = 0; i < N; i++) 
        {
            double val = fabs(x[i]);
            if (val > norm_max_fine) 
            {
                norm_max_fine = val;
            }
        }
        
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < N; i++) 
        {
            x[i] = x[i] / norm_max_fine;
        }
        
        time2 = omp_get_wtime();
        printf(" Fine-Grained Max-Norm Time:   %.5f sec\n", time2 - time1);

        // 4. Coarse-Grained Max-Norm
    
        init_vector(x, N);
        double norm_max_coarse = 0.0;
        
        time1 = omp_get_wtime();
        
        #pragma omp parallel num_threads(num_threads)
        {
            int my_rank = omp_get_thread_num();
            int N_per_thread = N / num_threads;
            int istart = my_rank * N_per_thread;
            int iend = (my_rank == num_threads - 1) ? N : (my_rank + 1) * N_per_thread; // uneven division for the last thread

            double max_thread = 0.0;
            for (int i = istart; i < iend; i++) 
            {
                double val = fabs(x[i]);
                if (val > max_thread) 
                {
                    max_thread = val;
                }
            }

            #pragma omp critical
            {
                if (max_thread > norm_max_coarse) 
                {
                    norm_max_coarse = max_thread;
                }
            }

            #pragma omp barrier // barrier 

            for (int i = istart; i < iend; i++) 
            {
                x[i] = x[i] / norm_max_coarse;
            }
        }
        
        time2 = omp_get_wtime();
        printf(" Coarse-Grained Max-Norm Time: %.5f sec\n\n", time2 - time1);
    }

    free(x);
    return 0;
}