/*
 * Image Blurring exercise with OpenMP
 * This program generates a random N x N image and applies a 3x3 averaging filter to blur it.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main() 
{
    // Image dimension of N x N
    // N=1000 will be tried first
    int N = 1000; 
    
    // Using 1D arrays to represent 2D matrices so that the processor can access the memory
    // more efficiently across the rows.
    double* img = (double*)malloc(N * N * sizeof(double));
    double* blur = (double*)malloc((N - 2) * (N - 2) * sizeof(double));

    // Generating random image
    srand(time(NULL));
    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            img[i * N + j] = rand() % 256;
        }
    }

    // Save the original image to a file to compare with blurred
    printf("Saving original image...\n");
    FILE* f_orig = fopen("original.dat", "w");
    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            fprintf(f_orig, "%.1f ", img[i * N + j]);
        }
        fprintf(f_orig, "\n");
    }
    fclose(f_orig);

    printf("Image size: %d x %d\n", N, N);
    printf("Starting blur...\n\n");

    int thread_counts[] = {1, 2, 4, 8};

    // Run the blurring process for 1, 2, 4, and 8 threads
    for (int t = 0; t < 4; t++) 
    {
        int num_threads = thread_counts[t];
        
        double start_time = omp_get_wtime();

        // Parallelize the outer loop
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 1; i < N - 1; i++) 
        {
            for (int j = 1; j < N - 1; j++) 
            {
                double sum = 0.0;
                
                // 3x3 averaging filter
                for (int di = -1; di <= 1; di++) 
                {
                    for (int dj = -1; dj <= 1; dj++) 
                    {
                        sum += img[(i + di) * N + (j + dj)];
                    }
                }
                
                // Store in the smaller (N-2) x (N-2) blurred image
                blur[(i - 1) * (N - 2) + (j - 1)] = sum / 9.0;
            }
        }

        double end_time = omp_get_wtime();
        printf(" Threads: %d | Walltime: %f sec\n", num_threads, end_time - start_time);
    }

    // Save the blurred image to a file for Python visualization
    printf("\nSaving blurred image data...\n");
    FILE* f_blur = fopen("blurred.dat", "w");
    for (int i = 0; i < N - 2; i++) 
    {
        for (int j = 0; j < N - 2; j++) 
        {
            fprintf(f_blur, "%.1f ", blur[i * (N - 2) + j]);
        }
        fprintf(f_blur, "\n");
    }
    fclose(f_blur);

    free(img);
    free(blur);

    // Spawn the Python visualization script
    printf("Spawning Python visualization script...\n");
    system("/home/george/miniconda3/envs/agentic/bin/python plot.py");

    return 0;
}