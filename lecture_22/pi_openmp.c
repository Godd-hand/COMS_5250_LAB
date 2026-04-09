#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char* argv[]) 
{
    if (argc != 2) 
    {
        printf("Usage: ./pi_openmp.exe <N>\n"); 
        return 1;
    }

    long long N = strtoll(argv[1], NULL, 10);
    int thread_counts[] = {2, 4, 8};
    
    FILE* fp = fopen("openmp_pi_results.txt", "w");
    if (!fp) 
    {
        printf("Error opening output file.\n");
        return 1;
    }

    printf("==========================================\n");
    printf(" OpenMP Monte Carlo Pi (N = %lld)\n", N);
    printf("==========================================\n\n");

    for (int t = 0; t < 3; t++) 
    {
        int num_threads = thread_counts[t];
        omp_set_num_threads(num_threads);
        
        long long global_count = 0;
        double start_time = omp_get_wtime();

        // Coarse-grain parallelism
        #pragma omp parallel
        {
            int my_rank = omp_get_thread_num();
            unsigned int seed = 12345 + my_rank; 
            
            long long N_per_thread = N / num_threads;
            long long local_count = 0;

            for (long long i = 0; i < N_per_thread; i++) 
            {
                // Generate random numbers between 0.0 and 1.0
                double x = (double)rand_r(&seed) / RAND_MAX;
                double y = (double)rand_r(&seed) / RAND_MAX;

                if (x * x + y * y <= 1.0) 
                {
                    local_count++;
                }
            }

            #pragma omp critical
            global_count += local_count;
        }

        double end_time = omp_get_wtime();
        
        // Ensuring accurate division if N is not perfectly divisible by threads
        long long total_points = (N / num_threads) * num_threads;
        double pi_estimate = 4.0 * (double)global_count / (double)total_points;
        double walltime = end_time - start_time;

        printf(" Threads: %2d | Pi: %.9f | Time: %f sec\n", num_threads, pi_estimate, walltime);
        fprintf(fp, "%d %.9f %f\n", num_threads, pi_estimate, walltime);
    }

    fclose(fp);
    return 0;
}