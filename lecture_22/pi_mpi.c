#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]) 
{
    int comm_sz;
    int my_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (argc != 2) 
    {
        if (my_rank == 0) printf("Usage: mpirun -np <cores> ./pi_mpi.exe <N>\n");
        MPI_Finalize();
        return 1;
    }

    long long N = strtoll(argv[1], NULL, 10);
    long long N_local = N / comm_sz;
    long long local_count = 0;
    long long global_count = 0;

    unsigned int seed = 12345 + my_rank; // Unique seed for each process

    double start_time;
    if (my_rank == 0) 
    {
        start_time = MPI_Wtime();
    }

    // Coarse-grain
    for (long long i = 0; i < N_local; i++) 
    {
        double x = (double)rand_r(&seed) / RAND_MAX;
        double y = (double)rand_r(&seed) / RAND_MAX;

        if (x * x + y * y <= 1.0) 
        {
            local_count++;
        }
    }

    // Combining all local counts into the global count
    MPI_Reduce(&local_count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0) 
    {
        double end_time = MPI_Wtime();
        
        long long total_points = N_local * comm_sz;
        double pi_estimate = 4.0 * (double)global_count / (double)total_points;
        double walltime = end_time - start_time;
        
        printf(" MPI Cores: %2d | Pi: %.9f | Time: %f sec\n", comm_sz, pi_estimate, walltime);
        
        // Append to a single file since the MPI executable is called multiple times by the bash script
        FILE* fp = fopen("mpi_pi_results.txt", "a");
        if(fp) 
        {
            fprintf(fp, "%d %.9f %f\n", comm_sz, pi_estimate, walltime);
            fclose(fp);
        }
    }

    MPI_Finalize();
    return 0;
}