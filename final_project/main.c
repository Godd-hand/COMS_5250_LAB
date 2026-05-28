#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "grid.h"

int main(int argc, char *argv[]) {
    int rank, size;
    
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 1D Cartesian Domain Decomposition
    // Note: 1152 divides perfectly by all the core counts (2, 4, 8, 16, 32 and 64)  
    int global_nx = 1152, global_ny = 1152;
    int local_nx = global_nx / size;
    int ny = global_ny;
    int nx_padded = local_nx + 2; // Allocate 2 extra rows for the ghost Cells

    // Allocate memory for spatial domains and composite material grid
    double *grid = (double *)malloc(nx_padded * ny * sizeof(double));
    double *new_grid = (double *)malloc(nx_padded * ny * sizeof(double));
    double *alpha_grid = (double *)malloc(nx_padded * ny * sizeof(double));

    // Apply initial conditions and material properties
    initialize_grid(grid, nx_padded, ny, rank, size);
    initialize_alpha_grid(alpha_grid, nx_padded, local_nx, global_nx, global_ny, rank);

    // Discretization parameters
    double dx = 1.0, dy = 1.0;
    double dt = 1.0; 
    int num_steps = 200000;
    int output_interval = 500;

    double compute_time = 0.0;
    double start_time, end_time;

    // Main time-stepping loop
    for (int step = 0; step <= num_steps; step++) {
        start_time = MPI_Wtime();

        MPI_Request reqs[4];
        int req_count = 0;
        
        // Ghost cell exchange (non-blocking to prevent deadlocks)
         
        // Send top row, receive from rank above
        if (rank > 0) {
            MPI_Isend(&grid[1 * ny], ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(&grid[0 * ny], ny, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &reqs[req_count++]);
        }
        // Send bottom row, receive from rank below
        if (rank < size - 1) {
            MPI_Isend(&grid[local_nx * ny], ny, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &reqs[req_count++]);
            MPI_Irecv(&grid[(local_nx + 1) * ny], ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &reqs[req_count++]);
        }
        
        // Wait for ghost cells to arrive before doing math updates
        MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
        
        // Execute numerical scheme
        update_grid(grid, new_grid, alpha_grid, nx_padded, ny, dt, dx, dy);

        // Pointer swap instead of copying data
        double *temp = grid; grid = new_grid; new_grid = temp;

        end_time = MPI_Wtime();
        compute_time += (end_time - start_time);

        // Save data to disk (outside the timer to avoid I/O bottlenecks in scaling metrics)
        if (step % output_interval == 0) {
            write_output(&grid[1 * ny], local_nx, ny, rank, step);
        }
    }

    if (rank == 0) {
        printf("Pure compute time: %f seconds using %d ranks.\n", compute_time, size);
    }

    free(grid); 
    free(new_grid);
    free(alpha_grid);
    MPI_Finalize();
    return 0;
}