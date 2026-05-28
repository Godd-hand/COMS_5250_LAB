#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"

// Set up the initial 5000C heat spike in the center

void initialize_grid(double *grid, int nx_padded, int ny, int rank, int size) {
    int local_nx = nx_padded - 2;
    int global_nx = local_nx * size; 
    int global_ny = ny;

    int block_size = 100;
    int start_x = (global_nx / 2) - (block_size / 2);
    int end_x = (global_nx / 2) + (block_size / 2);
    int start_y = (global_ny / 2) - (block_size / 2);
    int end_y = (global_ny / 2) + (block_size / 2);

    // Set background to room temp (Dirichlet BC)
    for (int i = 0; i < nx_padded; i++) {
        for (int j = 0; j < ny; j++) {
            grid[i * ny + j] = 20.0;
        }
    }

    // Map global coordinates to figure out which rank owns the center block
    for (int i = 1; i <= local_nx; i++) {
        int global_x = rank * local_nx + (i - 1);
        if (global_x >= start_x && global_x < end_x) {
            for (int j = start_y; j < end_y; j++) {
                grid[i * ny + j] = 5000.0;
            }
        }
    }
}

// Creates the spatially variable diffusivity map
// Inner radius traps heat (0.02), outer area conducts normally (0.1)

void initialize_alpha_grid(double *alpha_grid, int nx_padded, int local_nx, int global_nx, int global_ny, int rank) {
    int global_cx = global_nx / 2;
    int global_cy = global_ny / 2;
    double target_radius = 150.0; // Radius of the insulative target zone
    
    for (int i = 0; i < nx_padded; i++) {
        for (int j = 0; j < global_ny; j++) {
            // Map the local MPI array index to the true physical coordinate
            int global_x = rank * local_nx + (i - 1); 
            int global_y = j;
            
            // Calculate radial distance from the center of the domain
            double dist = sqrt(pow(global_x - global_cx, 2) + pow(global_y - global_cy, 2));
            
            int idx = i * global_ny + j;
            if (dist < target_radius) {
                alpha_grid[idx] = 0.02; // insulative target zone
            } else {
                alpha_grid[idx] = 0.1;  // ambient zone with normal thermal conductivity
            }
        }
    }
}

// 2D Explicit FTCS scheme

void update_grid(double *grid, double *new_grid, double *alpha_grid, int nx_padded, int ny, double dt, double dx, double dy) {
    // Loop ignores the outer ghost cell rows (i=0 and i=nx_padded-1)
    for (int i = 1; i < nx_padded - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            int idx = i * ny + j;
            
            // Calculate face-averaged thermal diffusivity between adjacent control volumes
            double alpha_right = 0.5 * (alpha_grid[(i + 1) * ny + j] + alpha_grid[idx]);
            double alpha_left  = 0.5 * (alpha_grid[(i - 1) * ny + j] + alpha_grid[idx]);
            double alpha_up    = 0.5 * (alpha_grid[i * ny + (j + 1)] + alpha_grid[idx]);
            double alpha_down  = 0.5 * (alpha_grid[i * ny + (j - 1)] + alpha_grid[idx]);

            // Apply the conservative explicit stencil
            new_grid[idx] = grid[idx] + 
                (dt / (dx * dx)) * (alpha_right * (grid[(i + 1) * ny + j] - grid[idx]) - alpha_left * (grid[idx] - grid[(i - 1) * ny + j])) + 
                (dt / (dy * dy)) * (alpha_up * (grid[i * ny + (j + 1)] - grid[idx]) - alpha_down * (grid[idx] - grid[i * ny + (j - 1)]));
        }
    }
}

// Dump raw binary data so python can stitch it later
void write_output(double *grid, int local_nx, int ny, int rank, int step) {
    char filename[50];
    sprintf(filename, "output_rank%d_step%06d.dat", rank, step);
    FILE *fp = fopen(filename, "wb");
    if (fp != NULL) {
        fwrite(grid, sizeof(double), local_nx * ny, fp);
        fclose(fp);
    }
}