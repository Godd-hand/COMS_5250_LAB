#ifndef GRID_H
#define GRID_H

// Initializes the initial temperature spike
void initialize_grid(double *grid, int nx_padded, int ny, int rank, int size);

// Initializes the spatially variable thermal diffusivity
void initialize_alpha_grid(double *alpha_grid, int nx_padded, int local_nx, int global_nx, int global_ny, int rank);

// Main math loop using the explicit Forward-Time Central-Space (FTCS) scheme
void update_grid(double *grid, double *new_grid, double *alpha_grid, int nx_padded, int ny, double dt, double dx, double dy);

// Writes the computational domain excluding the ghost cells
void write_output(double *grid, int local_nx, int ny, int rank, int step);

#endif