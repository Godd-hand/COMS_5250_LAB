#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NX 1152
#define NY 1152

void initialize_grid(double *grid) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            grid[i * NY + j] = 20.0; 
        }
    }
    
    int block_size = 100;
    int start_x = (NX / 2) - (block_size / 2);
    int end_x = (NX / 2) + (block_size / 2);
    int start_y = (NY / 2) - (block_size / 2);
    int end_y = (NY / 2) + (block_size / 2);

    for (int i = start_x; i < end_x; i++) {
        for (int j = start_y; j < end_y; j++) {
            grid[i * NY + j] = 5000.0; 
        }
    }
}

void initialize_alpha_grid(double *alpha_grid) {
    int cx = NX / 2;
    int cy = NY / 2;
    double target_radius = 150.0;
    
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            double dist = sqrt(pow(i - cx, 2) + pow(j - cy, 2));
            int idx = i * NY + j;
            if (dist < target_radius) {
                alpha_grid[idx] = 0.02; // Insulative target zone
            } else {
                alpha_grid[idx] = 0.1;  // Conductive ambient Zone
            }
        }
    }
}

void update_grid(double *grid, double *new_grid, double *alpha_grid, double dt, double dx, double dy) {
    for (int i = 1; i < NX - 1; i++) {
        for (int j = 1; j < NY - 1; j++) {
            int idx = i * NY + j;
            
            double local_alpha = alpha_grid[idx];
            double cx = (local_alpha * dt) / (dx * dx);
            double cy = (local_alpha * dt) / (dy * dy);

            new_grid[idx] = grid[idx] + 
                            cx * (grid[(i + 1) * NY + j] - 2.0 * grid[idx] + grid[(i - 1) * NY + j]) + 
                            cy * (grid[i * NY + (j + 1)] - 2.0 * grid[idx] + grid[i * NY + (j - 1)]);
        }
    }
}

void write_output(double *grid, int step) {
    char filename[50];
    sprintf(filename, "output_rank0_step%06d.dat", step); 
    FILE *fp = fopen(filename, "wb");
    if (fp != NULL) {
        fwrite(grid, sizeof(double), NX * NY, fp);
        fclose(fp);
    }
}

int main() {
    double *grid = (double *)malloc(NX * NY * sizeof(double));
    double *new_grid = (double *)malloc(NX * NY * sizeof(double));
    double *alpha_grid = (double *)malloc(NX * NY * sizeof(double));

    initialize_grid(grid);
    initialize_alpha_grid(alpha_grid);

    double dx = 1.0, dy = 1.0, dt = 1.0; 
    int num_steps = 200000, output_interval = 500;

    double compute_time = 0.0;
    clock_t start_time, end_time;

    for (int step = 0; step <= num_steps; step++) {
        start_time = clock();
        update_grid(grid, new_grid, alpha_grid, dt, dx, dy);
        
        double *temp = grid;
        grid = new_grid;
        new_grid = temp;
        end_time = clock();

        compute_time += (double)(end_time - start_time) / CLOCKS_PER_SEC;

        if (step % output_interval == 0) {
            write_output(grid, step);
        }
    }

    printf("Serial pure compute time: %f seconds.\n", compute_time);

    free(grid); free(new_grid); free(alpha_grid);
    return 0;
}