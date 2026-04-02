/*
* Main.c for lab assignment 20. 
* Adaptive quadrature for 1D electrostatic potential and Bessel function.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adaptive_quad.h"
#ifdef _OPENMP
#include <omp.h>
#endif


const double PI = 3.14159265358979323846;

// Integrand for 1D Electrostatic Potential
double V_integrand(double x, double x0) 
{
    // Avoiding division by zero
    if (fabs(x - x0) < 1e-14) 
    {
        return 0.0;
    }
    return (exp(-x * x) - exp(-x0 * x0)) / fabs(x - x0);
}

// Integrand for Bessel function
double J_integrand(double theta, double x) 
{
    return cos(x * sin(theta));
}

int main(int argc, char* argv[]) 
{
    // test parameters
    double TOL = 1e-6;
    int thread_counts[] = {1, 2, 4, 8, 16};
    
    int num_pts_V = 100;
    int num_pts_J = 200;

    printf("==========================================\n");
    printf(" Adaptive Quadrature (Version 1)\n");
    printf(" TOL: %e\n", TOL);
    printf("==========================================\n\n");

    // Updated: Timing 1D Electrostatic Potential
    
    printf("Timing Electrostatic Potential ...\n");
    
    for (int t = 0; t < 5; t++) 
    {
        int threads = thread_counts[t];
        double start_time = omp_get_wtime();

        for (int i = 1; i < num_pts_V; i++) 
        {
            double x0 = -1.0 + i * (2.0 / num_pts_V);
            
            // Integrate the first part 
            double integral = AdaptiveInt_OMP_v1(-1.0, 1.0, TOL, V_integrand, x0, threads);
            
            // Adding the analytical part
            double analytical_part = exp(-x0 * x0) * log((1.0 - x0) / (x0 - (-1.0)));
            double V_x0 = integral + analytical_part;
            
            if (V_x0 == 0.0) continue;
        }

        double end_time = omp_get_wtime();
        printf(" Threads: %2d | Walltime: %f sec\n", threads, end_time - start_time);
    }

    // Timing Bessel Function

    printf("\nTiming Bessel function ...\n");
    
    // x belongs to [0, 50]
    for (int t = 0; t < 5; t++) 
    {
        int threads = thread_counts[t];
        double start_time = omp_get_wtime();

        for (int i = 0; i <= num_pts_J; i++) 
        {
            double x = i * (50.0 / num_pts_J);
            
            // Integrate from 0 to pi using the OpenMP adaptive quadrature
            double integral = AdaptiveInt_OMP_v1(0.0, PI, TOL, J_integrand, x, threads);
            
            double J_x = (1.0 / PI) * integral;
            
            if (J_x == 0.0) continue;
        }

        double end_time = omp_get_wtime();
        printf(" Threads: %2d | Walltime: %f sec\n", threads, end_time - start_time);
    }

    // Updated: Generating Data Files outside of the timing
    
    printf("\nGenerating output files ...\n");

    FILE* fp_v = fopen("V_out.dat", "w");
    if (!fp_v) 
    {
        printf("Error opening V_out.dat\n"); return 1; 
    }

    for (int i = 1; i < num_pts_V; i++) 
    {
        double x0 = -1.0 + i * (2.0 / num_pts_V);
        double integral = AdaptiveInt_OMP_v1(-1.0, 1.0, TOL, V_integrand, x0, 16);
        double analytical_part = exp(-x0 * x0) * log((1.0 - x0) / (x0 - (-1.0)));
        double V_x0 = integral + analytical_part;
        
        fprintf(fp_v, "%.6f %.10f\n", x0, V_x0);
    }
    fclose(fp_v);

    FILE* fp_j = fopen("J_out.dat", "w");
    if (!fp_j) 
    { 
        printf("Error opening J_out.dat\n"); return 1; 
    }

    for (int i = 0; i <= num_pts_J; i++) 
    {
        double x = i * (50.0 / num_pts_J);
        double integral = AdaptiveInt_OMP_v1(0.0, PI, TOL, J_integrand, x, 16);
        double J_x = (1.0 / PI) * integral;
        
        fprintf(fp_j, "%.6f %.10f\n", x, J_x);
    }
    fclose(fp_j);

    printf("\nComputations finished. Run plot.py to visualize results.\n");
    
    return 0;
}