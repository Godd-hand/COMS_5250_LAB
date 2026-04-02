#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "adaptive_quad.h"

// Simpson's rule
double Q_simpson(double a, double b, double (*f)(double, double), double param) 
{
    double c = 0.5 * (a + b);
    return (b - a) / 6.0 * (f(a, param) + 4.0 * f(c, param) + f(b, param));
}

// recursive adaptive quadrature - serial
double AdaptiveIntSerial(double a, double b, double TOL, double (*f)(double, double), double param) 
{
    double Qab = Q_simpson(a, b, f, param);
    
    double c = 0.5 * (a + b);
    double Qac = Q_simpson(a, c, f, param);
    double Qcb = Q_simpson(c, b, f, param);

    // Error estimation is given as 1/15
    double error_est = (1.0 / 15.0) * fabs(Qac + Qcb - Qab);

    if (error_est < TOL) 
    {
        return Qac + Qcb;
    } 
    else 
    {
        return AdaptiveIntSerial(a, c, 0.5 * TOL, f, param) + 
               AdaptiveIntSerial(c, b, 0.5 * TOL, f, param);
    }
}

// OpenMP Version 1: Dividing the domain across the threads 
double AdaptiveInt_OMP_v1(double a, double b, double TOL, double (*f)(double, double), double param, int thread_count) 
{
    double total_I = 0.0;
    double width = (b - a) / (double)thread_count;
    double TOL_local = TOL / (double)thread_count;

    #pragma omp parallel num_threads(thread_count)
    {
        int my_rank = 0;
        
        #ifdef _OPENMP
        my_rank = omp_get_thread_num();
        #endif

        double a_local = a + my_rank * width;
        double b_local = a_local + width;

        // Perform serial adaptive integration on each sub-domain
        double I_thread = AdaptiveIntSerial(a_local, b_local, TOL_local, f, param);

        #pragma omp critical
        total_I += I_thread;
    }
    
    return total_I;
}