/*Lab assignment 18*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Specifying PI since M_PI is troublesome sometimes
const double PI = 3.14159265358979323846;

// Phi(x) using erf() from math.h 
double Phi(double x) 
{
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

// Black-Scholes Exact Formula
double BlackScholes(double S0, double K, double r, double sigma, double T_mat) 
{
    double d1 = (log(S0 / K) + (r + (sigma * sigma) / 2.0) * T_mat) / (sigma * sqrt(T_mat));
    double d2 = d1 - sigma * sqrt(T_mat);
    
    return S0 * Phi(d1) - K * exp(-r * T_mat) * Phi(d2);
}

/* The Simpson's rule integrand for the European call option price integral
*/

// The integrand: (S - K) * p(S)
double integrand(double S, double K, double S0, double r, double sigma, double T_mat) 
{
    // By examination, we see that if S = K, C is max, 
    // but since we integrate from K to Smax, S-K has to be >= 0
    double mu = log(S0) + (r - (sigma * sigma) / 2.0) * T_mat;
    
    // Probability density function p(S)
    double pS = (1.0 / (S * sigma * sqrt(2.0 * PI * T_mat))) * exp(-(pow(log(S) - mu, 2.0)) / (2.0 * sigma * sigma * T_mat));
                
    return (S - K) * pS;
}

// Composite Simpson's Rule using OpenMP Parallel #4 (Reduction)
double CompSimpson(double a, double b, int N, int thread_count, 
                   double K, double S0, double r, double sigma, double T_mat) 
{
    double h = (b - a) / (double)N;
    
    // Evaluate endpoints
    double T = integrand(a, K, S0, r, sigma, T_mat) + integrand(b, K, S0, r, sigma, T_mat);

    // Parallelize the loop using reduction on T
    #pragma omp parallel for num_threads(thread_count) reduction(+: T)
    for (int i = 1; i < N; i++) 
    {
        double x = a + i * h;
        // Even indices multiply by 2, odd indices multiply by 4
        double weight = (i % 2 == 0) ? 2.0 : 4.0;
        T += weight * integrand(x, K, S0, r, sigma, T_mat);
    }
    
    return (h / 3.0) * T;
}

int main() 
{
    // Test case parameters
    double S0 = 100.0;
    double K = 100.0;
    double r = 0.05;
    double T_mat = 1.0;
    double sigma = 0.2;
    
    // Numerical parameters
    int N = 10000000; 
    if (N % 2 != 0) N++; // Simpson's rule requires an even N
    
    double Smax = 5.0 * S0; // Truncate at 5*S0
    
    printf("=========================================================\n");
    printf(" European Call Option Pricing (Black-Scholes vs Simpson's)\n");
    printf("=========================================================\n");
    printf(" S0 = %.2f, K = %.2f, r = %.2f, T = %.2f, sigma = %.2f\n", S0, K, r, T_mat, sigma);
    printf(" Grid points (N) = %d\n", N);
    printf(" Integration Limits: [%.2f, %.2f]\n\n", K, Smax);
    
    // 1. Calculate Exact Price
    double exact_price = BlackScholes(S0, K, r, sigma, T_mat);
    printf(" Exact Black-Scholes Price: %.8f\n", exact_price);
    printf("---------------------------------------------------------\n");

    // 2. Calculate Numerical Price using OpenMP and observe thread scaling
    int thread_counts[] = {1, 2, 4, 8, 16};
    
    for (int t = 0; t < 5; t++) 
    {
        int threads = thread_counts[t];
        
        double time1 = omp_get_wtime();
        
        // Evaluate integral
        double integral = CompSimpson(K, Smax, N, threads, K, S0, r, sigma, T_mat);
        
        // Option price 
        double numerical_price = exp(-r * T_mat) * integral;
        
        double time2 = omp_get_wtime();
        
        double err = fabs(numerical_price - exact_price);
        
        printf(" Threads: %2d | Price: %.8f | Err: %10.5e | Time: %.5f sec\n", 
               threads, numerical_price, err, time2 - time1);
    }
    printf("=========================================================\n");

    return 0;
}