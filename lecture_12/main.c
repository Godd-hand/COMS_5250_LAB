#include <stdio.h>
#include "bisector.h"

int main() 
{
    double initial_a = 1.0;
    double initial_b = 2.0;
    double tolerance = 1e-6;

    printf("\nBisection Method using Stack \n");
    printf("Function: f(x) = x^3 - x - 2\n");
    printf("Initial interval: [%.4f, %.4f]\n\n", initial_a, initial_b);

    run_bisection(initial_a, initial_b, tolerance);

    return 0;
}