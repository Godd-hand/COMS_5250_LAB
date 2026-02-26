#include <stdio.h>
#include <math.h>
#include "node.h"
#include "bisector.h"

double f(double x) 
{
    return (x * x * x) - x - 2.0;
}

void run_bisection(double initial_a, double initial_b, double tolerance) 
{
    node* stack = NULL;
    int stack_size;

    if (f(initial_a) * f(initial_b) >= 0) 
    {
        printf("Error: f(a) and f(b) must have opposite signs.\n");
        return;
    }

    Push(initial_a, initial_b, &stack);

    int iteration = 0;
    
    while (1) // Loop until root is found or stack is empty
    {
        GetStackSize(stack, &stack_size);
        if (stack_size == 0)
        {
            break;
        }

        double a, b;
        Pop(&stack, &a, &b);

        double mid = (a + b) / 2.0;
        double f_mid = f(mid);

        printf("Iter %2d: Interval [%.6f, %.6f], Mid = %.6f, f(Mid) = %.6f\n", 
               iteration, a, b, mid, f_mid);

        if (fabs(f_mid) < tolerance || (b - a) / 2.0 < tolerance) 
        {
            printf("\nRoot found at x = %.6f after %d iterations.\n\n", mid, iteration);
            break;
        }

        if (f(a) * f_mid < 0) 
        {
            Push(a, mid, &stack);
        } 
        else 
        {
            Push(mid, b, &stack);
        }

        iteration++;
    }

    DeleteStack(&stack);
}