/* 
Lab Assignment 7

Develop a script "main program" with the name lab.c: 
Computing factorial of an integer n, exponential of a real number x, and logarithm of a real number y
*/

# include <stdio.h>
# include <math.h>

// Function to compute the factorial of an integer n

int main()
{
    // Declare variables
    const int n = 8; // Example integer for factorial
    int nfactorial = 1;
    int i;
    // Compute factorial using a loop
    for (i = 1; i <= n; i++)
    {
        nfactorial *= i;
    }
    printf("Factorial of %d is %d\n", n, nfactorial);

    // Declare variables for exponential and logarithm
    const double x = 2.0; // for exponential
    const double y = 10.0; // for logarithm

    printf("Exponential of %f is %f\n", x, exp(x));
    printf("Logarithm of %f is %f\n", y, log(y));

    return 0;
}

/*
    long long factorial(int n) 
{
    if (n < 0) 
    {
        printf("Error: Factorial is not defined for negative numbers.\n");
        return 0;
    }
    
    long long result = 1;
    
    for (int i = 1; i <= n; i++) 
    {
        result *= i;
    }
    return result;
}

int main()
{
    // Factorial
    int n = 8; // integer for factorial
    long long nfactorial = factorial(n);

    // Exponential
    double x = 2.0; // real number for exponential
    double exp_result = exp(x);

    // Logarithm
    double y = 10.0; // real number for logarithm
    double log_result = log(y);

    // Print results
    printf("Factorial of %d is %lld\n", n, nfactorial);
    printf("Exponential of %f is %f\n", x, exp_result);
    printf("Logarithm of %f is %f\n", y, log_result);

    return 0;
}
*/
