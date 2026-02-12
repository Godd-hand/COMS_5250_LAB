/*

Develop source codes ("main program") with the name lab8.c: define two functions 
(1) to compute the factorial of an integer n, and 
(2) to compute the exponential of a real number x.

    a. For computing factorial(n), define a function with recursive calls.
    b. For computing exponential(x), define a function with Taylor expansions (using the python examples),
        where the constant e = 2.718281828459, x_0 = round(x) is the closest integer to x and use the pow(.,.) from math.h
        for computing powers. Use your own factorial function.
    c. Generate a set of points x = 0, 0.02, 0.04,...,0.98,1 and compute the exponential of these points using your function.
    d. Output the data to a data file, and visualize the exp(x) function with the data in Python using (system("python3 plot.py")).

*/

# include <stdio.h>
# include <math.h>
# include <stdlib.h>

// Function to compute the factorial of an integer n using recursion
long long factorial(int n)
{
    if (n == 0 || n == 1) 
    {
        return 1;
    }
    
    return n * factorial(n - 1);
}

// Function to compute the exponential of a real number x using Taylor expansion
double exponential(double x)
{
     if (x == 0) 
    {
        return 1.0; // e^0 is 1
    }
    const double e = 2.718281828459; 
    int x0 = round(x);
    double sum = 1.0; // Start with the first term of the series (x^0 / 0!)
    double term = 1.0; // To store each term of the series
    int n = 1; // Start with the first term (x^1 / 1!) since we already added the x^0 / 0! term to the sum

    // using the created factorial function to compute the terms of the series

    while (fabs(term) > 1e-14) // specifying tolerance
    {
        term = pow(x, n) / factorial(n); 
        sum += term; 
        n++; 
    }

    return sum;
}

int main()
{
    // Generating the points and computing the exponential values.
    // The points are generated from 0 to 1 with a step of 0.02.
    double start = 0.0;
    double end = 1.0;
    double step = 0.02;

    for (double x = start; x <= end; x += step) 
    {
        double exp_value = exponential(x);
        printf("%f %f\n", x, exp_value); // Output the data to the console
    }

    // Writing the data to a file
    FILE *file = fopen("exp_data.txt", "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file for writing.\n");
        return 1;
    }
    for (double x = start; x <= end; x += step) 
    {
        double exp_value = exponential(x);
        fprintf(file, "%f %f\n", x, exp_value); // Write the data to the file
    }
    fclose(file);

    // Visualize the exp(x) function with the data in Python
    system("/home/george/miniconda3/envs/agentic/bin/python3 plot.py exp_data.txt"); 
    return 0;
    
}