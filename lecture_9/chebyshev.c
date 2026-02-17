/*
Expansion with Chebyshev polynomials (0 <= N <= 5).
Use a switch statement.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void SamplePoly(const int N, const int NumPts, const double b[], const double x[], double y[]);
void WritePoly(const int NumPts, const double x[], const double y[]);

int main()
{
    const int Nmax = 5;
    int N;

    // Input Polynomial Degree
    printf("\n Input polynomial degree (0-%i): ", Nmax);
    scanf("%i", &N);

    // Validate Input
    if (N < 0 || N > Nmax)
    {
        printf(" Invalid value N=%i.\n", N);
        printf(" N must satisfy: 0 <= N <= %i\n\n", Nmax);
        exit(1);
    }
    printf("\n");

    // read-in Coefficients
    double b[Nmax + 1];
    // Initialize all to 0.0 just in case
    for(int k=0; k<=Nmax; k++) b[k] = 0.0;

    for (int i = 0; i <= N; i++)
    {
        printf(" Set b[%i]: ", i);
        scanf("%lf", &b[i]);
    }
    printf("\n");

    // Set x-coordinates 
    const int NumPts = 21;
    double x[NumPts];
    double y[NumPts];

    for (int i = 0; i < NumPts; i++)
    {
        // Generates points from -1.0 to 1.0 
        x[i] = -1.0 + i * (2.0 / (1.0 * (NumPts - 1)));
    }

    // Calculate Polynomial Values
    SamplePoly(N, NumPts, b, x, y);

    //  Write Data to File
    WritePoly(NumPts, x, y);

    // Call Python to Plot
    printf(" Calling Python script to plot...\n");
    system("/home/george/miniconda3/envs/agentic/bin/python PlotPoly.py");

    return 0;
}

/* 
Function to evaluate the polynomial at all x points, using a switch statement
*/
void SamplePoly(const int N, const int NumPts, const double b[], const double x[], double y[])
{
    for (int i = 0; i < NumPts; i++)
    {
        const double a = x[i];
        double phi;
        
        // Initialize with the 0th term: b[0] * phi_0(x) where phi_0(x) = 1
        y[i] = b[0];

        switch(N)
        {
            case 5:
                // phi_5 = 16x^5 - 20x^3 + 5x 
                phi = 16.0*pow(a, 5) - 20.0*pow(a, 3) + 5.0*a;
                y[i] += b[5] * phi;
            case 4:
                // phi_4 = 8x^4 - 8x^2 + 1 
                phi = 8.0*pow(a, 4) - 8.0*pow(a, 2) + 1.0;
                y[i] += b[4] * phi;
            case 3:
                // phi_3 = 4x^3 - 3x 
                phi = 4.0*pow(a, 3) - 3.0*a;
                y[i] += b[3] * phi;
            case 2:
                // phi_2 = 2x^2 - 1 
                phi = 2.0*pow(a, 2) - 1.0;
                y[i] += b[2] * phi;
            case 1:
                // phi_1 = x 
                phi = a;
                y[i] += b[1] * phi;
                break; 
            case 0:
                break;
            default:
                printf("\n Error:\n");
                exit(1);
        }
    }
}

/*
write x and y arrays to 'poly.data'
*/
void WritePoly(const int NumPts, const double x[], const double y[])
{
    FILE *fid = fopen("poly.data", "w");
    if (fid == NULL)
    {
        printf("Error opening file\n");
        exit(1);
    }

    for (int i = 0; i < NumPts; i++)
    {
        fprintf(fid, "%lf %lf\n", x[i], y[i]);
    }
    fclose(fid);
    printf(" Data written to poly.data\n");
}