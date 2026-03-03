/*
Generate normal distribution noise
*/

#include <stdlib.h>
#include <math.h>
#include "noise.h"

#define PI 3.14159265358979323846

// Generates normal distribution noise using Box-Muller Transform
double generate_normal(double sigma) {
    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;

    if (u1 <= 0.0) u1 = 1e-10;

    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    return z0 * sigma;
}