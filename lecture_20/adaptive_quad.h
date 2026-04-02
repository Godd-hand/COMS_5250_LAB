#ifndef ADAPTIVE_QUAD_H
#define ADAPTIVE_QUAD_H

// integration functions
double Q_simpson(double a, double b, double (*f)(double, double), double param);
double AdaptiveIntSerial(double a, double b, double TOL, double (*f)(double, double), double param);

// Implementation of OpenMP Version 1
double AdaptiveInt_OMP_v1(double a, double b, double TOL, double (*f)(double, double), double param, int thread_count);

#endif