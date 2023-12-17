#include "bisectroot.c"

#include <stdio.h>

double f(double x)
{
    return x * x - 5.0 * x + 6.0;
}

int main()
{
    printf("f(x) = x^2 - 5x + 6\n");
    double a = 2.5;
    double b = 10.0;
    double tolerance = 1e-6;
    double root = bisectroot(f, a, b, tolerance);
    printf("[%.1f, %.1f]: %f\n", a, b, root);
    
    a = -10.0;
    b = 2.5;
    root = bisectroot(f, a, b, tolerance);
    printf("[%.1f, %.1f]: %f\n", a, b, root);

    return 0;
}