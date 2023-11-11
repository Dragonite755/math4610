#include <stdio.h>
#include <math.h>

double secantroot(double (*f)(double), double a, double b, double tolerance)
{
    const int maxiter = 15;
    double x0 = a;
    double x1 = b;
    double f1 = f(x1);
    for (int i = 0; i < maxiter; i++)
    {
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f(x0));
        x0 = x1;
        x1 = x2;
        f1 = f(x1);
        if (fabs(x1 - x0) < tolerance)
        {
            return x1;
        }
    }
    printf("Root not found");
    return NAN;
}