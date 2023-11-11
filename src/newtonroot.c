#include <stdio.h>
#include <math.h>

/*
	Approximate the root of a function using Newton's method
	
	Input:
		f: function to approximate the root of
		df: derivative of f
		x: starting point for x
		tolerance: determines precision of calculated root
		
	Output:
		Approximated root or NaN if a root within the given tolerance was not found
*/
double newtonroot(double (*f)(double), double (*df)(double), double x, double tolerance)
{
    const int maxiter = 15;
    double x1 = x;
    for (int i = 0; i < maxiter; i++)
    {
        double x2 = x1 - f(x1) / df(x1);
        if (fabs(x2 - x1) <= tolerance)
        {
            return x2;
        }
        x1 = x2;
    }
    // If root is not found
    printf("Root not found");
    return NAN;
}    