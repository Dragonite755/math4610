#include <stdio.h>
#include <math.h>

/*
	Approximates a fixed point of a function g
	
	Input:
		g: function to find fixed point of
		x: starting x value
		tolerance: tolerance value (used as a stopping condition)
		
	Output:
		approximation of fixed point or NaN if not found
*/
double fixedpntiter(double (*g)(double), double x, double tolerance)
{
    const int maxiter = 15;
    double x1 = x;
    for (int i = 0; i < maxiter; i++)
    {
        double x2 = g(x1);
        if (fabs(x2 - x1) <= tolerance)
        {
            return x2;
        }
        x1 = x2;
    }
    printf("Failed to find fixed point");
    return NAN;
}