#include <stdio.h>
#include <math.h>

/*
	Utilizes the bisection method to solve for a root of a function n a given interval
	
	Input:
		f: function to find root of
		a: left-point of interval to find root in
		b: right-point of interval to find root in
		tolerance: determines the precision of the root
		
	Output:
		Approximate root or NaN if no root was found within the given tolerance
*/
double bisectroot(double (*f)(double), double a, double b, double tolerance)
{
    double x1 = a;
    double x2 = b;
    
    // Check that both values are at opposite signs
    double f1 = f(x1);
    double f2 = f(x2);
    if (f1 * f2 >= 0)
    {
        return NAN; // No root guaranteed within the interval
    }
    
    // Iterate
    double xmid = (x1 + x2) / 2.0; // Midpoint
    double fmid = f(xmid);
    while (x2 - x1 > tolerance)
    {
        if (fmid == 0)
        {
            return fmid;
        }
        else if (f1 * fmid < 0)
        {
            x2 = xmid;
            f2 = fmid;
        }
        else
        {
            x1 = xmid;
            f1 = fmid;
        }
        xmid = (x1 + x2) / 2.0; // Midpoint
        fmid = f(xmid);
    }
    return xmid;
}