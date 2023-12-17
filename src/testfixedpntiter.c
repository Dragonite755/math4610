#include "fixedpntiter.c"

#include <stdio.h>
#include <math.h>

double f(double x)
{
    return x * x - 5 * x + 6;
}

double g1(double x)
{
    return x + f(x);
}

double g2(double x)
{
    return x - f(x);
}

int main()
{
    double tolerance = 1e-6;
    printf("f(x) = x^2 - 5x + 6\n\n");
    printf("g(x) = x + f(x)\n");
    for(int i = -15; i <= 15; i++)
    {
        double x0 = i / 3.0;
        printf("x0 = %4.1f: ", x0);
        double root = fixedpntiter(g1, x0, tolerance);
        if (!isnan(root))
        {
            printf("root at x = %f\n", root);
        }
        else
        {
            printf("\n");
        }
    }
    
    printf("\ng(x) = x - f(x)\n");
    for(int i = -15; i <= 15; i++)
    {
        double x0 = i / 3.0;
        printf("x0 = %4.1f: ", x0);
        double root = fixedpntiter(g2, x0, tolerance);
        if (!isnan(root))
        {
            printf("root at x = %f\n", root);
        }
        else
        {
            printf("\n");
        }
    }

    return 0;
}