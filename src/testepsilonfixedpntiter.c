#include "fixedpntiter.c"

#include <stdio.h>

double f(double x)
{
    return x * x - 5.0 * x + 6.0;
}

double g(double x)
{
    return x - f(x) / 10.0;
}

int main()
{
    printf("f(x) = x^2 - 5x + 6\n");
    printf("g(x) = x - 0.1 f(x)\n");
    const double tolerance = 1e-3;
    for (int i = 80; i <= 140; i++)
    {
        double x0 = i / 10.0;
        printf("x0 = %4.1f: ", x0);
        double root = fixedpntiter(g, x0, tolerance);
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