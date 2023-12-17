#include "newtonroot.c"

#include <stdio.h>
#include <math.h>

double f(double x)
{
    return x * x - 5 * x + 6;
}

double df(double x)
{
    return 2 * x - 5;
}

int main()
{
    double tolerance = 1e-6;
    printf("f(x) = x^2 - 5x + 6\n");
    for (int i = -15; i <= 15; i++)
    {
        double x0 = i / 3.0;
        printf("x0 = %4.1f: ", x0);
        double root = newtonroot(f, df, x0, tolerance);
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