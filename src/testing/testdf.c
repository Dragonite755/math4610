#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "forwarddf.c"
#include "backwarddf.c"
#include "centraldf.c"

/*
	Tests approximations of derivatives utilizing difference quotients
*/
int main()
{
	srand(time(NULL));
	
	printf("Forward Difference Quotient:\n");
	double x = -5 + (double) rand() / RAND_MAX * 10;
	printf("\tApproximating sin(x) at x = %f\n", x);
	printf("\td/dx sin(%f) = %.10f\n\n", x, cos(x)); // d/dx sin(x) = cos(x)
	double h = 1.0;
	for (int i = 0; i < 10; i++)
	{
		printf("\th = %.10f : %.10f\n", h, forwarddf(sin, x, h));
		h /= 10;
	}
	
	printf("Backward Difference Quotient:\n");
	x = -5 + (double) rand() / RAND_MAX * 10;
	printf("\tApproximating exp(x) at x = %f\n", x);
	printf("\td/dx exp(%f) = %.10f\n\n", x, exp(x)); // d/dx e^x = e^x
	h = 1.0;
	for (int i = 0; i < 10; i++)
	{
		printf("\th = %.10f : %.10f\n", h, backwarddf(exp, x, h));
		h /= 10;
	}
	
	printf("Central Difference Quotient:\n");
	x = (double) rand() / RAND_MAX * 10;
	printf("\tApproximating sqrt(x) at x = %f\n", x);
	printf("\td/dx sqrt(%f) = %.10f\n\n", x, 1 / (2 * sqrt(x))); // d/dx sqrt(x) = 1 / (2 sqrt(x))
	h = 1.0;
	for (int i = 0; i < 10; i++)
	{
		printf("\th = %.10f : %.10f\n", h, backwarddf(sqrt, x, h));
		h /= 10;
	}
	
	return 0;
}