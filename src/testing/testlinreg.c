#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "linreg.c"

/*
	Tests linear regression by creating a line and x points, slightly offsetting
	y values, and comparing it with the initial line
*/
int main()
{
	// Create random line y = ax + b to approximate with linear regression
	srand(time(NULL));
	int n = 100;   // Enough points to smoothen out irregularities of offsets
	double a = -10 + (double) rand() / RAND_MAX * 20; // a in [-10, 10]
	double b = -10 + (double) rand() / RAND_MAX * 20; // b in [-10, 10]
	
	// Generate x values and y values slightly offset from the line
	double x[n];
	double y[n];
	for (int i = 0; i < n; i++)
	{
		x[i] = -10 + (double) rand() / RAND_MAX * 20; // x values in [-10, 10]
		y[i] = a * x[i] + b - 1 + (double) rand() / RAND_MAX * 2; // Offset by +- 1
	}
	
	// Test that linear regression line approximates initial line
	printf("Exact line:\n\ty = %fx + %f\n", a, b);
	
	double c[2];
	linreg(n, x, y, c);
	printf("Linear Regression:\n\ty = %fx + %f\n", c[0], c[1]);
	
	return 0;
}