#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backsub.c"
#include "matvec.c"

void main() {
	// Generate 10x10 upper-triangular system of equations
	// ax = b
	const int n = 10;
	double x[n];
	double a[n][n];
	for (int i = 0; i < n; i++)
	{
		x[i] = 1.0;
		for (int j = i; j < n; j++)
		{
			a[i][j] = (double) rand() / RAND_MAX;
		}
	}
	
	double b[n];
	matvec(n, a, x, b);
	
	// Solve ax = b
	printf("Solution should be all 1s:\n[ ");
	backsub(n, a, b, x); // Re-use x vector, which will be filled with new values
	for (int i = 0; i < n - 1; i++)
	{
		printf("%.2f ", x[i]);
	}
	printf("]");
}