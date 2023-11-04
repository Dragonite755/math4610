#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backsub.c"
#include "matvec.c"

/*
	Tests back-substitution on an upper-triangular system of equations
*/
void main() {
	// Generate 10x10 upper-triangular system of equations
	// ay = b
	// y = [all 1s]
	const int n = 10;
	double y[n];
	double a[n][n];
	for (int i = 0; i < n; i++)
	{
		y[i] = 1.0;
		for (int j = i; j < n; j++)
		{
			a[i][j] = (double) rand() / RAND_MAX;
		}
	}
	
	double b[n];
	matvec(n, a, y, b);
	
	// Solve ax = b
	// y = x should be the solution
	printf("Solution should be all 1s:\n[ ");
	double x[n];
	backsub(n, a, b, x);
	for (int i = 0; i < n; i++)
	{
		printf("%.2f ", x[i]);
	}
	printf("]\n");
}