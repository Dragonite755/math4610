#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backsub.c"
#include "gausselim.c"

/*
	Tests the routines required for solving a system of equations via Gaussian elimination
*/
int main()
{
	// Generate strictly diagonal dominant matrix a
	// System of equations: ay = b
	const int n = 10;
	double y[n] // All 1s
	double a[n][n]
	for (int i = 0; i < n; i++)
	{
		y[i] = 1.0;
		double sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			a[i][j] = (double) rand() / RAND_MAX;
			sum += a[i][j];
		}
		a[i][i] = sum + (double) rand() / RAND_MAX; // Establish strict diagonal dominance on the row
	}
	
	double b[n];
	matvec(n, a, y, b);
	
	// Solve ax = b
	double x[n];
	gausselim(n, a, b);
	backsub(n, a, b, x);
	
	// Solution should be x = y
	printf("Solution should be all 1s:\n[ ");
	for (int i = 0; i < n; i++)
	{
		printf("%.2f ", x[i]);
	}
	printf("]\n");

	return 0;
}