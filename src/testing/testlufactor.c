#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backsub.c"
#include "forwardsub.c"
#include "lufactor.c"
#include "matvec.c"

/*
	Tests the routines required for solving a system of equations via Gaussian elimination
*/
int main()
{
	// Generate strictly diagonal dominant matrix a
	// System of equations: az = b
	const int n = 10;
	double z[n]; // All 1s
	double a[n][n];
	for (int i = 0; i < n; i++)
	{
		z[i] = 1.0;
		double sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			a[i][j] = (double) rand() / RAND_MAX;
			sum += a[i][j];
		}
		a[i][i] = sum + (double) rand() / RAND_MAX; // Establish strict diagonal dominance on the row
	}
	
	double b[n];
	matvec(n, a, z, b);
	
	// Lu factor a into a = lu
	lufactor(n, a);
	
	// Solve ly = b
	double y[n];
	forwardsub(n, a, b, y);
	
	// Solve ux = y
	double x[n];
	backsub(n, a, y, x);
	
	// Solution should be x = z
	printf("Solution should be all 1s:\n[ ");
	for (int i = 0; i < n; i++)
	{
		printf("%.2f ", x[i]);
	}
	printf("]\n");

	return 0;
}