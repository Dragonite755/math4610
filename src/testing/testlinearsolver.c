#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backsub.c"
#include "gaussjordan.c"
#include "l1norm.c"

/*
	Tests linear solver routines
*/
int main()
{
	// Randomly generate an LU factorization, then un-factor it to ensure a diagonally
	// dominant matrix
	srand(time(NULL));
	const int n = 10;
	int i = n;
	double l[n][n];
	double u[n][n];
	const double minvalue = -10;
	const double maxvalue = 10;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				l[i][j] = 0.0;
				u[i][j] = minvalue + (double) rand() / RAND_MAX * (maxvalue - minvalue);
			}
			else if (i == j)
			{
				l[i][j] = 1.0;
				u[i][j] = minvalue + (double) rand() / RAND_MAX * (maxvalue - minvalue);
			}
			else {
				l[i][j] = minvalue + (double) rand() / RAND_MAX * (maxvalue - minvalue);
				u[i][j] = 0.0;
			}
		}
	}
	
	double a[n][n];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = 0.0;
			for (int k = 0; k < n; k++)
			{
				a[i][j] += l[i][k] * u[k][j];
			}
		}
	}
	
	// Create test ay = b
	double y[n]; // Exact solution x (all 1s)
	double b[n];
	
	for (int i = 0; i < n; i++)
	{
		y[i] = 1.0;
	}
	for (int i = 0; i < n; i++)
	{
		b[i] = 0.0;
		for (int j = 0; j < n; j++)
		{
			b[i] += a[i][j] * y[j];
		}
	}
	
	// Test gaussian linear solver to solve ax = b
	gaussjordan(n, a, b);
	
	double x[n]; // Solved solution
	backsub(n, a, b, x);
	
	// Display data
	printf("Testing Gauss-Jordan elimination and back substitution\n");
	printf("a:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t");
		for (int j = 0; j < n; j++)
		{
			printf("%7.3f ", a[i][j]);
		}
		printf("\n");
	}
	
	printf("b:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%10.6f\n", b[i]);
	}
	
	printf("Exact x:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%15.12f\n", y[i]);
	}
	
	printf("Solved x:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%15.12f\n", x[i]);
	}
	
	// Compute error
	double error[n];
	for (int i = 0; i < n; i++)
	{
		error[i] = x[i] - y[i];
	}
	
	printf("Error:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%15.12f\n", error[i]);
	}
	
	return 0;
}