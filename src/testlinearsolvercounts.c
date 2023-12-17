#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backsubcount.c"
#include "forwardsubcount.c"
#include "gausselimcount.c"
#include "lufactorcount.c"
#include "matvec.c"

/*
	Tests the complexity (measured in number of mathematical operations) between Gaussian elimination and LU factorization
	in solving a system of equations
*/
int main()
{
	// Increment the value of n to find the breaking-even point
	srand(time(NULL));
	printf("Counting number of operations to solve n x n system of equations:\n\n");
	printf(" n || elim  | back  | Gauss || fact  | forw  | back  |    LU\n------------------------------------------------------------\n");
	for (int n = 2; n <= 50; n++)
	{
		printf("%2d || ", n);
		double a[n][n]; // Create separate matrices and vectors for each strategy
		double z[n]; // Name z to reserve x and y for LU factorization
		
		// Create strictly diagonal matrix a
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				a[i][j] = (double) rand() / RAND_MAX;
			}
			z[i] = 1.0;
		}
		
		// System of equations az = b
		double b[n];
		matvec(n, a, z, b);
		
		// Gaussian elimination
		double a_gauss[n][n]; // Make copies of a and b to preserve them for the lu method
		double b_gauss[n];
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				a_gauss[i][j] = a[i][j];
			}
			b_gauss[i] = b[i];
		}
		
		double x_gauss[n];
		int op_count = 0; // Count for a single routine in the method
		int op_total = 0; // Total operations for current method
		
		gausselimcount(n, a_gauss, b_gauss, &op_count);
		printf("%5d | ", op_count);
		op_total += op_count;
		op_count = 0;
		
		backsubcount(n, a_gauss, b_gauss, x_gauss, &op_count);
		op_total += op_count;
		printf("%5d | ", op_count);
		printf("%5d || ", op_total);
		
		// LU factorization
		op_count = 0;
		op_total = 0;
		
		lufactorcount(n, a, &op_count);
		printf("%5d | ", op_count);
		op_total += op_count;
		op_count = 0;
		
		double y[n];
		forwardsubcount(n, a, b, y, &op_count);
		printf("%5d | ", op_count);
		op_total += op_count;
		op_count = 0;
		
		double x_lu[n];
		backsubcount(n, a, y, x_lu, &op_count);
		printf("%5d | ", op_count);
		op_total += op_count;
		printf("%5d\n", op_total);
	}
	
	return 0;
}