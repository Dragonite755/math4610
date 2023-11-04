#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "l1norm.c"
#include "l2norm.c"
#include "linfnorm.c"
#include "l1dist.c"
#include "l2dist.c"
#include "linfdist.c"

/*
	Tests all norm length and distance routines
*/
int main() 
{
	// Generate random test vectors
	const int nmin = 1;        // min value of size of test vectors
	const int nmax = 5;        // max value of size of test vectors
	const double vmin = -10.0; // min value of elements in test vectors
	const double vmax = 10.0;  // max value of elements in test vectors
	
	srand(time(NULL));
	int n = nmin + rand() % nmax; // size of test vectors
	
	double v1[n]; // test vector
	double v2[n]; // test vector
	for (int i = 0; i < n; i++)
	{
		v1[i] = vmin + (double) rand() / RAND_MAX * (vmax - vmin);
		v2[i] = vmin + (double) rand() / RAND_MAX * (vmax - vmin);
	}
	
	// Display test vectors
	printf("Randomly-Generated Test Vectors:\n-----\nv1:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%f\n", v1[i]);
	}
	printf("\nv2:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%f\n", v2[i]);
	}
	
	// Test norms and distances
	printf("\nTesting norms and distances:\n-----\n");
	
	printf("L1\n");
	printf("\tNorm v1         : %f\n", l1norm(n, v1));
	printf("\tNorm v2         : %f\n", l1norm(n, v2));
	printf("\tDistance v1, v2 : %f\n", l1dist(n, v1, v2));
	
	printf("L2\n");
	printf("\tNorm v1         : %f\n", l2norm(n, v1));
	printf("\tNorm v2         : %f\n", l2norm(n, v2));
	printf("\tDistance v1, v2 : %f\n", l2dist(n, v1, v2));
	
	printf("L-Infinity\n");
	printf("\tNorm v1         : %f\n", linfnorm(n, v1));
	printf("\tNorm v2         : %f\n", linfnorm(n, v2));
	printf("\tDistance v1, v2 : %f\n", linfdist(n, v1, v2));
	
	return 0;
}