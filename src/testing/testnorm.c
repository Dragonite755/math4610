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
	Tests all norm and distance routines
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
	
	printf("Test vectors:\n");
	printf("1: \n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%f", v1[i]);
	}
	printf("\n2: \n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%f", v2[i]);
	}
	
	return 0;
}