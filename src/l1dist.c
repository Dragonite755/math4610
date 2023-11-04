#include <math.h>

/*
	Computes the L1 distance between two vectors
	
	Input:
		n: length of v1 and v2
		v1: input vector to calculate distance from
		v2: input vector to calculate distance from
		
	Output:
		return: distance between v1 and v2
*/
double l1distance(double v1[], double v2[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(v1[i] - v2[i]);
	}
	return sum;
}