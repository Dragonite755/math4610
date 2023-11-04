#include <math.h>

/*
	Computes the L1 norm of a vector
	
	Input:
		n: length of v	
		v: vector to calculate norm of
		
	Output:
		return: norm of v
*/
double l1norm(int n, double v[])
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(v[i]);
	}
	return sum;
}