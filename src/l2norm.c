#include <math.h>

/*
	Computes the L2 norm of a vector
	
	Input:
		n: length of v
		v: vector to calculate norm of
		
	Output:
		return: norm of v
*/
double l2norm(int n, double v[])
{
	double sum = 0.0;   // norm^2
	for (int i = 0; i < n; i++)
	{
		sum += v[i] * v[i];
	}
	return sqrt(sum);
}