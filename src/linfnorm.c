#include <math.h>

/*
	Computes the L-infinity norm of a vector
	
	Input:
		n: length of v
		v: vector to calculate norm of
		
	Output:
		return: norm of v
*/
double linfnorm(int n, double v[])
{
	double max = 0.0;
	for (int i = 0; i < n; i++)
	{
		double abs = fabs(v[i]);
		if (abs > max)
		{
			max = abs;
		}
	}
	return max;
}