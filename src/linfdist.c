#include <math.h>

/*
	Computes the L-infinity distance between two vectors
	
	Input:
		n: length of v1 and v2
		v1: input vector to calculate distance from
		v2: input vector to calculate distance from
		
	Output:
		return: distance between v1 and v2
*/
double linfdistance(double v1[], double v2[], int n)
{
	double max = 0.0;
	for (int i = 0; i < n; i++)
	{
		double abs = fabs(v1[i] - v2[i]);
		if (abs > max)
		{
			max = abs;
		}
	}
	return max;
}