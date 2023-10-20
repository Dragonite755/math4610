#include <math.h>

double l1norm(double v[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(v[i]);
	}
	return sum;
}