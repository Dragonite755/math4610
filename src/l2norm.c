#include <math.h>

double l2norm(double v[], int n)
{
	double sum = 0.0;   // norm^2
	for (int i = 0; i < n; i++)
	{
		sum += v[i] * v[i]);
	}
	return sqrt(sum);
}