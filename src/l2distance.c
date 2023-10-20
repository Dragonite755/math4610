#include <math.h>

double l2distance(double v1[], double v2[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		double x = v1[i] - v2[i];
		sum += x * x;
	}
	return sqrt(sum);
}