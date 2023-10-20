#include <math.h>

double l1distance(double v1[], double v2[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(v1[i] - v2[i]);
	}
	return sum;
}