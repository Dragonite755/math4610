#include <math.h>

double linfdistance(double v1[], double v2[], int n)
{
	double maxAbs = 0.0;
	for (int i = 0; i < n; i++)
	{
		double abs = fabs(v1[i] - v2[i]);
		if (abs > maxAbs)
		{
			maxAbs = abs;
		}
	}
	return maxAbs;
}