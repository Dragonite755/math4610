float dmaceps()
{
	double eps = 1.0;
	double prevEps;
	while (1.0 + eps != 1.0)
	{
		prevEps = eps;
		eps /= 2.0;
	}
	return prevEps;
}