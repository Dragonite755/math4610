float smaceps()
{
	float eps = 1.0f;
	float prevEps;
	while (1.0f + eps != 1.0f)
	{
		prevEps = eps;
		eps /= 2.0f;
	}
	return prevEps;
}