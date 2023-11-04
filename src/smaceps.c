/*
	Computes the machine epsilon in single-precision
	
	Input:
		
	Output:
		return: Single-precision machine epsilon value
*/
float smaceps()
{
	float eps = 1.0f;
	float prevEps;   // Stores the last eps such that 1.0f + eps != 1.0f
	while (1.0f + eps != 1.0f)
	{
		prevEps = eps;
		eps /= 2.0f;
	}
	return prevEps;
}