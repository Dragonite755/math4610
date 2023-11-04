/*
	Computes the machine epsilon in double-precision
	
	Input:
		
	Output:
		return: Double-precision machine epsilon value
*/
double dmaceps()
{
	double eps = 1.0;
	double prevEps;   // Stores the last eps such that 1.0 + eps != 1.0
	while (1.0 + eps != 1.0)
	{
		prevEps = eps;
		eps /= 2.0;
	}
	return prevEps;
}