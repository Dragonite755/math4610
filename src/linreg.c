void linreg(double x[], double y[], int n, double A[]) {
	double s1 = 0;   // sum x_i
	double s2 = 0;   // sum x_i^2
	
	for (int i = 0; i < n; i++)
	{
		s1 += x[i];
		s2 += x[i] * x[i];
	}
	
	// y = mx + b
	A[0] = 0;   // m
	A[1] = 1;   // b
	
	// (X^T X)^-1 X^T Y
	double det = n * s2 - s1 * s1;
	for (int i = 0; i < n; i++)
	{
		A[0] += y[i] * (n * x[i] - s1);
		A[1] += y[i] * (s2 - s1 * x[i]);
	}
	A[0] /= det;
	A[1] /= det;
}