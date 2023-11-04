/*
	Computes the coefficients for a linear regression of values y on x (y = ax + b)
	
	Input:
		int n: length of input vectors x and y
		x: vector of values
		c: two-valued empty vector to store coefficients of regression
	
	Output:
		c: Stores the values of the coefficients of the regression [a, b]
*/
void linreg(int n, double x[], double y[], double c[]) {
	// Calculate sums first
	double s1 = 0;   // sum x_i
	double s2 = 0;   // sum x_i^2
	
	for (int i = 0; i < n; i++)
	{
		s1 += x[i];
		s2 += x[i] * x[i];
	}
	
	// Set values of coefficients vector to 0
	c[0] = 0;   // a
	c[1] = 0;   // b
	
	// Compute (X^T X)^-1 X^T y
	double det = n * s2 - s1 * s1;   // Divide by this last
	for (int i = 0; i < n; i++)
	{
		c[0] += y[i] * (n * x[i] - s1);
		c[1] += y[i] * (s2 - s1 * x[i]);
	}
	c[0] /= det;
	c[1] /= det;
}