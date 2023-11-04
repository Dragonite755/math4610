/*
	Solve au = b for x, where u is upper triangular

	Input:
		n: size of matrix and vectors
		u: upper triangular matrix
		b: vector
		x: vector (empty)
	
	Output:
		x: vector solution to ax = b
*/
void backsub(int n, double u[][n], double b[], double x[])
{
    x[n - 1] = b[n - 1] / u[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += u[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / u[i][i];
    }
}