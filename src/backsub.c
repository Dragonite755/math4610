/*
	Solve ax = b for x, where a is upper triangular

	Input:
		n: size of matrix and vectors
		a: upper triangular matrix
		b: vector
		x: vector (empty)
	
	Output:
		x: vector solution to ax = b
*/
void backsub(int n, double a[][n], double b[], double x[])
{
    x[n - 1] = b[n - 1] / a[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}