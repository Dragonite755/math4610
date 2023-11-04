/*
	Utilizes Gauss-Jordan elimination to reduce ax = b into upper triangular form (modifying both a and b)

	Input:
		n: size of a and b
		a: square matrix in ax = b
		b: vector in ax = b
		
	Output:
		a: reduced to upper-triangular matrix
		b: changed to maintain ax = b
*/
void gaussjordan(int n, double a[][n], double b[])
{
    // Make a upper triangular (modifying b accordingly)
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = a[i][k] / a[k][k];
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
}