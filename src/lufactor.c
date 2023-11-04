/*
	LU factorizes the matrix a, where the upper-right and diagonal elements are the upper triangular matrix
	and the lower-left elements are those of the lower triangular matrix

	Input:
		n: size of matrix
		a: n x n matrix
*/
void lufactor(int n, double a[][n])
{
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = a[i][k] / a[k][k];
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= factor * a[k][j];
            }
            a[i][k] = factor;
        }
    }
}