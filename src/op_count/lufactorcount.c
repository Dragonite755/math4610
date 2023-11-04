/*
	LU factorizes the matrix a, where the upper-right and diagonal elements are the upper triangular matrix
	and the lower-left elements are those of the lower triangular matrix

	Input:
		n: size of matrix
		a: n x n matrix
	
	Output:
		op_count: incrememented by the number of mathematical operations performed
*/
void lufactorcount(int n, double a[][n], int* op_count)
{
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = a[i][k] / a[k][k];
			*op_count += 1;
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= factor * a[k][j];
				*op_count += 2;
            }
            a[i][k] = factor;
        }
    }
}