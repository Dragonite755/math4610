/*
	Utilizes Gauss-Jordan elimination to reduce ax = b into upper triangular form (modifying both a and b)

	Input:
		n: size of a and b
		a: square matrix in ax = b
		b: vector in ax = b
		
	Output:
		a: reduced to upper-triangular matrix
		b: changed to maintain ax = b
		op_count: incrememented by the number of mathematical operations taken
*/
void gausselimcount(int n, double a[][n], double b[], int* op_count)
{
    // Make a upper triangular (modifying b accordingly)
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
            b[i] -= factor * b[k];
			*op_count += 2;
        }
    }
}