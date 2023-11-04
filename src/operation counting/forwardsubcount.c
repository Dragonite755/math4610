/*
	Solve ly = b for y, where l is a lower triangular matrix with diagonal values 1. Additionally counts the number of operations needed.

	Input:
		n: size
		l: lower triangular matrix of size n x n
		b: vector of size n
		y: vector of size n
		op_count: will be incremented by the number of mathematical operations taken
*/
void forwardsub(int n, double l[][n], double b[], double y[], int* op_count)
{
    y[0] = b[0];
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += l[i][j] * y[j];
			*op_count += 1;
        }
        y[i] = b[i] - sum;
		*op_count += 1;
    }
}