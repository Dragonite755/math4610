/*
	Solve ly = b for y, where l is a lower triangular matrix with diagonal values 1

	Input:
		n: size
		l: lower triangular matrix of size n x n
		b: vector of size n
		y: vector of size n
*/
void forwardsub(int n, double l[][n], double b[], double y[])
{
    y[0] = b[0];
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += l[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
}