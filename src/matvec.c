/*
	Multiply a vector by a matrix
	
	Input:
		n: size of square matrix a and vector n
		m: square (n x n) matrix
		v: vector
		
	Output:
		p: product m * v
*/
void matvec(int n, double m[][n], double v[], double p[])
{
	for (int i = 0; i < n; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			sum += m[i][j] * v[j];
		}
		p[i] = sum;
	}
}