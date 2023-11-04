/*
Solve ly = b for y, where l is a lower triangular matrix

n: size
l: lower triangular matrix
b: vector
y: vector
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