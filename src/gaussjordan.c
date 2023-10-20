/*
Utilizes Gauss-Jordan elimination on ax = b to put a into upper triangular form, modifying b accordingly

n: size
a: matrix
b: vector
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