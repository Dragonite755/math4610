double backwardDf(double (*f)(double), double x, double h)
{
	return (f(x) - f(x - h)) / h;
}	