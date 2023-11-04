/*
	Computes the backward difference quotient to approximate a derivative
	
	Input:
		f: function to approximate the derivative of
		x: value of x at which to approximate f'(x)
		h: increment for the approximation
		
	Output:
		return: approximation of f'(x)
*/
double backdf(double (*f)(double), double x, double h)
{
	return (f(x) - f(x - h)) / h;
}	