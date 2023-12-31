/*
	Computes the forward difference quotient to approximate a derivative
	
	Input:
		f: function to approximate the derivative of
		x: value of x at which to approximate f'(x)
		h: increment for the approximation
		
	Output:
		return: approximation of f'(x)
*/
double forwarddf(double (*f)(double), double x, double h)
{
	return (f(x + h) - f(x)) / h;
}