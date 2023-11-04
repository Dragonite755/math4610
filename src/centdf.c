/*
	Computes the central difference quotient to approximate a derivative
	
	Input:
		f: function to approximate the derivative of
		x: value of x at which to approximate f'(x)
		h: increment for the approximation
		
	Output:
		return: approximation of f'(x)
*/
double centdf(double (*f)(double), double(x), double h)
{
	return ((f(x + h) - f(x - h)) / (2 * h));
}