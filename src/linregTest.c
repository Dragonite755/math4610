#include <stdio.h>
#include "linreg.c"


int main() {
	double x[] = {2, 3, 5, 7, 9};
	double y[] = {4, 5, 7, 10, 15};
	int n = sizeof(x) / sizeof(x[0]);
	double A[2];
	linreg(x, y, n, A);
	
	printf("y = %fx + %f", A[0], A[1]);
}