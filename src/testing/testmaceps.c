#include <stdio.h>
#include "smaceps.c"
#include "dmaceps.c"

int main()
{
	printf("Single-precision machine epsilon: %.50f\n", smaceps());
	printf("Double-precision machine epsilon: %.50f\n", dmaceps());
	
	return 0;
}