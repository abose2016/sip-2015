#include <stdio.h>
#include <vector>
#include <gsl/gsl_interp.h>

addtwonumbers(int x = 5, int y = 7)
{
	int z = 0;
	printf("\nEnter first number: ");
	scanf("%d", &x);
	printf("\nEnter second number: ");
	scanf("%d", &y);
	z = x + y;
	printf("\nResult is %d\n", z);
} 
