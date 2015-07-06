#include <iostream>
#include <gsl/gsl_math.h>

void gslDemo1()
{
	double x = 5.0;
	int y = 2;
	double z = gsl_pow_int(x, y);

	std::cout<<z << endl;
}
