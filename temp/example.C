#include <iostream>
#include <vector>
#include <gsl/gsl_interp.h>

int main()
{
    int y = 4;
    std::vector< double > vx;
    vx.push_back(3);
    y = vx.size();

    std::cout <<y;
}
