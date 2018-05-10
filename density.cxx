#include "density.h"
#include "math.h"

densityfunction::densityfunction(double sigma){
	sigma_=sigma;
}

densityfunction::~densityfunction(){}

double densityfunction::eval(double x){							// Density function formula.
	return 1+sigma_*x;
}

double densityfunction::diff(double x){
	return (eval(x+epsilon_)-eval(x-epsilon_))/(2*epsilon_);	// Numerical derivative of
}																// the density function.
