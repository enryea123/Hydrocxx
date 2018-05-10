#ifndef __density_h__
#define __density_h__
#include "vector.h"

class densityfunction{
	public:
		densityfunction(double);
		~densityfunction();
		double eval(double);					// Evaluate
		double diff(double);					// Differentiate
	private:
		double sigma_;
		static const double epsilon_=5e-8;		// Small epsilon for the numerical derivative.
};

#endif
