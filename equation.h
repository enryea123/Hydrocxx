#ifndef __equation_h__
#define __equation_h__
#include "vector.h"
#include "density.h"

class equation{
	public:
		equation(double,double,double);
		~equation();
		void change_parameters(double,double,double);
		double get_P_coeff(double,double);
		double get_Q_coeff(double,double);
		linearvector eval(double,double,linearvector);
	private:
		double kappa_squared_, gravity_;
		densityfunction* RHO;
};

#endif
