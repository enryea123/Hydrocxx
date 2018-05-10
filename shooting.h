#ifndef __shooting_h__
#define __shooting_h__
#include "equation.h"
#include "rungekutta.h"

class shooting{
	public:
		shooting(int,double,int,int,double,double);
		~shooting();

		void change_parameters(int,double,int,int,double,double);
		void shoot_now(equation*,bool);					// The bool is to decide if to print.
		void bruteforce();								// not implemented

		double get_eigenvalues(int);
		void print_eigenvalues(int);

	private:
		int n_, n_shoot_, extrasteps_, eig_order_max_;
		double stepsize_, x0_, omega_squared_guess_, gravity_;
		double *vector_W, *vector_F, *eigenvalues;
		rungekutta *RK;
};

#endif
