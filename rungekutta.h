#ifndef __rungekutta_h__
#define __rungekutta_h__
#include "vector.h"
#include "equation.h"

class solvingmethod{									// Base solving method class with the pure
	public:												// virtual method calculate_step, which must
		solvingmethod();								// be implemented in the derived classes,
		~solvingmethod();								// plus some useful methods and variables.

		void change_parameters(int,double,double);
		double getposition(int);
		double getsolution(int);
		double* getposition();
		double* getsolution();

		double* get_logerror_ex1();						// Two examples of error compared with
		double* get_logerror_ex2();						// an analitical solution.

		virtual linearvector calculate_step(double,linearvector&,double,equation*,double) const=0;
		void solver(double,equation*);
	protected:											// Protected and not private so the
		int n_;											// derived classes can modify them.
		double stepsize_, x0_, *vector_x, *vector_u, *logerror_ex1, *logerror_ex2;
};

class rungekutta: public solvingmethod{					// RK4 method implementation.
	public:
		rungekutta(int,double,double);
		~rungekutta();
		virtual linearvector calculate_step(double,linearvector&,double,equation*,double) const;
};

class euler: public solvingmethod{						// Euler method implementation.
	public:
		euler(int,double,double);
		~euler();
		virtual linearvector calculate_step(double,linearvector&,double,equation*,double) const;
};

#endif
