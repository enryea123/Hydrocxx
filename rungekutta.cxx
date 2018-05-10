#include "rungekutta.h"
#include "math.h"

solvingmethod::solvingmethod(){}

solvingmethod::~solvingmethod(){
	delete[] vector_x, vector_u, logerror_ex1, logerror_ex2;
}

void solvingmethod::change_parameters(int n, double stepsize, double x0){
	n_=n;
	stepsize_=stepsize;
	x0_=x0;
	delete[] vector_x, vector_u, logerror_ex1, logerror_ex2;
	vector_x = new double[n_];
	vector_u = new double[n_];
	logerror_ex1 = new double[n_];
	logerror_ex2 = new double[n_];
}

double solvingmethod::getposition(int i){
	return vector_x[i];
}
double solvingmethod::getsolution(int i){
	return vector_u[i];
}

double* solvingmethod::getposition(){
	return vector_x;
}
double* solvingmethod::getsolution(){
	return vector_u;
}

double* solvingmethod::get_logerror_ex1(){
	return logerror_ex1;
}

double* solvingmethod::get_logerror_ex2(){
	return logerror_ex2;
}

void solvingmethod::solver(double omega_squared, equation* EQ){
	double x0=x0_;																		// Reset x0 to the starting value.
	linearvector initial_data(2);
	initial_data.setx(0,0.);															// Setting the initial data xi(0)=0.
	initial_data.setx(1,EQ->get_P_coeff(initial_data.getx(0),omega_squared));			// Here I am setting phi(0)=P(0,omega^2)=omega^2,
																						// so that xi'(0) = phi(0)/P(0,omega^2) = 1.
	for (int i=0;i<n_;i++){
		vector_x[i]=x0;																	// Setting the position and solution (u) vector,
		vector_u[i]=initial_data.getx(0);												// to get them use the get methods.

		logerror_ex1[i]=log10( fabs(initial_data.getx(0)-sinh(x0)) );					// This is the error only for k^2=1,sigma=0
		logerror_ex2[i]=log10( fabs(initial_data.getx(0)-log(2.*x0+1.)/2.) );			// This is the error only for k^2=0,sigma=2

		initial_data = calculate_step(x0,initial_data,stepsize_,EQ,omega_squared);		// Using iteration on initial_data.
		x0=x0+stepsize_;																// Moving the starting point ahead.
	}
}

//------------------------------------------------RK4------------------------------------------------//

rungekutta::rungekutta(int n, double stepsize, double x0){
	n_=n;
	stepsize_=stepsize;
	x0_=x0;
	vector_x = new double[n_];
	vector_u = new double[n_];
	logerror_ex1 = new double[n_];
	logerror_ex2 = new double[n_];
}

rungekutta::~rungekutta(){}

linearvector rungekutta::calculate_step(double x,linearvector& initial_data, double h, equation* EQ, double omega_squared) const{
	linearvector v(2);
	linearvector k1=EQ->eval(x,omega_squared,initial_data);						// calculating the k_i
	linearvector k2=EQ->eval(x+h/2,omega_squared,k1*(h/2)+initial_data);
	linearvector k3=EQ->eval(x+h/2,omega_squared,k2*(h/2)+initial_data);		// h=stepsize
	linearvector k4=EQ->eval(x+h,omega_squared,k3*h+initial_data);
	v=((k2+k3)*2+k1+k4)*(h/6) + initial_data;									// Runge-Kutta
	return v;
}

//-----------------------------------------------Euler-----------------------------------------------//

euler::euler(int n, double stepsize, double x0){
	n_=n;
	stepsize_=stepsize;
	x0_=x0;
	vector_x = new double[n_];
	vector_u = new double[n_];
	logerror_ex1 = new double[n_];
	logerror_ex2 = new double[n_];
}

euler::~euler(){}

linearvector euler::calculate_step(double x,linearvector& initial_data, double h, equation* EQ, double omega_squared) const{
	linearvector v(2);
	v=(EQ->eval(x,omega_squared,initial_data))*h + initial_data;				// Euler method
	return v;
}
