#include "equation.h"

equation::equation(double kappa_squared, double gravity, double sigma){
	kappa_squared_=kappa_squared;
	gravity_=gravity;
	RHO = new densityfunction(sigma);
}

equation::~equation(){
	delete RHO;
}

void equation::change_parameters(double kappa_squared, double gravity, double sigma){
	kappa_squared_=kappa_squared;
	gravity_=gravity;
	delete RHO;
	RHO = new densityfunction(sigma);
}

double equation::get_P_coeff(double x,double omega_squared){
	return omega_squared*RHO->eval(x);
}

double equation::get_Q_coeff(double x,double omega_squared){
	return kappa_squared_*(RHO->eval(x)*omega_squared+gravity_*RHO->diff(x));
}

linearvector equation::eval(double x,double omega_squared,linearvector vec){
	linearvector v(2);
	v.setx(0,vec.getx(1)/get_P_coeff(x,omega_squared));
	v.setx(1,vec.getx(0)*get_Q_coeff(x,omega_squared));
	return v;
}
