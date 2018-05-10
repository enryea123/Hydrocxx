#include <iostream>
#include <iomanip>
#include "math.h"
#include "shooting.h"
using namespace std;

shooting::shooting(int n, double stepsize, int extrasteps, int eig_order_max, double x0, double gravity){

	n_shoot_=10000;				// The bigger it is, the faster and more accurate is the program.
	n_=n;
	stepsize_=stepsize;
	extrasteps_=extrasteps;
	eig_order_max_=eig_order_max;
	x0_=x0;
	gravity_=gravity;

	omega_squared_guess_=10000.*(-gravity_);		// Initial guess for omega. It is safer to make it
													// depend on gravity, in case of g>100*gravity_sun
	vector_W = new double[n_shoot_];
	vector_F = new double[n_shoot_];
	eigenvalues = new double [n_shoot_];

	// The number of eigenvalues that it can find depends on gravity: for gravity_earth or sun 
	// it can find up to 200-300 of them, for gravity_moon/10 it starts to fail at around 17-18.
}

shooting::~shooting(){
	delete[] vector_W, vector_F, eigenvalues;
}

void shooting::change_parameters(int n, double stepsize, int extrasteps, int eig_order_max, double x0, double gravity){
	n_=n;
	stepsize_=stepsize;
	extrasteps_=extrasteps;
	eig_order_max_=eig_order_max;
	x0_=x0;
	gravity_=gravity;

	delete[] vector_W, vector_F, eigenvalues;
	vector_W = new double[n_shoot_];
	vector_F = new double[n_shoot_];
	eigenvalues = new double [n_shoot_];
}

double shooting::get_eigenvalues(int i){
	return eigenvalues[i];
}

void shooting::print_eigenvalues(int i){
	cout<<fixed<<setprecision(int(1+log10(1/stepsize_)));
	if (i+1<10)												// Just add an extra space for better formatting.
		cout<<"Omega_"<<i+1<<"  = "<<eigenvalues[i]<<endl;
	else
		cout<<"Omega_"<<i+1<<" = "<<eigenvalues[i]<<endl;
}

void shooting::shoot_now(equation* EQ, bool print_values){

	int count=1;											// Reset all the parameters when called.
	double ratio=1., omega_squared=omega_squared_guess_;
	RK = new rungekutta(n_, stepsize_, x0_);				// Initialize rungekutta.

	for(int eig_order_=0;eig_order_<eig_order_max_;eig_order_++){
		for (int i=0;i<n_shoot_;i++){

			if (omega_squared<ratio/50.){
				i=0;
				count++;
				ratio /= ((eig_order_+1.)*10.);				// Decreasing the ratio.
			}

			RK->solver(omega_squared,EQ);
			vector_W[i]=omega_squared;
			vector_F[i]=RK->getsolution(n_-extrasteps_);

			if (vector_F[i]==0 || ( i>0 && vector_F[i-1]*vector_F[i]<0 ) ){	// Rough eigenvalues finder

				double temp_eigenvalue=vector_W[i];

				for (int k=0;k<int(1+log10(1/stepsize_));k++){	// Accurate eigenvalues finder	k_max-1=number of digits
					for (int l=0;l<100;l++){

						double temp_omega_squared = temp_eigenvalue*(1.+ratio*pow(10,-k-1)*(l/50.-1) );	// Refine loop

						RK->solver(temp_omega_squared,EQ);
						vector_W[l]=temp_omega_squared;
						vector_F[l]=RK->getsolution(n_-extrasteps_);

						if (vector_F[l]==0 || (l>0 && vector_F[l-1]*vector_F[l]<0)){	// If a zero is found, the corresponding
							temp_eigenvalue=vector_W[l];								// omega^2 is set as an eigenvalue.
							break;
						}
					}
					if (fabs(RK->getsolution(n_-extrasteps_))<stepsize_)
						break;									// break if abs(F(w^2))<1e-3
				}	// end of accurate eigenvalues finder

				eigenvalues[eig_order_]=temp_eigenvalue;		// Fill eigenvalues vector
				if(print_values)
					print_eigenvalues(eig_order_);				// Print eigenvalue of order eig_order_
				break;
			}
			omega_squared *= 1.-1./((eig_order_/5.+1.)*pow(10,count));
		}
	}

	delete RK;
}	// end of shoot_now method

void shooting::bruteforce(){
	// not implemented, might be used to find eigenvalues of very high order, but it would be very slow
}
