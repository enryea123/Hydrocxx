// Written by Enrico Amedeo Albano, student number r0650552.

// Compiled and tested on Debian 8.6.0-amd64, with g++ 4.9.2-2, 
// using root-system 5.34.19+dfsg-1.2 (http://root.cern.ch).
// To compile, type "make" in the terminal, to execute type "./program".
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>										// Library to check the elapsed time.
#include "TApplication.h"								// All the T***.h headers are for Root Cern.
#include "TROOT.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TAttLine.h"
#include "equation.h"
#include "rungekutta.h"
#include "shooting.h"
#define g_0 -9.80665
#define gravity_earth	1	  * g_0
#define gravity_sun		28.02 * g_0
#define gravity_mars	0.376 * g_0
#define gravity_moon	0.165 * g_0
using namespace std;

int main(){

	cout<<"Calculating..."<<endl<<endl;
	TApplication theApp("App",0,0);						// Initialize Root Cern.
	clock_t time=clock();								// Start counting elapsed time.							

	double stepsize=1e-3, x0=0., x1=1.;					// Stepsize and boundaries.
	int extrasteps=50;									// Extra steps to avoid the computation
														// to stop exactly at x_1.
	int n=int((x1-x0)/stepsize+extrasteps+0.5);			// Lenght of the vectors (number of
														// points), +0.5 to avoid truncation.
	double gravity=gravity_earth, omega_squared=1.;		// Parameters (omega_guess is defined
	double sigma=2, kappa_squared=1.;					// inside shooting.cxx !!!).

	int eig_order_max=25;								// Number of eigenvalues wanted. Put
														// a value between 1 and 250.

	equation *EQU;										// Creating the objects needed.
	shooting *SHO;
	rungekutta *RK;
	euler *EU;

	EQU = new equation(kappa_squared,gravity,sigma);	// Omega^2 is not passed as a parameter.
	SHO = new shooting(n,stepsize,extrasteps,eig_order_max,x0,gravity);
	RK = new rungekutta(n,stepsize,x0);
	EU = new euler(n,stepsize,x0);

	SHO->shoot_now(EQU,true);							// Use the shooting method now!
														// true="print eigenvalues"

//---------------------------------------------------------------------------------------------------//
//---------------------------------------------PLOTTING----------------------------------------------//
//---------------------------------------------------------------------------------------------------//

//Most of the code in this section is related to the style of the plots, so it won't be commented much.
	char legend[50];
//--------------------------------Some eigenvalues plotted together----------------------------------//

	TCanvas *c1 = new TCanvas("c1","Some eigenvalues",800,600);
	c1->SetGrid();
//	c1->Divide(3,2);
//	TPad *p1 = (TPad *)(c1->cd(k+1));
//	p1->SetGrid();
	TGraph *g1[5];
	TLegend *l1 = new TLegend(0.95,0.7,0.8,0.95);
	l1->SetTextAlign(32);

	for(int z=0;z<4;z++){
		RK->solver(SHO->get_eigenvalues(z),EQU);		// Plots the first 4 omega^2 and the last one.
		g1[z] = new TGraph(n,RK->getposition(),RK->getsolution());
		sprintf(legend,"#omega_{%d}^{2} = %0.4f",z+1,SHO->get_eigenvalues(z));
		l1->AddEntry(g1[z],legend,"l");
	}
	RK->solver(SHO->get_eigenvalues(eig_order_max-1),EQU);
	g1[4] = new TGraph(n,RK->getposition(),RK->getsolution());
	sprintf(legend,"#omega_{%d}^{2} = %0.4f",eig_order_max,SHO->get_eigenvalues(eig_order_max-1));
	l1->AddEntry(g1[4],legend,"l");

	g1[0]->GetYaxis()->SetNdivisions(5);
	g1[0]->SetTitle(" ");
	g1[0]->GetXaxis()->SetRangeUser(0.,1.);
	g1[0]->GetYaxis()->SetRangeUser(-0.25,0.25);
	g1[0]->GetXaxis()->SetTitle("x");
	g1[0]->GetYaxis()->SetTitle("#xi(x)");

	g1[0]->SetLineColor(2);
	g1[0]->Draw("AL");
	for(int z=1;z<5;z++){
		g1[z]->SetLineColor(2+z);
		g1[z]->Draw("same");
	}

	l1->Draw();
	c1->Modified();
	c1->Update();
	c1->SaveAs("eigenvalues.eps");

	c1->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");

	gROOT->SetBatch(kTRUE);

//----------------------------------------- Omega vs Kappa ------------------------------------------//

	cout<<endl<<"Calculating Omega^2(k^2)..."<<endl;
	eig_order_max=5;

	TCanvas *c2 = new TCanvas("c2","Omega vs kappa",800,600);
	c2->SetGrid();
	TGraph *g2[eig_order_max];
	TLegend *l2 = new TLegend(0.95,0.68,0.88,0.95);
	l2->SetTextAlign(12);

	int loops=20;
	double v_k[loops], v_w1[loops], v_w2[loops], v_w3[loops], v_w4[loops], v_w5[loops];

	SHO->change_parameters(n,stepsize,extrasteps,eig_order_max,x0,gravity);

	for (int i=0;i<loops;i++){
		v_k[i]=5.*pow(i,1.5)+0.01;						// Fills the k^2 vector.

		EQU->change_parameters(v_k[i],gravity,sigma);
		SHO->shoot_now(EQU,false);

		v_w1[i]=SHO->get_eigenvalues(0);				// Fills 5 vector, each corresponding
		v_w2[i]=SHO->get_eigenvalues(1);				// to an eigenvalue.
		v_w3[i]=SHO->get_eigenvalues(2);
		v_w4[i]=SHO->get_eigenvalues(3);
		v_w5[i]=SHO->get_eigenvalues(4);
	}

	g2[0] = new TGraph(loops,v_k,v_w1);
	g2[1] = new TGraph(loops,v_k,v_w2);
	g2[2] = new TGraph(loops,v_k,v_w3);
	g2[3] = new TGraph(loops,v_k,v_w4);
	g2[4] = new TGraph(loops,v_k,v_w5);

	for (int i=0;i<eig_order_max;i++){
		sprintf(legend,"#omega_{%d}^{2}",i+1);
		l2->AddEntry(g1[i],legend,"l");
	}

	g2[0]->GetXaxis()->SetNdivisions(5);
	g2[0]->GetYaxis()->SetNdivisions(5);
	g2[0]->SetTitle("");
	g2[0]->GetXaxis()->SetRangeUser(0.,200.);
//	g2[0]->GetYaxis()->SetRangeUser(-2.,2.);
	g2[0]->GetXaxis()->SetTitle("k^{2}");
	g2[0]->GetYaxis()->SetTitle("#omega^{2}(k^{2})");

	g2[0]->SetLineColor(2);
	g2[0]->Draw("AL");
	for(int z=1;z<5;z++){
		g2[z]->SetLineColor(2+z);
		g2[z]->Draw("same");
	}

	l2->Draw();
	c2->Modified();
	c2->Update();
	c2->SaveAs("omega-vs-kappa.eps");

	c2->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");

//---------------------------------------- Numerical Errors -----------------------------------------//

	TCanvas *c3 = new TCanvas("c3","Numerical error",800,600);
	c3->SetGrid();
	TGraph *g3[6];
	TLegend *l3 = new TLegend(0.95,0.68,0.53,0.95);		// For more information about the
	l3->SetTextAlign(32);								// errors, read rungekutta.cxx.

	omega_squared=1; sigma=0; kappa_squared=1; gravity=gravity_earth;
	EQU->change_parameters(kappa_squared,gravity,sigma);

	RK->solver(omega_squared,EQU);
	g3[0] = new TGraph(n,RK->getposition(),RK->get_logerror_ex1());
	sprintf(legend,"RK4, h = %0.3f, k^{2} = %0.2f, #sigma = %0.2f",stepsize,kappa_squared,sigma);
	l3->AddEntry(g3[0],legend,"l");

	EU->solver(omega_squared,EQU);
	g3[1] = new TGraph(n,EU->getposition(),EU->get_logerror_ex1());
	sprintf(legend,"Euler, h = %0.3f, k^{2} = %0.2f, #sigma = %0.2f",stepsize,kappa_squared,sigma);
	l3->AddEntry(g3[1],legend,"l");

	omega_squared=1; sigma=2; kappa_squared=0; gravity=gravity_earth;
	EQU->change_parameters(kappa_squared,gravity,sigma);

	RK->solver(omega_squared,EQU);
	g3[2] = new TGraph(n,RK->getposition(),RK->get_logerror_ex2());
	sprintf(legend,"RK4, h = %0.3f, k^{2} = %0.2f, #sigma = %0.2f",stepsize,kappa_squared,sigma);
	l3->AddEntry(g3[2],legend,"l");

	EU->solver(omega_squared,EQU);
	g3[3] = new TGraph(n,EU->getposition(),EU->get_logerror_ex2());
	sprintf(legend,"Euler, h = %0.3f, k^{2} = %0.2f, #sigma = %0.2f",stepsize,kappa_squared,sigma);
	l3->AddEntry(g3[3],legend,"l");

	double stepsize_new=10.*stepsize;
	int n_new=int((x1-x0)/stepsize_new+extrasteps+0.5);

	RK->change_parameters(n_new,stepsize_new,x0);		// Decreasing the precision of the solving
	EU->change_parameters(n_new,stepsize_new,x0);		// methods to see how the error changes.

	omega_squared=1; sigma=0; kappa_squared=1; gravity=gravity_earth;
	EQU->change_parameters(kappa_squared,gravity,sigma);

	RK->solver(omega_squared,EQU);
	g3[4] = new TGraph(n_new,RK->getposition(),RK->get_logerror_ex1());
	sprintf(legend,"RK4, h = %0.2f, k^{2} = %0.2f, #sigma = %0.2f",stepsize_new,kappa_squared,sigma);
	l3->AddEntry(g3[4],legend,"l");

	EU->solver(omega_squared,EQU);
	g3[5] = new TGraph(n_new,EU->getposition(),EU->get_logerror_ex1());
	sprintf(legend,"Euler, h = %0.2f, k^{2} = %0.2f, #sigma = %0.2f",stepsize_new,kappa_squared,sigma);
	l3->AddEntry(g3[5],legend,"l");

	g3[0]->GetXaxis()->SetNdivisions(5);
	g3[0]->GetYaxis()->SetNdivisions(5);
	g3[0]->SetTitle("");
	g3[0]->GetXaxis()->SetRangeUser(0.02,1.);
	g3[0]->GetYaxis()->SetRangeUser(-20.,5.);
	g3[0]->GetXaxis()->SetTitle("x");
	g3[0]->GetYaxis()->SetTitle("Log_{10}[#xi(x)-f(x)]");

	g3[0]->SetLineColor(2);
	g3[0]->Draw("AL");
	for(int z=1;z<6;z++){
		g3[z]->SetLineColor(2+z);
		g3[z]->Draw("same");
	}

	l3->Draw();
	c3->Modified();
	c3->Update();
	c3->SaveAs("errors.eps");

	c3->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");

//---------------------------------------- Sigma variation ------------------------------------------//

	TCanvas *c4 = new TCanvas("c4","Sigma variation",800,600);
	c4->SetGrid();
	TGraph *g4[5];
	TLegend *l4 = new TLegend(0.95,0.7,0.65,0.95);
	l4->SetTextAlign(12);

	eig_order_max=2; kappa_squared=1.; gravity=gravity_earth;

	SHO->change_parameters(n,stepsize,extrasteps,eig_order_max,x0,gravity);
	RK->change_parameters(n,stepsize,x0);
	EU->change_parameters(n,stepsize,x0);

	for(int i=0;i<5;i++){
		sigma=pow(i-1,4)+0.01;
		if (i==0)
			sigma = -19.;
		EQU->change_parameters(kappa_squared,gravity,sigma);

		SHO->shoot_now(EQU,false);
		RK->solver(SHO->get_eigenvalues(eig_order_max-1),EQU);
		g4[i] = new TGraph(n,RK->getposition(),RK->getsolution());
		sprintf(legend,"#omega_{%d}^{2} = %0.4f , #sigma = %0.2f",
				eig_order_max,SHO->get_eigenvalues(eig_order_max-1),sigma);
		l4->AddEntry(g4[i],legend,"l");
	}

	g4[0]->GetYaxis()->SetNdivisions(5);
	g4[0]->SetTitle(" ");
	g4[0]->GetXaxis()->SetRangeUser(0.,1.);
	g4[0]->GetYaxis()->SetRangeUser(-0.2,0.5);
	g4[0]->GetXaxis()->SetTitle("x");
	g4[0]->GetYaxis()->SetTitle("#xi(x)");

	g4[0]->SetLineColor(2);
	g4[0]->Draw("AL");
	for(int z=1;z<5;z++){
		g4[z]->SetLineColor(2+z);
		g4[z]->Draw("same");
	}

	l4->Draw();
	c4->Modified();
	c4->Update();
	c4->SaveAs("sigma.eps");

	c4->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");

//---------------------------------------------------------------------------------------------------//
//---------------------------------------- END OF PLOTTING ------------------------------------------//
//---------------------------------------------------------------------------------------------------//

	cout<<endl<<"Done! Elapsed time: "
		<<((double)clock()-time)/CLOCKS_PER_SEC<<"s"<<endl<<endl;	// Print the elapsed time.

	theApp.Run();													// Run Root Cern to show the plots.

	delete EQU, SHO, RK, EU;										// Deallocate memory.
	return 0;
}
