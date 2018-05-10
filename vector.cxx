#include <iostream>
#include <cstdlib>
#include "vector.h"
using namespace std;

vector::vector(int n){
	n_=n;
	v_=new double[n];
}

vector::vector(const vector& v){
	n_=v.n_;
	v_=new double[n_];
	for(int i=0;i<n_;i++)
		v_[i]=v.v_[i];
}

vector::~vector(){
	delete[] v_;
}

vector& vector::operator=(const vector& v){
	n_=v.n_;
	if(v_)
		delete[] v_;
	v_=new double[n_];
	for(int i=0;i<n_;i++)
		v_[i]=v.v_[i];
	return *this;
}

void vector::setx(int i,double x){
	if(i<n_)
		v_[i]=x;
	else{
		cout<<"Error!"<<endl;
		exit(-1);
	}
}

double vector::getx(int i){
	if(i<n_)
		return v_[i];
	else{
		cout<<"Error!"<<endl;
		exit(-1);
	}
}

linearvector::linearvector(int n):vector(n){}
linearvector::~linearvector(){}

linearvector linearvector::operator+ (linearvector& a){		// vectors sum
	linearvector v(getn());
	if(getn()!=a.getn()){
		cout<<"Error!"<<endl;
		exit(-1);
	}
	for(int i=0;i<getn();i++)
		v.setx(i,getx(i)+a.getx(i));
	return v;
}

double linearvector::operator*(linearvector& a){	// vectors cross product
	double v;
	for(int i=0;i<getn();i++)
		v+=getx(i)*a.getx(i);
	return v;
}

linearvector linearvector::operator*(double b){		// vectors scalar product
	linearvector v(getn());
	for(int i=0;i<getn();i++)
		v.setx(i,b*getx(i));
	return v;
}
