#ifndef __vector_h__
#define __vector_h__

class vector{
	public:
		vector()		{n_=0;v_=0;}
		vector(int);
		vector(const vector&);
		~vector();
		vector& operator=(const vector&);
		int getn()		{return n_;}
		void setx(int,double);
		double getx(int);
	protected:
		int n_;
		double* v_;
};

class linearvector: public vector{
	public:
		linearvector(int);
		~linearvector();
		linearvector operator+(linearvector&);
		double operator*(linearvector&);
		linearvector operator*(double);
};

#endif
