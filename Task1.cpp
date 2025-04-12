#include "C:\Arsenal2019New\Arsenal.h"

double A=100,alf=0.05,T=50,beta=0;

double f(double x)
{
	return A*exp(-alf*x)*cos(pi2/T*x+beta);
}

double Lagrange_0(double f(double x),Vect& X,double x)
{
	int i,j,n=X->m;
	double S=0,P;
	DO(i,1,n)
	{
		P=1;
		DO(j,1,n)if(j!=i)P*=(x-X(j))/(X(i)-X(j));
		S+=f(X(i))*P;
	}
	return S;
}

//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
int main()
{
	fstream out("Результаты.txt",ios::out|ios::in|ios::trunc);
	fstream graph("Для построения графиков.txt",ios::out|ios::trunc);
	cout.precision(15);	out.precision(15); graph.precision(15);

	int n,k,count=0,Ltab=30,O=12;
	double x,xmin=0,xmax=100,dx=5,fx,gx,Eerr=0,Eerr_=0;
	n=(xmax-xmin)/dx+1;
	Vect X(n);
	DO(k,1,n)X(k)=xmin+(k-1)*dx;
	Do(x,xmin,xmax,dx/10)
	{
		fx=f(x); gx=Lagrange_0(f,X,x);
		Eerr+=fabs(fx-gx); if(fabs(fx)>1e-6)Eerr_+=fabs((fx-gx)/fx);

		if(count%Ltab==0)
		{
			emp(2);
			Look(cout,"t","fx","gx","fx-gx","err_rel,%",O);
		}
		emp();
		Look(cout,x,fx,gx,fx-gx,(fx-gx)/fx*100,O);

		Look(graph,x,fx,gx,fx-gx,(fx-gx)/fx*100,O); //Для графиков
		emp(graph);

		count++;
	}
	Eerr/=count; Eerr_=Eerr_/count*100;
	Ecrir(cout,"\n\n      Eerr=",Eerr,"  Eerr_,%=",Eerr_);

	Ecrir(cout,"\a");
	out.close(); graph.close();
	emp(5);
}
//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
