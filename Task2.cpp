#include "C:\Arsenal2019New\Arsenal.h"

double A=100,alf=0.05,T=50,beta=0;

double f(double x)
{
	return A*exp(-alf*x)*cos(pi2/T*x+beta);
}

Vect Power(int L,double x)
{
	int k;
	VSTACK(L,P);
	DO(k,1,L)P(k)=pow(x,k-1);
	return P;
}

//Полином Чебышева степени n
//для любого значения аргумента x
double PolCheb_0(int n,double x)
{
	if(n==0)return 1;
	if(n==1)return x;
	int i;
	double a=1,b=x,c;
	DO(i,2,n)
	{
		c=2*x*b-a;
		a=b; b=c;
	}
	return c;
}

//Вектор полиномов Чебышева.
//L - количество полиномов Чебышева;
//x - значение аргумента
Vect Cheb_0(int L,double x)
{
	int i;
	VSTACK(L,p);
	DO(i,1,L)p(i)=PolCheb_0(i-1,x);
	return p;
}

Vect& Approx_0(Vect Basis(int L,double x),int L,Vect& X,Vect& q)
{
	int k,n=X->m;
	Matr Fi(L,n);
	VSTACK(L,c);
	//DO(k,1,n)Fi(Basis(L,X(k)),k);
	c=~(Fi*!Fi)*Fi*q;
	return c;
}


Matr& Matr::operator*(Matr& B) //M17. Перемножение матриц
{
	if (n != B->m)
	{
		Ecrir(cout, "\n\n Matr::operator*(Matr& B): dimensionality of ");
		Ecrir(cout, "multiplicands are not coordinated!\n");
		Ecrir(cout, " -----> A->n=", n, "   B->m=", B->m); END
	}
	int l = B->n, i, j, k;
	STACKM(m, l)
	DO(i, 1, m)DO(j, 1, l)
		DO(k, 1, n)(*H)(i, j) += A[i - 1][k - 1] * B(k, j);
	return *H;
}

//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
int main()
{
	fstream out("Результаты.txt",ios::out|ios::in|ios::trunc);
	fstream graph("Для построения графиков.txt",ios::out|ios::trunc);
	cout.precision(15);	out.precision(15); graph.precision(15);

	int n,L=12,k,count=0,Ltab=30,O=12;
	double x,xmin=0,xmax=100,dx=5,fx,gx,Eerr=0,Eerr_=0;
	n=(xmax-xmin)/dx+1;
	Vect X(n),q(n),c(L), a;
	DO(k,1,n){X(k)=xmin+(k-1)*dx; q(k)=f(X(k));}

	c = Approx_0(Power,L,X,q);
	//c = Approx_0(Cheb_0,L,X,q);
	VecOut(cout,"c=",c,O);
	Do(x,xmin,xmax,dx/10)
	{
		fx=f(x); gx= c*Power(L,x);
		//fx=f(x); gx=c*Cheb_0(L,x);
		Eerr+=fabs(fx-gx); if(fabs(fx)>1e-6)Eerr_+=fabs((fx-gx)/fx);

		if(count%Ltab==0)
		{
			emp(2);
			Look(cout,"x","fx","gx","fx-gx","err_rel,%",O);
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
