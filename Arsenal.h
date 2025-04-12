//-------------------------------------------------
//		 Соглашение о записи идентификаторов
//		------------------------------------- 
//Запись вида xOOy есть обозначение объединения
//(прямого сложения) векторов x и y.
//Запись вида d_f_dx обозначает дифференцирование
//функции f по переменной x
//-------------------------------------------------
#include <fstream.h>
#include <time.h>
#include <math.h>
#include <iomanip.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <stdarg.h>
#include <new.h>
#include <malloc.h>
#include <process.h>
#include <signal.h>
#include <dos.h>
#include <float.h>
#include <complex>
//#include <locate.h> //Для VS2019, VS2022

extern double pi,rdn,pi_4,pi_2,pi2;

extern double deltaStar1,alfaStar1, 
	   deltaStar2,alfaStar2,
	   deltaStar3,alfaStar3,
	   deltaStar4,alfaStar4,
	   deltaStar5,alfaStar5,
	   deltaStar6,alfaStar6,
	   deltaStar7,alfaStar7,
	   deltaStar8,alfaStar8,
	   deltaStar9,alfaStar9,
	   deltaStar10,alfaStar10,
	   deltaStar11,alfaStar11;

extern double rE_,aE_,bE_,acm_,mu_,muM_,muS_,omE_,Apr_,vC_,ro0_;

//xxxxxxxxxxxxxxxxxxxxxxxx     МАКРООПРЕДЕЛЕНИЯ И МАКРОСЫ     xxxxxxxxxxxxxxxxxxxxxxxx  
#define	MSTACK(m,n,C) static bool flag=1; \
					  static Matr* addr; Matr& C=*new Matr(m,n); \
					  if(flag){addr=&C; flag=0;} else{free((void*)addr); addr=&C;}

#define	VSTACK(m,x) static bool flag=1; \
					static Vect* addr; Vect& x=*new Vect(m); \
					if(flag){addr=&x; flag=0;} else{free((void*)addr); addr=&x;}

#define	STRSTACK(m,x) static bool flag=1; \
					  static Str* addr; Str &x=*new Str(m); \
					  if(flag){addr=&x; flag=0;}else{free((void*)addr); addr=&x;}

#define	CHARSTACK(m,x) static bool flag=1; \
					   static char* addr; char &x=*new char[m+1]; \
					   if(flag){addr=&x; flag=0;}else{free((void*)addr); addr=&x;}

#define DO(t,t1,t2) for(t=t1; t<=t2; t++)

#define OD(t,t1,t2) for(t=t1; t>=t2; t--)
 
#define Do(t,t1,t2,dt) if(t1<=t2)for(t=t1; t<=t2; t+=dt)

#define dO(t,t1,t2,dt) for(t=t1; t>=t2; t+=dt)

#define SQ_(x) (x)*(x)

#define EQ(x,x0) x>x0-1e-9 && x<x0+1e-9

#define END emp(3); exit(1);

#define MARKEND Ecrir(cout,"\n			MARKER_END"); emp(3); exit(1);

//Для StatCharPol (модификация входной функции на случай определения СХ её ОЛ)
#define FUNC(x,f) f=Func(x); if(type)f=f-fx0-dfdxT*(x-x0);
//xxxxxxxxxxxxxxxxxxxxxxxx  END  МАКРООПРЕДЕЛЕНИЯ И МАКРОСЫ  END  xxxxxxxxxxxxxxxxxxxxxxxx  

class Vect;	//Предварительное объявление класса для
	        //использования объектов типа Vect
	        //в определении класса Matr
            
class Matr
{
public:
	static int numbMatr,numbMatrMx;
	int m,n; double **A;
	Matr(){}	        //Конструктор без параметров
	Matr(int p,int q);	//Прототип "основного" конструктора
	~Matr();	        //Прототип деструктора
	void M(int p,int q,Matr& X); //Функция, заменяющая основной конструктор класса Matr,
						         //для объекта, созданного конструктором без параметров 
	Matr M(int p,int q); //Функция, заменяющая основной конструктор класса Matr
						 //для объекта, созданного конструктором без параметров
	double& operator()(int i,int j);
	Matr operator()(int i1,int i2,int j1,int j2);
	Vect operator()(int i1,int i2,int i);
	Vect operator()(int Ns);
	Matr& operator()(int k,int n0,Vect& a);
	Matr& operator()(Vect& a,int m0,int k);
	Vect operator[](int Nc);
	Matr operator()(Vect& s);
	Matr operator[](Vect& c);
	Matr operator()(Vect& s,Vect& c);
	Matr operator-(Matr& B);
	Matr operator+(Matr& B);
	Matr operator*(Matr& B);
	Vect operator*(Vect& b);
    Matr operator-(double q);
	Matr operator+(double q);
	Matr operator*(double q);
	Matr operator/(double q);
	Matr operator/(Matr& B);
	Matr operator!();
	Matr operator-();
	Matr operator~();
	Matr& operator=(double q);
	Matr& operator=(Matr& B);
	Matr& operator()(int k,Vect& a);
	Matr& operator()(Vect& a,int k);
	Matr& operator()(Matr& B,int m0,int n0);
	Matr *operator->(){return this;}
	void del();
	Matr RootMat();
	Matr Pow(double q);
	Matr& Dmat();
	Matr& Dmat(Vect& a);
	Matr& Dmat(double q);
	Matr Grevill();
	Matr GrevillNew();
	Matr SubMatr();
	Matr SubMatr(int rang,Vect& Am,Vect& An);
	Matr& TrCol(int i,int j);
	Matr& TrStr(int i,int j);
	bool Negativ();
	bool TestZero();
	Vect MatVecCol();
	Vect MatVecStr();
	Matr& VecMatCol(Vect& a);
	Matr& VecMatStr(Vect& a);
	Vect SimMatVec();
	Matr& VecSimMat(Vect& a);
	int Silvestr();
	Matr Cut(int I,int J);
	double Max();
	double Min();
};

//Прототипы глобальных операторов для класса Matr
Matr operator-(double q,Matr& A);
Matr operator+(double q,Matr& A);
Matr operator*(double q,Matr& A);
Matr operator/(double q,Matr& A);
Matr operator&(Matr& A,Matr& B);

class Vect
{																			
public:																
	static int numbVect,numbVectMx;			
	int m;
	double* a;
	Vect(){}	  //Конструктор без параметров			 
	Vect(int p);  //Прототип "основного" конструктора
	~Vect();      //Прототип деструктора
	void V(int p,Vect& x); //Функция, заменяющая основной конструктор класса Vect,
						   //для объекта, созданного конструктором без параметров 
	Vect V(int p); //Функция, заменяющая основной конструктор класса Vect,
				   //для объекта, созданного конструктором без параметров
	double& operator()(int i);										 
	Vect operator()(int i1,int i2);
	Vect operator()(Vect& b);
	Vect operator-(Vect& x);
	Vect operator+(Vect& x);
	Vect operator*(Matr& A);
	Matr operator^(Vect& b);
	double operator%(Vect& b);
	Vect operator-(double q);
	Vect operator+(double q);
	Vect operator*(double q);
	Vect operator/(double q);
	Vect operator/(Vect& b);
	double operator*(Vect& x);
	double operator!();
	Vect& operator=(double q);
	Vect& operator=(Vect& b);
	Vect& operator()(Vect& x,int m0); //хочешь так,
	Vect& operator()(int m0,Vect& x); //а хочешь этак
	Vect& operator()(int m0,double x);
	Vect& operator()(int m1,int m2,double x); //хочешь так,
	Vect& operator()(double x,int m1,int m2); //а хочешь этак
	Vect operator-();
	Vect *operator->(){return this;}
	Vect& Pow(double q);
	Vect& Dvec(Matr& A);
	bool TestZero();
	Vect& TrStr(int i,int j);
	double Min();
	double Max();
	void del();
};

//Прототипы глобальных операторов для класса Vect
Vect operator-(double q,Vect& x);
Vect operator+(double q,Vect& x);
Vect operator*(double q,Vect& x);
Vect operator/(double q,Vect& x);
Vect operator&(Vect& x,Vect& y);

class Tensor													
{																			
public:
	static int numbTensor,numbTensorMx;	
	int m1,m2,m3; double ***T;
	Tensor(){};						//Конструктор без параметров
	Tensor(int n1,int n2,int n3);	//Прототип "основного" конструктора
	~Tensor();			            //Прототип деструктора
	double& operator()(int k1,int k2,int k3);		
	Vect operator()(char *c,int k2,int k3);
	Vect operator()(int k1,char *c,int k3);
	Vect operator()(int k1,int k2,char *c);
	Tensor &operator()(Vect& x,int k2,int k3);
	Tensor &operator()(int k1,Vect& x,int k3);
	Tensor &operator()(int k1,int k2,Vect& x);
	Matr operator()(int k1,char *c2,char *c3);
	Matr operator()(char *c1,int k2,char *c3);
	Matr operator()(char *c1,char *c2,int k3);		
	Tensor &operator()(int k1,char *c2,char *c3,Matr& A);
	Tensor &operator()(char *c1,int k2,char *c3,Matr& A);
	Tensor &operator()(char *c1,char *c2,int k3,Matr& A);
	Tensor *operator->() {return this;}
	void del();	
};						
           
class Prob
{
public:
	int n;
	double rand,Level;
    Prob(){}								 //Конструктор без параметров
	Prob(double rand_);						 //Прототип конструктора
	Prob(int n_);							 //Прототип конструктора
	Prob(int n_,double rand_);				 //Прототип конструктора
	Prob(int n_,double rand_,double Level_); //Прототип конструктора
	~Prob(){}								 //Деструктор
	double Rrandom();						 //Датчик равн. распр. ПСЧ на [0,1]
	double Nrandom();						 //Датчик норм. распр. ПСЧ - N(0,1)
	Vect RandVec(Vect& x0,Matr& P);			 //Датчик норм. вектора N(Ex,Px)
	double Erf(double X);					 //Интеграл ошибок
	double ProbEll(int n,double r);			 //Вероятность попадания нормально распределенного
											 //n-мерного вектора X (достаточные статистики которого x0,P)
											 //										       T -1		  1/2
											 //в эллипсоид с приведенным радиусом r=[(x-x0) P  (x-x0)]
    double RadEll(int n_,double alfa);	     //Определение приведенного радиуса n-мерного
											 //эллипсоида, вероятностная мера которого равна alfa
    double Func(double ro);					 //Модуль отклонения вероятностной меры n-мерного эллипсоида
											 //с приведенным радиусом ro от заданного уровня Level
	double GoldSect(int n,double a,double b);
    Prob *operator->(){return this;}
};
            
class Filt
{
public:
	int static count1, //Номер шага прогноза (предсказания)
			   count2; //Номер шага обработки наблюдений
	int m,n,   //Размерности вектора наблюдения и вектора состояния
		Type_Line,	//Type_Line=0 - линеаризация ПФН без обновления g(x),H(x)
					//в фильтрах с последовательной обработкой компонент ВН,
					//Type_Line=1 - линеаризация ПФН с обновлением g(x)
					//в процессе обработки текущей компоненты ВН,
					//Type_Line=2 - линеаризация ПФН с обновлением g(x),H(x)
					//в процессе обработки текущей компоненты ВН
		N_Pxg, //Номер алгоритма определения априорной матрицы Pxg
		Nit, //Количество итераций в алгоритмах ОРБА
		Tu;  //Номер алгоритма определения u-параметра в ОРБА2
	bool Txs, //Тип учета ковариации ВС с ОЛ в алгоритмах ОРБА
		 Txg, //Тип учета ковариации ВС с ННФ в алгоритмах ОРБА
		 Tsg, //Тип учета ковариации ОЛ с ННФ в алгоритмах ОРБА
		 Control;  //Управление контролем достоверности оценок ВС
	double muFi,   //Масштаб шага для численного дифференцирования (прогноз)
		   muH,    //Масштаб шага для численного дифференцирования (коррекция)
		   Ro1,    //Масштаб шага для задания узлов интерполяции
				   //в алгоритмах определения статистик ОЛ
		   Ro2,    //Дополнительный масштаб шага для задания узлов интерполяции
				   //в алгоритме определения статистик ОЛ третьего порядка
		   uNit,   //Предельное число итераций
				   //для определения параметра u в ОРБА2
		   uEps,   //Точность определения параметра u в ОРБА2
				   //(при uEps>0 выход из процесса по значению функции,
				   //при uEps<0 выход из процесса по значению аргумента)
		   LevVal; //Уровень доверительной вероятности
	double m_omega,m_q0;
	Vect *m_psi,*m_q;
	Matr *m_Teta;
	Filt(){Init();};             //Конструктор без параметров
	Filt(int m_,int n_);         //Прототип конструктора
	Filt(int m_,int n_,int Tu_); //Прототип конструктора
	~Filt();					 //Прототип деструктора
	void Init();
	void Riccati(Vect& Deta,Matr& H,Matr& P);
	void Kalman(Vect& Dksi,Vect& fx0,Matr& Fi,Vect& x0,Vect& x,Matr& P);
	void Kalman(Vect& Deta,Vect& gx0,Matr& H,
				Vect& yR,Vect& x0,Vect& x,Matr& P);
	void Scalar_APR(Vect Func(Vect& x),Vect& Dksi,Vect& c,Matr& Pc,
					Vect& df,Matr& Pf,Vect& x0,Vect& x,Matr& P);
	void Scalar_APR(Vect Func(Vect& x),Matr& Fi,Vect& Dksi,Vect& c,Matr& Pc,
					Vect& df,Matr& Pf,Vect& x0,Vect& x,Matr& P);
	void Scalar_APS(Vect Func(Vect& x),Vect& Deta,Vect& s,Vect& dg,
					Matr& Ps,Matr& Pg,Vect& yR,Vect& x0,Vect& x,Matr& P);
	void Scalar_APS(Vect Func(Vect& x),Matr& H,Vect& Deta,Vect& s,Vect& dg,
					Matr& Ps,Matr& Pg,Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA_APR(Vect Func(Vect& x),Vect& Dksi,
				  Vect& c,Matr& Pc,Matr& Pxc,
				  Vect& df,Matr& Pf,Matr& Pxf,
				  Matr& Pcf,Vect& x0,Vect& x,Matr& P);
	void CRBA_APR(Vect Func(Vect& x),Matr& Fi,Vect& Dksi,
					Vect& c,Matr& Pc,Matr& Pxc,
				    Vect& df,Matr& Pf,Matr& Pxf,
					Matr& Pcf,Vect& x0,Vect& x,Matr& P);
	void CRBA1_APS(Vect Func(Vect& x),Vect& Deta,Vect& s,Matr& Ps,
				   Matr& Pxs,Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
				   Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA1_APS(Vect Func(Vect& x),Matr& H,
				   Vect& Deta,Vect& s,Matr& Ps,Matr& Pxs,
				   Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
				   Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA1_APS_IT(Vect Func(Vect& x),Vect& Deta,Vect& s,Matr& Ps,
					  Matr& Pxs,Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
					  Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA1_APS_IT(Vect Func(Vect& x),Matr& H,Vect& Deta,Vect& s,Matr& Ps,
					  Matr& Pxs,Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
					  Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA1(Vect Func(Vect& x),Matr Matrix(Vect& x),
			   Vect& Deta,Vect& s,Matr& Ps,Matr& Pxs,
			   Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
			   Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA1(Vect Func(Vect& x),Matr& H,Vect& Deta,Vect& s,Matr& Ps,
			   Matr& Pxs,Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
		       Vect& yR,Vect& x0,Vect& x,Matr& P);
	void CRBA2(Vect Func(Vect& x),Matr Matrix(Vect& x),
			   Matr& H,Vect& Deta,Vect& s,Matr& Ps,Matr& Pxs,
			   Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
			   Vect& yR,Vect& xB,Vect& x,Matr& P);
	void CRBA2(Vect Func(Vect& x),Matr Matrix(Vect& x),
			   Vect& Deta,Vect& s,Matr& Ps,Matr& Pxs,
			   Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
			   Vect& yR,Vect& xB,Vect& x,Matr& P);
	void CRBA2(Vect Func(Vect& x),Vect& Deta,Vect& s,Matr& Ps,
			   Matr& Pxs,Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
			   Vect& yR,Vect& xB,Vect& x,Matr& P);
	void CRBA_1_APS(Vect Func(Vect& x),bool Txs,bool Txg,bool Tsg,
							  Matr& H,Vect& Deta,Vect& s,Matr& Ps,Matr& Pxs,
							  Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,
							  Vect& yR,Vect& x0,Vect& x,Matr& P);
	Vect MinFunc(int Nit,double eps,Vect& dx,Vect& x0);
	void Positiv(double R,double Q,double dy,Vect& A0,Vect& B0,Vect& C0,
				 Vect& A,Vect& B,Vect& C,Matr& wPxs,Matr& wPxg,Matr& wPsg,
				 Vect& xB,Vect& x,Vect& s,Vect& dg,Matr& P,Matr& Ps,
				 Matr& Pxs,Matr& Pg,Matr& Pxg,Matr& Psg);
	void Risk(Vect Func(Vect& x),Vect& x0,Vect& dx,Vect& Sf);
	void StatCharFun(Vect Func(Vect& x),int L,Vect& x,Matr& P,Vect& f);
	void StatCharFun(Vect Func(Vect& x),Matr& dfdxT,
					 int L,Vect& x0,Vect& x,Matr& P,Vect& s);
	void StatCharFun(Vect Func(Vect& x),int N,Vect& x,
					 Matr& P,Vect& f,Matr& Pf);
	void StatCharFun(Vect Func(Vect& x),int N,Vect& x,
					 Matr& P,Vect& f,Matr& Pf,Matr& Pxf);
	void StatCharFun(Vect Func(Vect& x),Matr& dfdxT,int N,Vect& x0,
					 Vect& x,Matr& P,Vect& s,Matr& Ps);
	void StatCharFun(Vect Func(Vect& x),Matr& dfdxT,int N,Vect& x0,
					 Vect& x,Matr& P,Vect& s,Matr& Ps,Matr& Pxs);
	void TruncPol2(double Func(Vect& x),Vect& x0,Vect& dx,
				   Vect& A,Vect& B,Vect& C);
	Vect StatCharPol(Vect Func(Vect& x),Vect& x0,Vect& x,Matr& P);
	Vect StatCharPol(Vect& x0,Vect& x,Matr& P,Vect Func(Vect& x));
	void StatCharPol(Vect Func(Vect& x),Vect& x0,Vect& x,
					 Matr& P,Vect& f,Matr& Pf,Matr& Pxf);
	void StatCharPol(Vect Func(Vect& x),Vect& x0,
					 Vect& x,Matr& P,Vect& f,Matr& Pf);
	void StatCharPol(Vect Func(Vect& x),Matr& dfdxT,Vect& x0,
					 Vect& x,Matr& P,Vect& f,Matr& Pf,Matr& Pxf);
	void StatCharPol(Vect& x0,Vect& x,Matr& P,Vect& f,
					 Matr& Pf,Matr& Pxf,Vect Func(Vect& x));
	void StatCharPol(Vect& x0,Vect& x,Matr& P,
					 Vect& f,Matr& Pf,Vect Func(Vect& x));
	void StatCharPol(Matr& dfdxT,Vect& x0,Vect& x,Matr& P,
					 Vect& f,Matr& Pf,Matr& Pxf,Vect Func(Vect& x));
	void Muft(Vect& x0,Vect& x,Vect& fx0,Matr& A,Matr& B,Matr& C);
	void Muft(Vect& x0,Vect& x,Vect& fx0,Matr& A,Matr& B,
			  Matr& C,Matr& C1,Matr& C2,Matr& D,Matr& F);
	Matr SoftCovar(Matr& C,Matr& Px,Matr& Pz);
	Matr HardCovar(Matr& A,Matr& Px,Matr& Pz);
	Matr HardCovar(Matr& Pu,Matr& Pv);
	void Packing(Vect& x,Vect& s,Vect& dg,Vect& xCom);
	void Packing(Matr& P,Matr& Ps,Matr& Pxs,
				 Matr& Pg,Matr& Pxg,Matr& Psg,Matr& Pcom);
	Vect Param(int k,double Deta,double R,Vect& dy,Vect& A0,Vect& B0,
			   Vect& C0,Vect& A,Vect& B,Vect& C,Matr& P,Matr& Ps,
			   Matr& Pxs,Matr& Pg,Matr& Pxg,Matr& Psg,
			   Matr& H,Matr& Hk,double& Q,bool &bu);
	void Vcoef(double Deta,Vect& a0,Matr& a,
			   Vect& b0,Matr& b,Vect& C,Matr& Pcom);
	bool Vbest(Vect& v);
	double FuncTarget(Vect& v);
	bool Tester(Vect& xA,Vect& x,Matr& P);
	bool Tester(Vect& xA,Vect& x,Matr& PA,Matr& P);
	void See();
	Filt *operator->(){return this;}
};

class Matr;

class Str 
{						 
public:
	static int numbStr,numbStrMx;				 
	int m; char *s;														 
	Str(){};	//Конструктор без параметров
	Str(int p);	//Прототип "основного" конструктора
	~Str();	    //Прототип деструктора
	char &operator()(int ind);										 
	Str operator()(int ind1,int ind2);
	Str operator+(Str &x);									 
	Str operator+(char *x);										 
	Str &operator=(Str &x);
	Str &operator=(char *x);
	Str *operator->(){return this;}								 
	void del();
};

Str operator+(char *x,Str &y); //Прототип глобального оператора

class VStr	
{							
public:
	static int numbVStr,numbVStrMx;					
	int m,n; char **s;														
	VStr(){};	        //Конструктор без параметров
	VStr(int p,int q);	//Прототип "основного" конструктора
	~VStr();		    //Прототип деструктора
	Str operator()(int ind);										 
	char &operator()(int ind1,int ind2);					 
	Str operator()(int ind,int ind1,int ind2);	 
	VStr *operator->(){return this;}					 
	void del();	
};						

class MStr	
{							
public:
	static int numbMStr,numbMStrMx;
	int m1,m2,m3; char ***s;
	MStr(){};		            //Конструктор без параметров
	MStr(int n1,int n2,int n3); //Прототип "основного" конструктора
	~MStr();		            //Прототип деструктора
	Str operator()(int k1,int k2);								
	char &operator()(int k1,int k2,int k3);					
	Str operator()(int k1,int k2,int k3,int k4);	
	MStr *operator->(){return this;}						
	void del();									 
};														 

class TStr
{
public:
	static int numbTStr,numbTStrMx;
	int m1,m2,m3,m4; char ****s;
	TStr(){};                          //Конструктор без параметров:
	TStr(int n1,int n2,int n3,int n4);  //Прототип "основного" конструктора:
	~TStr();                            //Прототип деструктора;
	Str operator()(int k1,int k2,int k3);
	char &operator()(int k1,int k2,int k3,int k4);
	Str operator()(int k1,int k2,int k3,int k4,int k5);
	TStr *operator->(){return this;}
	void del();
};
           
class Ballist
{
public:
	bool Start,startMoon,startSun,startAtm;
	int Ngrav, //Порядок разложения ГПЗ в ряд по сферическим функциям
	    OnMoon,OnSun,OnLux,typeAtm, //Ключи учета возмущений
		kday,F0Atm,indexAtm,nday,day,month,year;
	double ud,		//Юлианская дата
		   BalCoef, //Баллистический коэффициент
		   A_dAtm,F107Atm,F81Atm,KpAtm,ud0Atm,
		   KFmoon[10][3],KFsun[11][3],
		   gamma, //Коэффициент  отражения (поглащения)
		   p, //Давление солнечной радиации вблизи Земли
		   Sm, //Площадь поперечного сечения КА
		   mass; //Масса КА
	Matr *Aatm,*Batm,*Catm,*Datm,*Eatm,*Fi1Atm,*Latm,*Natm;
	Ballist();
	Ballist(int N);
	~Ballist(); //Деструктор
	Vect ModProgn(double t,Vect& x);
	Vect ModProgn(Vect U(double t,Vect& x),double t,Vect& x);
	Vect sModProgn(double t,Vect& x);
	void GravJ2000(double t,Vect& c,Vect& d,Vect& x,Vect& dg);
	void GravGrinv(Vect& c,Vect& d,Vect& x,Vect& dg);
	Vect GravSphere(Vect& c,Vect& d,Vect& s);
	void SpeedMoon(double t,Vect& x,Vect& dg);
	void SpeedSun(double t,Vect& x,Vect& dg);
	void SpeedAtm(double t,Vect& x,Vect& dg);
	bool FactorLux(double t,Vect& x);
	void FactorLux(double t,Vect& x,double& dalfa);
	bool FactorLux(Vect& x,double t);
	void FactorLux(Vect& x,double t,double& deld);
	Vect PressLux(double t,Vect& x);
	double Atmos(double t,Vect& x);
	void Moon(double udt,Vect& r);
	void Sun(double udt,Vect& r);
	void Init();
	void Look();
	Ballist *operator->(){return this;}
};

//-----------------------------------------------------------------------------
//						Прототипы глобальных функций
Vect Abs(Vect& x);
Matr Abs(Matr& A);
Vect Adding(Vect& x,Vect& y);
Vect Adding(Vect& x,Vect& y,Vect& z);
Vect Adding(Vect& x,Matr& A);
Vect Adding(Matr& A,Vect& x);
Vect Adding(Vect& x, double s);
Vect Adding(double s, Vect& x);
Vect Adding(Vect& x,Matr& Lam,double c);
void Scatter(Vect& x,Vect& x1,Vect& x2);
void Scatter(Vect& x,Vect& x1,Vect& x2,Vect& x3);
void Scatter(Vect& xOOAOOc,Vect& x,Matr& A,double& c);
void Scatter(Vect& xOOAOOc,Vect& x,Matr& A);
void Scatter(Vect& AOOc,Matr& A,double& c);
double Angle(Vect& x,Vect& y);
Vect Approx(Vect Basis(int,double),int L,
			Vect& x,Vect& fx,int flag,double& eps);
Vect Approx(Vect Basis(int,double),int L,Vect& x,Vect& fx);
Vect Approx(double Bas(int,double),int L,Vect& T,Vect& q);
Matr Approx(Vect Basis(int,double),int L,Vect& Tb);
Matr Approx(int L,Vect& Tb);
Vect Approx(Vect Basis(int L,double t),int L,
			double ts,double p,Vect& t,Vect& q);
Vect Approx(Vect Basis(int L,double t),int L,int m,Vect& t,Vect& q);
Vect Approx(int L,int m,Vect& T,Vect& Q,Vect Basis(int L,double t));
double Approx(Vect Basis(int L,double t),int m,Vect& T,Vect& c,double t);
double Arcsin(double y,double x);
double Arctg(double y,double x);
double *Array(int m);
double **Array(int m,int n);
double ***Array(int m,int n,int l);
double Asin(double x,double y);
bool In_Paral(Vect& x,Vect& xmin,Vect& xmax);
bool Out_Paral(Vect& x,Vect& xmin,Vect& xmax);
//Определение нуля функции на отрезке. Дихотомия
double BiSec(double Func(double t),double t0,double tf,
			 double epsf,double epst,bool &kod);
//Определение нуля функции на отрезке. Дихотомия
bool BiSec(double Func(double t),double t0,double tf,
		   double epsf,double epst,double& t);
//Интегрирующий стартер для функции выхода
double Starter(Vect F(double t,Vect& x),double Fout(double t,Vect& x),
			   double eps,double& dt0,double t0,double t,Vect& x);
//Интерполирующий стартер для функции выхода
double Starter(double Fout(double t,Vect& x),
			   int L,double t,Vect& T,Matr& X,Vect& x);
//Обобщенная функция (интерполирование и интегрирование).
//Определением одной критической точки на [t0,tf] по функции
//выхода Fout(t,x). Функция f(t,x) - копия или упрощенный
//вариант функции F(t,x).
//Если eps>0, то для поиска КТ используется интегрирование,
//иначе - интерполирование
void BiSec(Vect F(double t,Vect& x),Vect f(double t,Vect& x),
		   double Fout(double t,Vect& x),double dt0,double t0,
		   double tf,double eps,double epsf,double epst,int L,
		   Vect& Tn,Matr& Xn,bool &kod,int &count,
		   Vect& x,double& T,Vect& X);
//Обобщенная функция (интерполирование и интегрирование).
//Определение всех критических точек на [t0,tf] по функции выхода
void BiSec(Vect F(double t,Vect& x),Vect f(double t,Vect& x),
		   double Fout(double t,Vect& x),bool &start,int L,
		   int Lout,double eps,double epsOut,double t0,
		   double tf,double dt,double dT,double& dt0,
		   Vect& x0,int &count,Vect& T0,Vect& Tf,
		   Matr& X0,Matr& Xf);
double Integral(double f(double x),double a,double b,double eps);
double Simpson(double f(double x),double a,double da,double fa,double fm,
			   double fb,int level,double absarea,double est,double eps);
//Интегрирование системы диф. уравнений
void BulSto(Vect F(double t,Vect& x),double eps,
			double& dt0,double t0,double tf,Vect& x);
void BulSto(Vect F(double t,Vect& x),double eps,
			double& dt0,double t0,double tf,Vect& s,Vect& x);
void BulSto(Vect F(double t,Vect& x),Vect& x,double eps,
			double& dt0,double t0,double tf);
void BulSto(Vect F(double t,Vect& x),double eps,double t0,double tf,Vect& x);
void BulSto(Vect F(double t,Vect& x),double FuncOut(double t,Vect& x),double eps,
			double epsOut,double& dt0,double t0,double tf,Vect& x,double& T,Vect& X);
void BulSto(Vect F(double t,Vect& x),double FuncOut(double t,Vect& x),
			double eps,double epsOut,double& dt0,double t0,
			double tf,int &count,Vect& x,Vect& T,Matr& X);
void BulSto(double FuncOut(double t,Vect& x),
			int L,double epsOut,double t0,double tf,
			int &count,Vect& Tnode,Matr& Xnode,
			Vect& x,Vect& T,Matr& X);
double BulSto(Vect F(double t,Vect& x),double FuncOut(double t,Vect& x),
			  double eps,double& dt0,double t0,double t,Vect& x);
double BulSto(double FuncOut(double t,Vect& x),
			  int L,double t,Vect& T,Matr& X,Vect& x);
Vect Polynom(int n,double x);
Vect Cheb(int n,double x);
int Cnm(int n,int m);
Matr Comb_2_m(int m);
void Comb(bool &prim,int n,int r,Vect& c);
Vect Comb(bool &prim,int n,int m);
Matr Comb(int n,int m);
void Comb(bool& prim,int& count,Vect& m,Vect& mm,Vect& a,Vect& b,Vect& c);
int Numb_Basis_Mult(int s,int n);
int Numb_All_Basis_Mult(int s,int n);
Matr Part_Basis_Mult(int s);
Matr All_Basis_Mult(int s,int n);
Matr Basis_Mult(int s,int n,int& m);
Vect Basis_Mult(double P(double dx,int k),Vect& dx,int s,int& m);
Vect All_Basis_Mult(double P(int k,double dx),int s,int L,Vect& dx);
double Polynom(double P(int k,double dx),int s,Vect& dx,Vect& A);
Vect M_N_K(double P(int k,double dx),
		   int s0,int L0,Vect& x0,Vect& xx,Vect& F);
Vect M_N_K(double P(int k,double dx),int s,int L,
		   Vect& x0,Vect& xx,Vect& F,Vect& X);
double Comp(double a,double b,double c);
double ConFr(Vect& x,Vect& y,double t);
Vect ConFr(Vect& x,Vect& y);
double CurSec();
Matr Cut(int I,Matr& P);
Matr Cut(Matr& P,int J);
Matr Cut(Matr& P,int I,int J);
Vect Cut(Vect& p,int k);
double dArcsin(double y,double x);
double dArctg(double y,double x);
double dAsin(double x,double y);
double dDist(Vect& q,Vect& r);
double dDist(Vect& q,Vect& r1,Vect& r2);
void DecBit(int x,char bit[]);
char *DecBit(int x);							
char *DecBase(int x,int base);
char *DecBase(double x,int base);
int NumberBase(int x,int base);
Vect DecSpher(Vect& d);
void Del(Matr& A);
void DerLegendr(int n,double x,Matr& dP);
double det(Matr& P);
Matr Except_one(int i,int j,Matr& A);
Matr df_dx(double f(Matr& x),Matr& x,Matr& dx);
Matr DifFun(Vect Func(Vect& x),int m,Vect& x,Vect& dx);
Matr Dif_Fun(Vect f(Vect& x),int m,Vect& dx,Vect& x);
Matr DifFunT(Vect Func(Vect& x),int m,Vect& x,Vect& dx);
Vect GradF(double f(Vect& x),Vect& x0,Vect& dx);
Vect Grad_Fun(double f(Vect& x),Vect& x,Vect& dx);
Vect DifVectScal(Vect& f(double t),int m,double t,double dt);
double DifFun(double Func(double x),double x,double dx);
double dnf_dxn(double f(double x),int n,double x,double dx);
double DivDif(int n,Vect& T,Vect& F,Vect& a);
void DivDif(Vect& T,Vect& F,Vect& a,Vect& b);
void DivDif(Vect& T,Vect& F,Vect& b);
void DivDif(double t,Vect& T,Vect& b,double& f,double& df);
Matr Dmat(Matr& A);
Matr Dmat(Vect& a);
Matr Dmat(int m,double a);
Vect Dvec(Matr& A);
Matr dTriang_dD(Matr& r,Vect& D,Vect& x0);
Matr dTriang_dx(double dx,Matr& r,Vect& D,Vect& x);
bool equ(double x,double x0);
bool equ(double x,double x0,double eps);
double Ephem(Vect& q,Vect& r,Vect& rV);
double Ephem(Vect& q,Vect& r,Vect& rb,Vect& rV,Vect& rVb);
Vect EphemErr(double ud,double t,Vect& Mvis,
			  Matr& xos,Vect& q,Matr& r,Matr& Pdr0);
double Erf(double X);
int Expon(double x);
long int Fact(int x);
double fact(int n);
double Forme(Vect& x,Matr& P,Vect& y);
double Forme(Vect& x,Matr& P);
Vect Fourier(Vect& f);
Vect Fourier(int k,double x);
void Equation(double a, double b, double c, double& x1, double& x2);
Vect Equation(double a, double b, double c);
double Pol(double f(Vect& x), Vect& x0, Vect& x, Vect& A, Vect& B, Vect& C);
void Polynom(double f(Vect& x), Vect& x0, Vect& dx, Vect& A, Vect& B, Vect& C);
double FRFunc(Vect& a,Vect& T,double t);
double Fifty_Fifty(double f(double x),int N,double a,double b);
double Fifty_Fifty(double f(double x),double eps,double a,double b);
double GoldSect(double Func(double x),double a,double b,int N);
double GoldSect(double Func(double t),double t0,double tf,double eps);
double GoldSect(double Func(double x),double a,double b,int N,double eps);
Matr Grev(Matr& B);
Vect Harmonic(int n,double t);
double IonSphere(Vect& q,Vect& r);
double IonSphere(Vect& q,Vect& r,Vect& rb);
void IonMinMax(double gammaMin,double& ionMin,double& ionMax);
Matr Grevill(Matr& A);
Matr GrevillNew(Matr& A);
Matr Inverse(Matr& C);
Matr Invert(Matr& C);
Matr Invers(Matr& A);
void Invers(Matr& A,Matr& B);
void Jacobi(Matr& B,Matr& S,Vect& lambda,double eps);
Vect Lagrange(Vect& T,double t);
double Lagrange(int L,double t,Vect& T,Vect& F);
double Legendre(int n,double x);
void Legendr(int n,double x,Matr& P);
void Legendr(int n,double x,Matr& P,Matr& dP);
Vect LinEqu(Matr& B,Vect& c);
Vect LinEqu(Matr& A,Vect& b,int nit,int &kod);
Matr MatDecSpher(double fi,double lam);
Matr MatDecSpher(Vect& x);
Matr MatrStatMom(int I,int J,int ni,int i1,int i2,int i3,
				 int j1,int j2,int j3,Matr& P);
Matr MatrStatMom(int I,int J,int i1,int i2,int i3,
				 int j1,int j2,int j3,Matr& P);

template <class Type> Type max(Type x1,Type x2){return x1> x2? x1: x2;}

template <class Type> Type max(Type x1,Type x2,Type x3)
{Type x=max(x1,x2); return max(x3,x);}

template <class Type> Type min(Type x1,Type x2,Type x3)
{Type x=min(x1,x2); return min(x3,x);}

template <class Type> Type min(Type x1,Type x2){return x1< x2? x1: x2;}

double Min(Vect& x);
double Max(Vect& x);

Matr MC(Vect& a);

void dy_Py(Vect& y,Vect& gx0,Vect& x0,Vect& x,Vect& s,Vect& dg,
		   Matr& H,Matr& P,Matr& Ps,Matr& Pxs,Matr& Pg,Matr& Pxg,
		   Matr& Psg,Matr& Deta,Vect& dy,Matr& Py);
void FAE_1(int k,Vect& xV,Vect& x,Vect& S0,double& Factor);
void FAE_2(int k,Vect& S0,Vect& S,double& Factor);
void FVE_1(int k,double ro_alfa,Vect& xV,Vect& x,Matr& P,double& Factor);
void FVE_2(int k,double ro_alfa,Vect& dy,Matr& Py,double& Factor);

void CRBA_APR(Vect Func(Vect& x),Matr& Fi,Vect& Dksi,Vect& c,Matr& Pc,
			  Matr& Pxc,Vect& df,Matr& Pf,Matr& Pxf,Vect& x0,Vect& x,Matr& P);
void CRBA_APS_1(Vect Func(Vect& x),bool Txs,bool Txg,bool Tsg, 
				Vect& Deta,Matr& H,Vect& s,Matr& Ps,Matr& Pxs,
				Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,Vect& yR,
				Vect& x0,Vect& x,Matr& P);
void CRBA1_APS_JOINT(Vect Func(Vect& x),bool Txs,bool Txg,bool Tsg, 
					 Matr& H,Vect& Deta,Vect& s,Matr& Ps,Matr& Pxs,
					 Vect& dg,Matr& Pg,Matr& Pxg,Matr& Psg,Vect& yR,
					 Vect& x0,Vect& x,Matr& P);	

int TestMinMax(double f(Vect& x),Vect g(Vect& x),
			   int N,int m,double delta,double eps,
			   double rand,Vect& x0);
Vect Metod_Prob(double f(Vect& x),Vect g(Vect& x),int N,
				int m,double& rand,Vect& x_min,Vect& x_max);
void Stochastic(double f(Vect& x),int N,double rand,
			    double& fmin,double& fmax,Vect& xmin,Vect& xmax);
Vect xmin_x_xmax(Vect& xmin,Vect& x,Vect& xmax);
Matr xmin_x_xmax(Matr& xmin,Matr& x,Matr& xmax);
Vect ProLinVar(Matr& A,Vect& b,Vect& p);
Vect ProLinVar(Vect& p,Vect& b,Matr& A);
double OptimStep(double f(Vect& x),double eps,
				 double hmin,double hmax,
				 Vect& x,Vect& s);
double OptimStepGold(double f(Vect& x),double eps,
				     double hmin,double hmax,
				     Vect& x,Vect& s);
double Min_fx_xk(double f(Vect& x),int k,double eps,
				 double xk_min,double xk_max,Vect& x);
double Min_fx_xk_Gold(double f(Vect& x),int k,double eps,
				      double xk_min,double xk_max,Vect& x);
int Inter_Direct(double f(Vect& x),Vect Gr(Vect& x),
				 double eps1,double eps2,
				 double hmin,double hmax,Vect& x);
void GradLow(double f(Vect& x),int N,double h,Vect& dx,Vect& x);
void GradLow(double f(Vect& x),int N,double eps,
			 double hmin,double hmax,Vect& dx,Vect& x);
void GradLow(double f(Vect& x),Vect Grad(Vect& x),
			 int N,double eps,double hmin,double hmax,Vect& x);
Vect Coord_Desc(double f(Vect& x),int& N,
				double eps1,double eps2,
				Vect& xmin,Vect& xmax);
void MePrGr(double f(Vect& x),int N,double eps,
			double hmin,double hmax,Vect& dx,
			Vect& xMin,Vect& xMax,Vect& x);
void MePrGr(double f(Vect& x),Vect& set,
			Vect& xMin,Vect& xMax,Vect& x);
void MePrGr(double f(Vect& x),Vect Grad(Vect& x),
			Vect& set,Vect& xMin,Vect& xMax,Vect& x);
void MePrGr(double f(Matr& X),int N,int M,
			double hmin,double hmax,double eps,
			Matr& dX,Matr& Xmin,Matr& Xmax,Matr& X);
void MePrGr(double f(Vect& x),Vect Grad(Vect& x),
			int N,int M,double hmin,double hmax,
			Vect& xMin,Vect& xMax,Vect& x);
void Packing(Matr& P,Matr& Ps,Matr& Pxs,
			 Matr& Pg,Matr& Pxg,Matr& Psg,Matr& Pcom);
Vect MinFunc(double Func(Vect& x),int nit,double eps,Vect& dx,Vect& x0);
Matr MinFunc(double Func(Matr& X),int Nit,double eps,Matr& dX,Matr& X0);
Vect MinFunc(double Func(Vect& x),int nit,double eps,Matr& xLim);

Matr MS(Vect& a);

void Muft(Vect& x0,Vect& x,Vect& fx0,Matr& A,Matr& B,Matr& C,Matr& Ac,Vect& q);
Matr Mult(Matr& A,Matr& B);
Vect Mult(Vect& a,Vect& b);
double Nrandom(double& x);
int Number(int n,double& rand,Vect& P);
int Number(int n,double& rand);
Vect Number(int m,int n,double& rand);
double NumbSect(double x,int m);
void Newton(double Func(Vect& x),Vect Grdnt(Vect& x),bool rgm,
			bool type_eps,int nit,double eps,Vect& dx,Vect& x);
void Perm(Vect& x);
Vect Perm(int n);
Vect Perm_(Vect& x);
double Harmon(int n,double t);
double PolPow(int n,double x);
double PolCheb(int L,double x);
double PolHerm(int m,double x);
double PolLag(int m,double x);
double PolLagw(int L,double x);
double PolLeg(int L,double x);
Vect BasisFunc(double P(int L,double t),int L,double t);
double Polynom(Vect& c,double x);
double mnPolynom(int m,int n,double x);
void PosSpots(double& rand,double R,double dL,Matr& Q);
void PosSpots(double& rand,double R,double dL,
			  double fimin,double fimax,Matr& Q);
Matr Pow(Matr& A,double q);
Vect Pow(Vect& a,double q);
double ProbEll(int n,double r);
double RadEll(int n,double alfa);
Vect RandVec(Vect& x0,Vect& sigma,double& rand);
Vect RandVec(double& rand,Vect& x0,Matr& P);
Vect RandVec(double& rand,Vect& xMin,Vect& xMax);
double BlancBruit(double& rand,double D,double dt);
Vect BlancBruit(bool& marker,double& rand,int N,
				Matr& D,double t,double t1,double t2);
int RangMatr(Matr& A,Vect& s,Vect& c,double eps);
int RangMatr(Matr& A,double eps);
//Частное нетривиальное решение
//линейной однородной системы уравнений
Vect NotTrivial(double eps,Matr& A);
//Определение собственных значений и собственных
//векторов матрицы A.
//Входные параметры:
//F		- функция, определяющая характеристическое
//		  уравнение det(A-lam*MatrixI(n))=0;
//Nit	- число итераций для решения характ. уравнения;
//eps	- наименьшая величина базисного минора,
//		  при которой он считается отличным от нуля
//A(n,n)- заданная матрица;
//Выходные параметры:
//lam	- вектор собственных значений.
//Возвращаемое значение:
//V(n,n)- матрица, i-й столбец которой - собствнный
//		  вектор, соответствующий lam(i), i=1,...,n
Matr OwnVector(double F(double lam),int Nit,double eps,
			   Matr& A,Vect& lam);
int RangMatr(double eps,Matr& C,Vect& row,Vect& col);
int RangMatr(double eps,Matr& C);
double Rrandom(double& x);
void Runge(double t,double tfin,Vect& x,Vect F(double t,Vect& x));
void Runge(Vect F(double t,Vect& x),double dt,double t,Vect& x);
void Runge(Vect F(double t,Vect& x),double eps,double t,double T,Vect& x);
void Runge(Matr F(double t,Matr& X),double t,double dt,Matr& X);
double Runge(Matr F(double t,Matr& X),double eps,
			 double t,double tfin,Matr& X);
double SigmMod(Vect& x,Matr& P);
double SigmMod(Vect& x,Matr& P,double& Emod);
double SigmMod(int N,Vect& x0,Matr& P,double& Edr);
double SigmMod(int N,Vect& x0,Matr& P);
int SizeFile(char *NameFile,int m);
int SizeFile(char *NameFile);
int SizeStr(char *NameFile);
Matr ContFile(char *Path);

template <class T> int sign(T x) 
{int sgn=0; if(x<0) sgn=-1; else if(x>0) sgn=1; return sgn;}

Vect SpherDec(Vect& s);
double Spline(Vect& T,Vect& F,double t);
double Spline(double t,Vect& T,Vect& F);
double Spline(int L,double t,Vect& T,Vect& F);
Matr Sqrt(Matr& A);
Vect Sqrt(Vect& a);
Vect Cos(Vect& x);
Vect Sin(Vect& x);
Matr& SortMatr(Matr& A);
Matr& SortMatr(int k,Matr& A);
void Sort1(Vect& x);
void Sort2(Vect& x);
Vect Sort3(Vect& x);
Vect Sort4(Vect& x);
Vect Sort5(Vect& x);
Vect Sort6(Vect& x);
Matr& Stack(Matr *addr,int m,int n);
Matr& Stack(int m,int n);
Vect& Stack(Vect *addr,int m);
Vect& Stack(int m);
//Интегрирующий стартер для функции выхода
double Starter(Vect F(double t,Vect& x),double FuncOut(double t,Vect& x),
			   double eps,double& dt0,double t0,double t,Vect& x);
double StartRand();
void StatCharFun(Vect Func(Vect& x),int L,Vect& x0,Vect& x,
				 Matr& P,Vect& s,Matr& Ps,Matr& Pxs);
void Statistics(int NumbTest,Vect& x,Vect& y,Vect& Ex,Vect& Ey,Matr& Pxy);
void Statistics(int NumbTest,double x,double y,double& Ex,double& Ey,double& Pxy);
void Statistics(int NumbTest,double x,double& Ex,double& Pxx);
double StatMom(int k1,int k2,int k3,int k4,Matr& P);
double StatMom(int k1,int k2,int k3,int k4,int k5,int k6,Matr& P);
void StMomErrApEx(double Basis(int,double),double t,Vect& Tb,
				  Matr& Eq,Matr& Sq,Matr& U,Vect& Ex,Vect& Sx);
void Submatrix(Matr& A,Matr& subA,Vect& Numb);
double Sum(Matr& A);
double Sum(Vect& a);
bool TestZero(Vect& x);
bool TestZero(Matr& A);
double Trace(Matr& A);
void TrEl(Vect& x,int i,int j);
Vect Triang(Matr& r,Vect& D,double& dD);
Vect Triang(Matr& r,Vect& D);
Vect Triang(Matr& r,Vect& D0,Vect& D,Vect& x0);
Vect Triang(Vect& D,Matr& r);
Vect Triang3(Matr& r,Vect& D);
Vect Triang5(Matr& r,Vect& D);
void Triangle(Matr& A,Matr& B,Matr& C);
double Unbalance(int M,int mA,int Ncor,double dL,Matr& dD,Matr& r,Vect& x);
Vect VN(int n);
Vect VS(int m,double s);
Vect VS(double s);
Vect V_ek(int n, int k);
Matr MatrixI(int n);
Vect VecProd(Vect& x,Vect& y);
Vect VerVS(Matr& q,Matr& r,Matr& rV);
Matr Zero(int m,int n);
Vect Zero(int m);
void ZeroFun(double F(double t),int Nit,double eps1,
			 double eps2,double eps3,double eta,Vect& c);
void ZeroFun(double F(double t),int Nit,double eps,Vect& c);

void emp();
void emp(int n);
char *Char(Str &x);
void Char(Str &x,char *s);						
void Copy(Str &s,VStr &vs,int k);
void Copy(char *s,VStr &vs,int k);
void Copy(char *s,VStr &vs,int k);
void Copy(char *s,MStr &Ms,int k1,int k2);
void Copy(Str &s,TStr &Ts,int k1,int k2,int k3);
void Copy(char *s,TStr &Ts,int k1,int k2,int k3);
char Lit(int x);
char *Liter(double x,int m);
char *Liter(char *c,int m);
void Liter(double x,int m,char *c);
void Liter(char *c,int m,char *s);
char *Line(int m);
char **Line(int m,int ns);
char ***Line(int m,int n,int ns);
char ****Line(int m,int n,int l,int ns);
void MatrEnt(fstream entry,Matr& A);
void MatrEnt(fstream entry,Vect& a);
void MatrEnt(fstream entry,Matr& A1,Matr& A2);
void MatrEnt(fstream entry,Vect& a1,Vect& a2);
void MatrEnt(fstream entry,Matr& A1,Matr& A2,Matr& A3);
void MatrEnt(fstream entry,Vect& a1,Vect& a2,Vect& a3);
void MatrEnt(fstream entry,Matr& A1,Matr& A2,Matr& A3,Matr& A4);
void MatrEnt(fstream entry,Vect& a1,Vect& a2,Vect& a3,Vect& a4);
void MatrEnt(fstream entry,Matr& A1,Matr& A2,Matr& A3,Matr& A4,Matr& A5);
void MatrEnt(fstream entry,Vect& a1,Vect& a2,Vect& a3,Vect& a4,Vect& a5);
char *Since(char *str1,char *str2);
int SizeStr(char *NameFile);

template <class TT,class T> void Read(TT entry,T &x) {entry>> x;}

template <class TT,class T1,class T2>
void Read(TT entry,T1 &x1,T2 &x2)
{entry>> x1>> x2;}

template <class TT,class T1,class T2,class T3>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3)
{entry>> x1>> x2>> x3;}

template <class TT,class T1,class T2,class T3,class T4>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4)
{entry>> x1>> x2>> x3>> x4;}

template <class TT,class T1,class T2,class T3,class T4,class T5>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5)
{entry>> x1>> x2>> x3>> x4>> x5;}

template <class TT,class T1,class T2,class T3,class T4,class T5,class T6>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6)
{entry>> x1>> x2>> x3>> x4>> x5>> x6;}

template <class TT,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,T7 &x7)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,class T10>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>> x10;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
		  class T10,class T11>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>> x10>> x11;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,T12 &x12)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>> x10>> x11>> x12;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
      class T10,class T11,class T12,class T13>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>> x10>> x11>> x12>> x13;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
   x10>> x11>> x12>> x13>> x14;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,class T15>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14,T15 &x15)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
   x10>> x11>> x12>> x13>> x14>> x15;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14,T15 &x15,T16 &x16)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
   x10>> x11>> x12>> x13>> x14>> x15>> x16;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
          class T15,class T16,class T17>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14,T15 &x15,T16 &x16,T17 &x17)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
   x10>> x11>> x12>> x13>> x14>> x15>> x16>> x17;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
          class T15,class T16,class T17,class T18>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14,T15 &x15,T16 &x16,T17 &x17,T18 &x18)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
  x10>> x11>> x12>> x13>> x14>> x15>> x16>> x17>> x18;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
          class T15,class T16,class T17,class T18,class T19>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14,T15 &x15,T16 &x16,
          T17 &x17,T18 &x18,T19 &x19)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
   x10>> x11>> x12>> x13>> x14>> x15>> x16>> x17>> x18>> x19;}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
          class T15,class T16,class T17,class T18,class T19,class T20>
void Read(TT entry,T1 &x1,T2 &x2,T3 &x3,T4 &x4,
		  T5 &x5,T6 &x6,T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
          T12 &x12,T13 &x13,T14 &x14,T15 &x15,T16 &x16,
          T17 &x17,T18 &x18,T19 &x19,T19 &x20)
{entry>> x1>> x2>> x3>> x4>> x5>> x6>> x7>> x8>> x9>>
   x10>> x11>> x12>> x13>> x14>> x15>> x16>> x17>> x18>> x19>> x20;}

//Вывод двумерного массива в поток out,
//|size|>=6 - число позиций, выделяемых под число.
//Если size<0, то имя массива указывается с каждым элементом.
//Иначе имя массива указывается один раз перед всеми его элементами,
//и после записи элементов точка с запятой не пишется.
//ncol - количество столбцов в таблице
template <class TT,class S> void MatrOut(TT out,S s,Matr& A,int size,int ncol)
{
	int Size=abs(size),m=A.m,n=A.n,i,j0,j,k,kmx,r,M,mn,ij,d;
	//ncol=7; if(Size<=6)ncol=15; if(Size>6 && Size<=9)ncol=10;
	if(Size>9 && Size<=12)ncol=8; M=min(n,ncol);
	kmx=n/M; r=n%M;
	//if(size>0)out<< "\n "<< s<< ":\n"; else out<< "\n";
	if(size>0)out<< "\n "<< s<< "\n"; else out<< "\n";
	mn=(int)(log10(m)+1)+(int)(log10(n)+1);
	for(k=1; k<=kmx; k++)
	{
		j0=(k-1)*M;
		for(i=1; i<=m; i++)
		{
			for(j=j0+1; j<=M+j0; j++)
			{
				ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-2;
				char *a=Line(Size+d); Liter(A(i,j),Size+d,a);
				if(size>0)out<< " ("<< i<<","<< j<< ")="<< a<< " ";
				if(size<0)out<< " "<< s<< "("<< i<<","<< j<< ")="<< a<< "; ";
				delete []a;
			}
			out<< "\n";
		}
		if(k<kmx || r>0)out<< "\n";
	}
	if(r==0)return;
	for(i=1; i<=m; i++)
	{
		for(k=j; k<j+r; k++)
		{
			ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-2;
			char *a=Line(Size+d); Liter(A(i,k),Size+d,a);
			if(size>0)out<< " ("<< i<<","<< k<< ")="<< a<< " ";
			if(size<0)out<< " "<< s<< "("<< i<<","<< k<< ")="<< a<< " ";
			delete []a;
		}
		out<< "\n";
	}
}

//Вывод двумерного массива в поток out,
//|size|>=6 - число позиций, выделяемых под число.
//Если size<0, то имя массива указывается с каждым элементом.
//Иначе имя массива указывается один раз перед всеми его элементами,
//и после записи элементов точка с запятой не пишется
template <class TT,class S> void MatrOut(TT out,S s,Matr& A,int size)
{
	int Size=abs(size),m=A.m,n=A.n,i,j0,j,k,kmx,r,ncol,M,mn,ij,d;
	ncol=7; if(Size<=6)ncol=15; if(Size>6 && Size<=9)ncol=10;
	if(Size>9 && Size<=12)ncol=8; M=min(n,ncol);
	kmx=n/M; r=n%M;
	//if(size>0)out<< "\n "<< s<< ":\n"; else out<< "\n";
	if(size>0)out<< "\n "<< s<< "\n"; else out<< "\n";
	mn=(int)(log10(m)+1)+(int)(log10(n)+1);
	for(k=1; k<=kmx; k++)
	{
		j0=(k-1)*M;
		for(i=1; i<=m; i++)
		{
			for(j=j0+1; j<=M+j0; j++)
			{
				ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-2;
				char *a=Line(Size+d); Liter(A(i,j),Size+d,a);
				if(size>0)out<< " ("<< i<<","<< j<< ")="<< a<< " ";
				if(size<0)out<< " "<< s<< "("<< i<<","<< j<< ")="<< a<< "; ";
				delete []a;
			}
			out<< "\n";
		}
		if(k<kmx || r>0)out<< "\n";
	}
	if(r==0)return;
	for(i=1; i<=m; i++)
	{
		for(k=j; k<j+r; k++)
		{
			ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-2;
			char *a=Line(Size+d); Liter(A(i,k),Size+d,a);
			if(size>0)out<< " ("<< i<<","<< k<< ")="<< a<< " ";
			if(size<0)out<< " "<< s<< "("<< i<<","<< k<< ")="<< a<< " ";
			delete []a;
		}
		out<< "\n";
	}
}

//Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT,class S>
void MatrOut(TT out,S s1,Matr& A1,S s2,Matr& A2,int size)
{
	MatrOut(out,s1,A1,size);
	MatrOut(out,s2,A2,size);
}

//Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT,class S>
void MatrOut(TT out,S s1,Matr& A1,S s2,Matr& A2,
			 S s3,Matr& A3,int size)
{
	MatrOut(out,s1,A1,size);
	MatrOut(out,s2,A2,size);
	MatrOut(out,s3,A3,size);
}

//Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT,class S>
void MatrOut(TT out,S s1,Matr& A1,S s2,Matr& A2,
			 S s3,Matr& A3,S s4,Matr& A4,int size)
{
	MatrOut(out,s1,A1,size);
	MatrOut(out,s2,A2,size);
	MatrOut(out,s3,A3,size);
	MatrOut(out,s4,A4,size);
}

//Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT,class S>
void MatrOut(TT out,S s1,Matr& A1,S s2,Matr& A2,
			 S s3,Matr& A3,S s4,Matr& A4,S s5,Matr& A5,int size)
{
	MatrOut(out,s1,A1,size);
	MatrOut(out,s2,A2,size);
	MatrOut(out,s3,A3,size);
	MatrOut(out,s4,A4,size);
	MatrOut(out,s5,A5,size);
}

//Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT,class S>
void MatrOut(TT out,S s1,Matr& A1,S s2,Matr& A2,
			 S s3,Matr& A3,S s4,Matr& A4,S s5,Matr& A5,
			 S s6,Matr& A6,int size)
{
	MatrOut(out,s1,A1,size);
	MatrOut(out,s2,A2,size);
	MatrOut(out,s3,A3,size);
	MatrOut(out,s4,A4,size);
	MatrOut(out,s5,A5,size);
	MatrOut(out,s6,A6,size);
}

//Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT,class S>
void MatrOut(TT out,S s1,Matr& A1,S s2,Matr& A2,
			 S s3,Matr& A3,S s4,Matr& A4,S s5,Matr& A5,
			 S s6,Matr& A6,S s7,Matr& A7,int size)
{
	MatrOut(out,s1,A1,size);
	MatrOut(out,s2,A2,size);
	MatrOut(out,s3,A3,size);
	MatrOut(out,s4,A4,size);
	MatrOut(out,s5,A5,size);
	MatrOut(out,s6,A6,size);
	MatrOut(out,s7,A7,size);
}

//Вывод двумерного массива в поток out,
//|size|>=6 - число позиций, выделяемых под число.
//Если size<0, то имя массива указывается с каждым элементом,
//и после каждого элемента пишется точка с запятой.
//Иначе имя массива указывается один раз перед всеми его элементами,
//и после записи элементов точка с запятой не пишется
template <class TT> void MatrOut(TT out,Matr& A,int size)
{
	int Size=abs(size),m=A.m,n=A.n,i,j0,j,k,
		kmx,r,ncol=9,M=min(n,ncol),mn,ij,d;
	kmx=n/M; r=n%M;	mn=(int)(log10(m)+1)+(int)(log10(n)+1);
	for(k=1; k<=kmx; k++)
	{
		j0=(k-1)*M;
		for(i=1; i<=m; i++)
		{
			for(j=j0+1; j<=M+j0; j++)
			{
				ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-3;
				char *a=Line(Size+d+2); Liter(A(i,j),Size+d+2,a); //04.04.20//
				if(size>0)out<< " ("<< i<<","<< j<< ")="<< a<< " ";
				if(size<0)out<< " "<< "("<< i<<","<< j<< ")="<< a<< "; ";
				delete []a;
			}
			out<< "\n";
		}
		if(k<kmx || r>0)out<< "\n";
	}
	if(r==0)return;
	for(i=1; i<=m; i++)
	{
		for(k=j; k<j+r; k++)
		{
			ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-3;
			char *a=Line(Size+d+2); Liter(A(i,j),Size+d+2,a); //04.04.20//
			if(size>0)out<< " ("<< i<<","<< k<< ")="<< a<< " ";
			if(size<0)out<< " "<< "("<< i<<","<< k<< ")="<< a<< " ";
			delete []a;
		}
		out<< "\n";
	}
}

//Вывод двумерного массива в поток out без индексов,
//size>=6 - число позиций, выделяемых под число.
template <class TT> void MatrOut(TT out,int size,Matr& A)
{
	int m=A->m,n=A->n,i,j0,j,k,kmx,r,ncol=10,M=min(n,ncol),mn,ij,d;
	kmx=n/M; r=n%M;	mn=(int)(log10(m)+1)+(int)(log10(n)+1);
	for(k=1; k<=kmx; k++)
	{
		j0=(k-1)*M;
		for(i=1; i<=m; i++)
		{
			for(j=j0+1; j<=M+j0; j++)
			{
				if(i==1 && j==j0+1)ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-3;
				char *a=Line(size+d+2); Liter(A(i,j),size+d+2,a); //04.04.20//
				out<< a<< " "; delete []a;
			}
			out<< "\n";
		}
		if(k<kmx || r>0)out<< "\n";
	}
	if(r==0)return;
	for(i=1; i<=m; i++)
	{
		for(k=j; k<j+r; k++)
		{
			if(i==1 && k==j)ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-3;
			char *a=Line(size+d+2); Liter(A(i,j),size+d+2,a); //04.04.20//
			out<< a<< " "; delete []a;
		}
		out<< "\n";
	}
}

//Вывод двумерного массива в поток out без индексов,
//size>=6 - число позиций, выделяемых под число.
template <class TT> void MatrOut(TT out,int size,int ncol,Matr& A)
{
	int m=A.m,n=A.n,i,j0,j,k,kmx,r,M=min(n,ncol),mn,ij,d;
	kmx=n/M; r=n%M; mn=(int)(log10(m)+1)+(int)(log10(n)+1);
	for(k=1; k<=kmx; k++)
	{
		j0=(k-1)*M;
		for(i=1; i<=m; i++)
		{
			for(j=j0+1; j<=M+j0; j++)
			{
				if(i==1 && j==j0+1)ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-3;
				char *a=Line(size+d+2); Liter(A(i,j),size+d+2,a); //04.04.20//
				out<< a<< " "; delete []a;
			}
			out<< "\n";
		}
		if(k<kmx || r>0)out<< "\n";
	}
	if(r==0)return;
	for(i=1; i<=m; i++)
	{
		for(k=j; k<j+r; k++)
		{
			if(i==1 && k==j)ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij-3;
			char *a=Line(size+d+2); Liter(A(i,j),size+d+2,a); //04.04.20//
			out<< a<< " "; delete []a;
		}
		out<< "\n";
	}
}

//Вывод трехмерного массива в поток out,
//|size|>=6 - число позиций, выделяемых под число.
//Массив выводится слоями - сечениями по третьему измерению.
//Если size<0, то имя массива указывается с каждым элементом,
//и после каждого элемента пишется точка с запятой.
//Иначе имя массива указывается один раз перед всеми его элементами,
//и после записи элементов точка с запятой не пишется
template <class TT,class S> void TensOut(TT out,S s,Tensor &T,int size)
{
	int Size=abs(size),m=T.m1,n=T.m2,m3=T.m3,i,j0,j,k,l,
		kmx,r,ncol=9,M=min(n,ncol),mn,ij,d;
	Matr A(m,n);
	kmx=n/M; r=n%M;
	//if(size>0)out<< "\n "<< s<< ":\n"; else out<< "\n";
	if(size>0)out<< "\n "<< s<< "\n"; else out<< "\n";
	mn=(int)(log10(m)+1)+(int)(log10(n)+1);
	DO(l,1,m3)
	{
		A=T("","",l);
		for(k=1; k<=kmx; k++)
		{
			j0=(k-1)*M;
			for(i=1; i<=m; i++)
			{
				for(j=j0+1; j<=M+j0; j++)
				{
					ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij+1;
					char *a=Line(Size+d+2); Liter(A(i,j),Size+d+2,a); //05.04.20//
					if(size>0)out<< " ("<< i<<","<< j<<","<< l<< ")="<< a<< " ";
					if(size<0)out<< " "<< s<< "("<< i<<","<< j<<","<< l<< ")="<< a<< "; ";
					delete []a;
				}
				out<< "\n";
			}
			if(k<kmx || r>0)out<< "\n";
		}
		if(r==0 && l==m3)return;
		for(i=1; i<=m; i++)
		{
			for(k=j; k<j+r; k++)
			{
				ij=(int)(log10(i)+1)+(int)(log10(j)+1); d=mn-ij+1;
				char *a=Line(Size+d+2); Liter(A(i,k),Size+d,a+2); //05.04.20//
				if(size>0)out<< " ("<< i<<","<< k<<","<< l<< ")="<< a<< " ";
				if(size<0)out<< " "<< s<< "("<< i<<","<< k<<","<< l<< ")="<< a<< " ";
				delete []a;
			}
		}
		out<< "\n";
	}
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s,T &x,int q)
{
	int i=0,m=x.m,ls=strlen(s),N=10; char *gap=Line(ls+2);
	//if(s!="")output<< "\n "<< s<< "="; else output<< "\n";
	if(s!="")output<< "\n "<< s; else output<< "\n";
	for(i=1; i<=m; i++)
	{
		char *v=Line(q); Liter(x(i),q,v);
		output<< v<< " "; if(fmod(i,N)==0){output<< "\n"<< gap;}
		delete []v;
	}
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,int q)
{
	VecOut(output,s1,x1,q);
	VecOut(output,s2,x2,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												S s5,T &x5,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
	VecOut(output,s5,x5,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												S s5,T &x5,
												S s6,T &x6,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
	VecOut(output,s5,x5,q); VecOut(output,s6,x6,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												S s5,T &x5,
												S s6,T &x6,
												S s7,T &x7,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
	VecOut(output,s5,x5,q); VecOut(output,s6,x6,q);
	VecOut(output,s7,x7,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												S s5,T &x5,
												S s6,T &x6,
												S s7,T &x7,
												S s8,T &x8,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
	VecOut(output,s5,x5,q); VecOut(output,s6,x6,q);
	VecOut(output,s7,x7,q); VecOut(output,s8,x8,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												S s5,T &x5,
												S s6,T &x6,
												S s7,T &x7,
												S s8,T &x8,
												S s9,T &x9,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
	VecOut(output,s5,x5,q); VecOut(output,s6,x6,q);
	VecOut(output,s7,x7,q); VecOut(output,s8,x8,q);
	VecOut(output,s9,x9,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class S,class T> void VecOut(TT output,S s1,T &x1,
												S s2,T &x2,
												S s3,T &x3,
												S s4,T &x4,
												S s5,T &x5,
												S s6,T &x6,
												S s7,T &x7,
												S s8,T &x8,
												S s9,T &x9,
												S s10,T &x10,
												int q)
{
	VecOut(output,s1,x1,q); VecOut(output,s2,x2,q);
	VecOut(output,s3,x3,q); VecOut(output,s4,x4,q);
	VecOut(output,s5,x5,q); VecOut(output,s6,x6,q);
	VecOut(output,s7,x7,q); VecOut(output,s8,x8,q);
	VecOut(output,s9,x9,q); VecOut(output,s10,x10,q);
}

//Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT,class T> void VectOut(TT output,T &x,int q)
{
	int i,j,m=x.m; emp(output);
	for(i=1; i<=m; i++)
	{
		char *v=Line(q); if(i<=m)Liter(x(i),q,v);
		for(j=0; j<q; j++)output<< v[j]; output<< " ";
		if(i<m)output<< "\n";
		delete []v;
	}	
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2>
void VectOut(TT output,T1 &x1,T2 &x2,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m;
	m=m1>m2? m1: m2; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		for(j=0; j<q; j++)output<< s1[j]; output<< " ";
		for(j=0; j<q; j++)output<< s2[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),*s3=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3; m=m>m4? m: m4; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),*s3=Line(q),*s4=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,m5=x5.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3; m=m>m4? m: m4; m=m>m5? m: m5; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),*s3=Line(q),*s4=Line(q),*s5=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4; delete []s5;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,class T6>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),*s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,m7=x7.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),*s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),*s7=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,m7=x7.m,m8=x8.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,m7=x7.m,m8=x8.m,m9=x9.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,
		m5=x5.m,m6=x6.m,m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,m5=x5.m,
		m6=x6.m,m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m,m11=x11.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,T12 &x12,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m,m11=x11.m,m12=x12.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
		 *s11=Line(q),*s12=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m,m11=x11.m,m12=x12.m,m13=x13.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m,m11=x11.m,m12=x12.m,
		m13=x13.m,m14=x14.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,T15 &x15,int q)
{
	int i=0,j=0,m=0,
		m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m,m11=x11.m,m12=x12.m,
		m13=x13.m,m14=x14.m,m15=x15.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; m=m>m15? m: m15; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q),
			 *s15=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		if(i<=m15)Liter(x15(i),q,s15); else for(j=0; j<q; j++)s15[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j]; output<< " ";
		for(j=0; j<q; j++) output<< s15[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14; delete []s15;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,T15 &x15,
			 T16 &x16,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,m4=x4.m,
		m5=x5.m,m6=x6.m,m7=x7.m,m8=x8.m,m9=x9.m,m10=x10.m,
		m11=x11.m,m12=x12.m,m13=x13.m,m14=x14.m,m15=x15.m,m16=x16.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; m=m>m15? m: m15; m=m>m16? m: m16; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q),
			 *s15=Line(q),*s16=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		if(i<=m15)Liter(x15(i),q,s15); else for(j=0; j<q; j++)s15[j]=' ';
		if(i<=m16)Liter(x16(i),q,s16); else for(j=0; j<q; j++)s16[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j]; output<< " ";
		for(j=0; j<q; j++) output<< s15[j]; output<< " ";
		for(j=0; j<q; j++) output<< s16[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14; delete []s15; delete []s16;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,class T17>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,T15 &x15,
			 T16 &x16,T17 &x17,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,
		m10=x10.m,m11=x11.m,m12=x12.m,
		m13=x13.m,m14=x14.m,m15=x15.m,
		m16=x16.m,m17=x17.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; m=m>m15? m: m15; m=m>m16? m: m16;
	m=m>m17? m: m17; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q),
			 *s15=Line(q),*s16=Line(q),
			 *s17=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		if(i<=m15)Liter(x15(i),q,s15); else for(j=0; j<q; j++)s15[j]=' ';
		if(i<=m16)Liter(x16(i),q,s16); else for(j=0; j<q; j++)s16[j]=' ';
		if(i<=m17)Liter(x17(i),q,s17); else for(j=0; j<q; j++)s17[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j]; output<< " ";
		for(j=0; j<q; j++) output<< s15[j]; output<< " ";
		for(j=0; j<q; j++) output<< s16[j]; output<< " ";
		for(j=0; j<q; j++) output<< s17[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14; delete []s15; delete []s16;
		delete []s17;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,class T17,
		  class T18>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,T15 &x15,
			 T16 &x16,T17 &x17,T18 &x18,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,
		m10=x10.m,m11=x11.m,m12=x12.m,
		m13=x13.m,m14=x14.m,m15=x15.m,
		m16=x16.m,m17=x17.m,m18=x18.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; m=m>m15? m: m15; m=m>m16? m: m16;
	m=m>m17? m: m17; m=m>m18? m: m18; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q),
			 *s15=Line(q),*s16=Line(q),
			 *s17=Line(q),*s18=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		if(i<=m15)Liter(x15(i),q,s15); else for(j=0; j<q; j++)s15[j]=' ';
		if(i<=m16)Liter(x16(i),q,s16); else for(j=0; j<q; j++)s16[j]=' ';
		if(i<=m17)Liter(x17(i),q,s17); else for(j=0; j<q; j++)s17[j]=' ';
		if(i<=m18)Liter(x18(i),q,s18); else for(j=0; j<q; j++)s18[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j]; output<< " ";
		for(j=0; j<q; j++) output<< s15[j]; output<< " ";
		for(j=0; j<q; j++) output<< s16[j]; output<< " ";
		for(j=0; j<q; j++) output<< s17[j]; output<< " ";
		for(j=0; j<q; j++) output<< s18[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14; delete []s15; delete []s16;
		delete []s17; delete []s18;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,class T17,
		  class T18,class T19>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,T15 &x15,
			 T16 &x16,T17 &x17,T18 &x18,T19 &x19,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,
		m10=x10.m,m11=x11.m,m12=x12.m,
		m13=x13.m,m14=x14.m,m15=x15.m,
		m16=x16.m,m17=x17.m,m18=x18.m,
		m19=x19.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; m=m>m15? m: m15; m=m>m16? m: m16;
	m=m>m17? m: m17; m=m>m18? m: m18; m=m>m19? m: m19;
	emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q),
			 *s15=Line(q),*s16=Line(q),
			 *s17=Line(q),*s18=Line(q),
			 *s19=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		if(i<=m15)Liter(x15(i),q,s15); else for(j=0; j<q; j++)s15[j]=' ';
		if(i<=m16)Liter(x16(i),q,s16); else for(j=0; j<q; j++)s16[j]=' ';
		if(i<=m17)Liter(x17(i),q,s17); else for(j=0; j<q; j++)s17[j]=' ';
		if(i<=m18)Liter(x18(i),q,s18); else for(j=0; j<q; j++)s18[j]=' ';
		if(i<=m19)Liter(x19(i),q,s19); else for(j=0; j<q; j++)s19[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j]; output<< " ";
		for(j=0; j<q; j++) output<< s15[j]; output<< " ";
		for(j=0; j<q; j++) output<< s16[j]; output<< " ";
		for(j=0; j<q; j++) output<< s17[j]; output<< " ";
		for(j=0; j<q; j++) output<< s18[j]; output<< " ";
		for(j=0; j<q; j++) output<< s19[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14; delete []s15; delete []s16;
		delete []s17; delete []s18; delete []s19;
	}
}

//Вывод значений векторов в поток output с разрядностью q
template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,class T17,
		  class T18,class T19,class T20>
void VectOut(TT output,T1 &x1,T2 &x2,T3 &x3,T4 &x4,T5 &x5,T6 &x6,
			 T7 &x7,T8 &x8,T9 &x9,T10 &x10,T11 &x11,
			 T12 &x12,T13 &x13,T14 &x14,T15 &x15,
			 T16 &x16,T17 &x17,T18 &x18,T19 &x19,
			 T20 &x20,int q)
{
	int i=0,j=0,m=0,m1=x1.m,m2=x2.m,m3=x3.m,
		m4=x4.m,m5=x5.m,m6=x6.m,
		m7=x7.m,m8=x8.m,m9=x9.m,
		m10=x10.m,m11=x11.m,m12=x12.m,
		m13=x13.m,m14=x14.m,m15=x15.m,
		m16=x16.m,m17=x17.m,m18=x18.m,
		m19=x19.m,m20=x20.m;
	m=m1>m2? m1: m2; m=m>m3? m: m3;m=m>m4? m: m4;
	m=m>m5? m: m5; m=m>m6? m: m6; m=m>m7? m: m7;
	m=m>m8? m: m8; m=m>m9? m: m9; m=m>m10? m: m10;
	m=m>m11? m: m11; m=m>m12? m: m12; m=m>m13? m: m13;
	m=m>m14? m: m14; m=m>m15? m: m15; m=m>m16? m: m16;
	m=m>m17? m: m17; m=m>m18? m: m18; m=m>m19? m: m19;
	m=m>m20? m: m20; emp(output);
	for(i=1; i<=m; i++)
	{
		char *s1=Line(q),*s2=Line(q),
			 *s3=Line(q),*s4=Line(q),
			 *s5=Line(q),*s6=Line(q),
			 *s7=Line(q),*s8=Line(q),
			 *s9=Line(q),*s10=Line(q),
			 *s11=Line(q),*s12=Line(q),
			 *s13=Line(q),*s14=Line(q),
			 *s15=Line(q),*s16=Line(q),
			 *s17=Line(q),*s18=Line(q),
			 *s19=Line(q),*s20=Line(q);
		if(i<=m1)Liter(x1(i),q,s1); else for(j=0; j<q; j++)s1[j]=' ';
		if(i<=m2)Liter(x2(i),q,s2); else for(j=0; j<q; j++)s2[j]=' ';
		if(i<=m3)Liter(x3(i),q,s3); else for(j=0; j<q; j++)s3[j]=' ';
		if(i<=m4)Liter(x4(i),q,s4); else for(j=0; j<q; j++)s4[j]=' ';
		if(i<=m5)Liter(x5(i),q,s5); else for(j=0; j<q; j++)s5[j]=' ';
		if(i<=m6)Liter(x6(i),q,s6); else for(j=0; j<q; j++)s6[j]=' ';
		if(i<=m7)Liter(x7(i),q,s7); else for(j=0; j<q; j++)s7[j]=' ';
		if(i<=m8)Liter(x8(i),q,s8); else for(j=0; j<q; j++)s8[j]=' ';
		if(i<=m9)Liter(x9(i),q,s9); else for(j=0; j<q; j++)s9[j]=' ';
		if(i<=m10)Liter(x10(i),q,s10); else for(j=0; j<q; j++)s10[j]=' ';
		if(i<=m11)Liter(x11(i),q,s11); else for(j=0; j<q; j++)s11[j]=' ';
		if(i<=m12)Liter(x12(i),q,s12); else for(j=0; j<q; j++)s12[j]=' ';
		if(i<=m13)Liter(x13(i),q,s13); else for(j=0; j<q; j++)s13[j]=' ';
		if(i<=m14)Liter(x14(i),q,s14); else for(j=0; j<q; j++)s14[j]=' ';
		if(i<=m15)Liter(x15(i),q,s15); else for(j=0; j<q; j++)s15[j]=' ';
		if(i<=m16)Liter(x16(i),q,s16); else for(j=0; j<q; j++)s16[j]=' ';
		if(i<=m17)Liter(x17(i),q,s17); else for(j=0; j<q; j++)s17[j]=' ';
		if(i<=m18)Liter(x18(i),q,s18); else for(j=0; j<q; j++)s18[j]=' ';
		if(i<=m19)Liter(x19(i),q,s19); else for(j=0; j<q; j++)s19[j]=' ';
		if(i<=m20)Liter(x19(i),q,s20); else for(j=0; j<q; j++)s20[j]=' ';
		for(j=0; j<q; j++) output<< s1[j]; output<< " ";
		for(j=0; j<q; j++) output<< s2[j]; output<< " ";
		for(j=0; j<q; j++) output<< s3[j]; output<< " ";
		for(j=0; j<q; j++) output<< s4[j]; output<< " ";
		for(j=0; j<q; j++) output<< s5[j]; output<< " ";
		for(j=0; j<q; j++) output<< s6[j]; output<< " ";
		for(j=0; j<q; j++) output<< s7[j]; output<< " ";
		for(j=0; j<q; j++) output<< s8[j]; output<< " ";
		for(j=0; j<q; j++) output<< s9[j]; output<< " ";
		for(j=0; j<q; j++) output<< s10[j]; output<< " ";
		for(j=0; j<q; j++) output<< s11[j]; output<< " ";
		for(j=0; j<q; j++) output<< s12[j]; output<< " ";
		for(j=0; j<q; j++) output<< s13[j]; output<< " ";
		for(j=0; j<q; j++) output<< s14[j]; output<< " ";
		for(j=0; j<q; j++) output<< s15[j]; output<< " ";
		for(j=0; j<q; j++) output<< s16[j]; output<< " ";
		for(j=0; j<q; j++) output<< s17[j]; output<< " ";
		for(j=0; j<q; j++) output<< s18[j]; output<< " ";
		for(j=0; j<q; j++) output<< s19[j]; output<< " ";
		for(j=0; j<q; j++) output<< s20[j];
		if(i<m)output<< "\n";
		delete []s1; delete []s2; delete []s3; delete []s4;
		delete []s5; delete []s6; delete []s7; delete []s8;
		delete []s9; delete []s10; delete []s11; delete []s12;
		delete []s13; delete []s14; delete []s15; delete []s16;
		delete []s17; delete []s18; delete []s19; delete []s20;
	}
}

template <class TT,class T>
void Ecrir(TT output,T x){output<< x;}

template <class TT,class T1,class T2>
void Ecrir(TT output,T1 x1,T2 x2){output<< x1<< " "<< x2;}

template <class TT,class T1,class T2,class T3>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3){output<< x1<< " "<< x2<< " "<< x3;}

template <class TT,class T1,class T2,class T3,class T4>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4)
{output<< x1<< " "<< x2<< " "<< x3<< " "<< x4;}

template <class TT,class T1,class T2,class T3,class T4,class T5>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5)
{output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5;}

template <class TT,class T1,class T2,class T3,class T4,class T5,class T6>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6)
{output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6;}

template <class TT,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "
		  << x4<< " "<< x5<< " "<< x6<< " "<< x7;
}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,
		   T4 x4,T5 x5,T6 x6,T7 x7,T8 x8)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4
	      << " "<< x5<< " "<< x6<< " "<< x7<< " "<< x8;
}

template <class TT,class T1,class T2,class T3,class T4,
		  class T5,class T6,class T7,class T8,class T9>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,
		   T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,T9 x9)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4
	      << " "<< x5<< " "<< x6<< " "<< x7<< " "<< x8<< " "<< x9;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,
		   T5 x5,T6 x6,T7 x7,T8 x8,T9 x9,T10 x10)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5
		  << " "<< x6<< " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		   T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10<< " "<< x11;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,class T12>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15,T16 x16)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15
		  << " "<< x16;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,class T17>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15,
		   T16 x16,T17 x17)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15
		  << " "<< x16<< " "<< x17;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,
		  class T17,class T18>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15,
		   T16 x16,T17 x17,T18 x18)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15
		  << " "<< x16<< " "<< x17<< " "<< x18;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,
		  class T17,class T18,class T19>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15,
		   T16 x16,T17 x17,T18 x18,T19 x19)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15
		  << " "<< x16<< " "<< x17<< " "<< x18<< " "<< x19;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,
		  class T17,class T18,class T19,class T20>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15,
		   T16 x16,T17 x17,T18 x18,T19 x19,T20 x20)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15
		  << " "<< x16<< " "<< x17<< " "<< x18<< " "<< x19<< " "<< x20;
}

template <class TT,class T1,class T2,class T3,class T4,class T5,
		  class T6,class T7,class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,class T16,
		  class T17,class T18,class T19,class T20,class T21>
void Ecrir(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,T7 x7,T8 x8,
		   T9 x9,T10 x10,T11 x11,T12 x12,T13 x13,T14 x14,T15 x15,
		   T16 x16,T17 x17,T18 x18,T19 x19,T20 x20,T21 x21)
{
	output<< x1<< " "<< x2<< " "<< x3<< " "<< x4<< " "<< x5<< " "<< x6
		  << " "<< x7<< " "<< x8<< " "<< x9<< " "<< x10
		  << " "<< x11<< " "<< x12<< " "<< x13<< " "<< x14<< " "<< x15
		  << " "<< x16<< " "<< x17<< " "<< x18<< " "<< x19<< " "<< x20
		  << " "<< x21;
}

template <class TT,class T1>
void Look(TT output,T1 x1,int q)
{
	int j; char *s1=Line(q); Liter(x1,q,s1);
	for(j=0; j<q; j++) output<< s1[j];
	delete []s1;
}

template <class TT,class T1,class T2>
void Look(TT output,T1 x1,T2 x2,int q)
{
	int j; char *s2=Line(q); Liter(x2,q,s2);
	Look(output,x1,q);
	output<< " "; for(j=0; j<q; j++) output<< s2[j];
	delete []s2;
}

template <class TT,class T1,class T2,class T3>
void Look(TT output,T1 x1,T2 x2,T3 x3,int q)
{
	int j; char *s3=Line(q); Liter(x3,q,s3);
	Look(output,x1,x2,q);
	output<< " "; for(j=0; j<q; j++) output<< s3[j];
	delete []s3;
}

template <class TT,class T1,class T2,class T3,class T4>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,int q)
{
	int j; char *s4=Line(q); Liter(x4,q,s4);
	Look(output,x1,x2,x3,q);
	output<< " "; for(j=0; j<q; j++) output<< s4[j];
	delete []s4;
}

template <class TT,class T1,class T2,class T3,class T4,class T5>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,int q)
{
	int j; char *s5=Line(q); Liter(x5,q,s5);
	Look(output,x1,x2,x3,x4,q);
	output<< " "; for(j=0; j<q; j++) output<< s5[j];
	delete []s5;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,T6 x6,int q)
{
	int j; char *s6=Line(q); Liter(x6,q,s6);
	Look(output,x1,x2,x3,x4,x5,q);
	output<< " "; for(j=0; j<q; j++) output<< s6[j];
	delete []s6;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,int q)
{
	int j; char *s7=Line(q); Liter(x7,q,s7);
	Look(output,x1,x2,x3,x4,x5,x6,q);
	output<< " "; for(j=0; j<q; j++) output<< s7[j];
	delete []s7;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,int q)
{
	int j; char *s8=Line(q); Liter(x8,q,s8);
	Look(output,x1,x2,x3,x4,x5,x6,x7,q);
	output<< " "; for(j=0; j<q; j++)output<< s8[j]; delete []s8;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,int q)
{
	int j; char *s9=Line(q); Liter(x9,q,s9);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,q);
	output<< " "; for(j=0; j<q; j++)output<< s9[j]; delete []s9;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,int q)
{
	int j; char *s10=Line(q); Liter(x10,q,s10);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,q);
	output<< " "; for(j=0; j<q; j++)output<< s10[j]; delete []s10;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,int q)
{
	int j; char *s11=Line(q); Liter(x11,q,s11);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,q);
	output<< " "; for(j=0; j<q; j++)output<< s11[j]; delete []s11;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,int q)
{
	int j; char *s12=Line(q); Liter(x12,q,s12);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,q);
	output<< " "; for(j=0; j<q; j++)output<< s12[j]; delete []s12;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,int q)
{
	int j; char *s13=Line(q); Liter(x13,q,s13);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,q);
	output<< " "; for(j=0; j<q; j++)output<< s13[j]; delete []s13;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,int q)
{
	int j; char *s14=Line(q); Liter(x14,q,s14);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,q);
	output<< " "; for(j=0; j<q; j++)output<< s14[j]; delete []s14;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,int q)
{
	int j; char *s15=Line(q); Liter(x15,q,s15);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,q);
	output<< " "; for(j=0; j<q; j++)output<< s15[j]; delete []s15;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,int q)
{
	int j; char *s16=Line(q); Liter(x16,q,s16);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,q);
	output<< " "; for(j=0; j<q; j++)output<< s16[j]; delete []s16;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,int q)
{
	int j; char *s17=Line(q); Liter(x17,q,s17);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,q);
	output<< " "; for(j=0; j<q; j++)output<< s17[j]; delete []s17;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,int q)
{
	int j; char *s18=Line(q); Liter(x18,q,s18);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,q);
	output<< " "; for(j=0; j<q; j++)output<< s18[j]; delete []s18;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18,class T19>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,T19 x19,int q)
{
	int j; char *s19=Line(q); Liter(x19,q,s19);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,
		 x18,q);
	output<< " "; for(j=0; j<q; j++)output<< s19[j]; delete []s19;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18,class T19,class T20>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,T19 x19,T20 x20,int q)
{
	int j; char *s20=Line(q); Liter(x20,q,s20);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,
		 x18,x19,q);
	output<< " "; for(j=0; j<q; j++)output<< s20[j]; delete []s20;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18,class T19,
		  class T20,class T21>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,T19 x19,T20 x20,T21 x21,int q)
{
	int j; char *s21=Line(q); Liter(x21,q,s21);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,
		 x18,x19,x20,q);
	output<< " "; for(j=0; j<q; j++)output<< s21[j]; delete []s21;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,T19 x19,T20 x20,T21 x21,
		  T22 x22,int q)
{
	int j; char *s22=Line(q); Liter(x22,q,s22);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,
		 x18,x19,x20,x21,q);
	output<< " "; for(j=0; j<q; j++)output<< s22[j]; delete []s22;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,T19 x19,T20 x20,T21 x21,
		  T22 x22,T23 x23,int q)
{
	int j; char *s23=Line(q); Liter(x23,q,s23);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,
		 x18,x19,x20,x21,x22,q);
	output<< " "; for(j=0; j<q; j++)output<< s23[j]; delete []s23;
}

template <class TT,class T1,class T2,class T3,
		  class T4,class T5,class T6,class T7,
		  class T8,class T9,class T10,class T11,
		  class T12,class T13,class T14,class T15,
		  class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,
		  class T24>
void Look(TT output,T1 x1,T2 x2,T3 x3,T4 x4,T5 x5,
		  T6 x6,T7 x7,T8 x8,T9 x9,T10 x10,T11 x11,
		  T12 x12,T13 x13,T14 x14,T15 x15,T16 x16,
		  T17 x17,T18 x18,T19 x19,T20 x20,T21 x21,
		  T22 x22,T23 x23,T24 x24,int q)
{
	int j; char *s24=Line(q); Liter(x24,q,s24);
	Look(output,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,
		 x18,x19,x20,x21,x22,x23,q);
	output<< " "; for(j=0; j<q; j++)output<< s24[j]; delete []s24;
}

//Универсальный вывод критических точек.
//T0,Tf - точки входа в критическую область и выхода из нее.
template <class C> void OutPoints(C out,int &N,int Lout,double t0,double ud0,Vect& T0,Vect& Tf)
{
	int i; double t00,tf0,t,tf,udf;
	t00=T0(1)-int(T0(1)/86400)*86400; //Гринв. время первой входной КТ
	tf0=Tf(1)-int(Tf(1)/86400)*86400; //Гринв. время первой выходной КТ
	emp(out); ud0=ud0+int((T0(1)-t0)/86400);
	if(t00<t0)ud0++; udf=ud0;
	if(Tf(1)<T0(1)){N--; T0=Cut(T0,1); Tf=Cut(Tf,1);}
	if(tf0<t00 && T0(1)<Tf(1))udf++;
	Ecrir(out," count=",N,"\n");
	DO(i,1,N)
	{
		if(i>1)
		{
			t=T0(i)-int(T0(i)/86400)*86400;
			tf=Tf(i)-int(Tf(i)/86400)*86400;
			if(t<t00)ud0++; if(tf<tf0)udf++;
			t00=t; tf0=tf;
		}
		if(fmod(i-1,50)==0)
		{
			if(i>1)emp(out);
			Look(out,"N","ulian1","ulian2","date1  ","date2  ","t1","t2",
				 "t2-t1","time1","time2","Length",Lout);
			emp(out);
		}
		if(!(i==1 && Tf(1)<T0(1)))
		{
			Look(out,i,ud0,udf,Char(DMY(ud0)),Char(DMY(udf)),T0(i),Tf(i),
				 Tf(i)-T0(i),Char(HMS(T0(i))),Char(HMS(Tf(i))),
				 Char(HMS(Tf(i)-T0(i))),Lout);
			emp(out);
		}
	}
}

template <class TT,class T1,class T2> void Write(TT &output,T1 x,T2 q)
{
	char *s=Line(q); Liter(x,q,s);
	for(int j=0; j<q; j++)output<< s[j];
	delete []s;
}

template <class TT,class T1,class T2,class T3,class T4>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2)
{
    Write(output,x1,q1);
	char *s2=Line(q2); Liter(x2,q2,s2);
	output<< " "; for(int j=0; j<q2; j++) output<< s2[j];
	delete []s2;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3)
{
    Write(output,x1,q1,x2,q2);
	char *s3=Line(q3); Liter(x3,q3,s3);
	output<< " "; for(int j=0; j<q3; j++) output<< s3[j];
	delete []s3;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,T8 q4)
{
    Write(output,x1,q1,x2,q2,x3,q3);
	char *s4=Line(q4); Liter(x4,q4,s4);
	output<< " "; for(int j=0; j<q4; j++) output<< s4[j];
	delete []s4;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4);
	char *s5=Line(q5); Liter(x5,q5,s5);
	output<< " "; for(int j=0; j<q5; j++) output<< s5[j];
	delete []s5;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5);
	char *s6=Line(q6); Liter(x6,q6,s6);
	output<< " "; for(int j=0; j<q6; j++) output<< s6[j];
	delete []s6;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6);
	char *s7=Line(q7); Liter(x7,q7,s7);
	output<< " "; for(int j=0; j<q7; j++) output<< s7[j];
	delete []s7;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7);
	char *s8=Line(q8); Liter(x8,q8,s8);
	output<< " "; for(int j=0; j<q8; j++) output<< s8[j];
	delete []s8;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8);
	char *s9=Line(q9); Liter(x9,q9,s9);
	output<< " "; for(int j=0; j<q9; j++) output<< s9[j];
	delete []s9;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,
		  x7,q7,x8,q8,x9,q9);
	char *s10=Line(q10); Liter(x10,q10,s10);
	output<< " "; for(int j=0; j<q10; j++) output<< s10[j];
	delete []s10;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,
		  x7,q7,x8,q8,x9,q9,x10,q10);
	char *s11=Line(q11); Liter(x11,q11,s11);
	output<< " "; for(int j=0; j<q11; j++) output<< s11[j];
	delete []s11;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,
		  x7,q7,x8,q8,x9,q9,x10,q10,x11,q11);
	char *s12=Line(q12); Liter(x12,q12,s12);
	output<< " "; for(int j=0; j<q12; j++) output<< s12[j];
	delete []s12;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,
		  x7,q7,x8,q8,x9,q9,x10,q10,x11,q11,x12,q12);
	char *s13=Line(q13); Liter(x13,q13,s13);
	output<< " "; for(int j=0; j<q13; j++) output<< s13[j];
	delete []s13;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,
		  x7,q7,x8,q8,x9,q9,x10,q10,x11,q11,x12,q12,x13,q13);
	char *s14=Line(q14); Liter(x14,q14,s14);
	output<< " "; for(int j=0; j<q14; j++) output<< s14[j];
	delete []s14;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,
		  x8,q8,x9,q9,x10,q10,x11,q11,x12,q12,x13,q13,x14,q14);
	char *s15=Line(q15); Liter(x15,q15,s15);
	output<< " "; for(int j=0; j<q15; j++) output<< s15[j];
	delete []s15;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,
		  x9,q9,x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15);
	char *s16=Line(q16); Liter(x16,q16,s16);
	output<< " "; for(int j=0; j<q16; j++) output<< s16[j];
	delete []s16;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32,
		  class T33,class T34>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16,
		   T33 x17,T34 q17)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,
		  x9,q9,x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15,x16,q16);
	char *s17=Line(q17); Liter(x17,q17,s17);
	output<< " "; for(int j=0; j<q17; j++) output<< s17[j];
	delete []s17;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32,
		  class T33,class T34,class T35,class T36>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16,
		   T33 x17,T34 q17,T35 x18,T36 q18)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,x9,q9,
		  x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15,x16,q16,x17,q17);
	char *s18=Line(q18); Liter(x18,q18,s18);
	output<< " "; for(int j=0; j<q18; j++) output<< s18[j];
	delete []s18;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32,
		  class T33,class T34,class T35,class T36,
		  class T37,class T38>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16,
		   T33 x17,T34 q17,T35 x18,T36 q18,T37 x19,T38 q19)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,x9,q9,
		  x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15,x16,q16,x17,q17,
		  x18,q18);
	char *s19=Line(q19); Liter(x19,q19,s19);
	output<< " "; for(int j=0; j<q19; j++) output<< s19[j];
	delete []s19;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32,
		  class T33,class T34,class T35,class T36,
		  class T37,class T38,class T39,class T40>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16,
		   T33 x17,T34 q17,T35 x18,T36 q18,T37 x19,T38 q19,
		   T39 x20,T40 q20)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,x9,q9,
		  x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15,x16,q16,x17,q17,
		  x18,q18,x19,q19);
	char *s20=Line(q20); Liter(x20,q20,s20);
	output<< " "; for(int j=0; j<q20; j++) output<< s20[j];
	delete []s20;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32,
		  class T33,class T34,class T35,class T36,
		  class T37,class T38,class T39,class T40,
		  class T41,class T42>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16,
		   T33 x17,T34 q17,T35 x18,T36 q18,T37 x19,T38 q19,
		   T39 x20,T40 q20,T41 x21,T42 q21)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,x9,q9,
		  x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15,x16,q16,x17,q17,
		  x18,q18,x19,q19,x20,q20);
	char *s21=Line(q21); Liter(x21,q21,s21);
	output<< " "; for(int j=0; j<q21; j++) output<< s21[j];
	delete []s21;
}

template <class TT,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9,
          class T10,class T11,class T12,class T13,class T14,
		  class T15,class T16,class T17,class T18,class T19,
		  class T20,class T21,class T22,class T23,class T24,
		  class T25,class T26,class T27,class T28,
		  class T29,class T30,class T31,class T32,
		  class T33,class T34,class T35,class T36,
		  class T37,class T38,class T39,class T40,
		  class T41,class T42,class T43,class T44>
void Write(TT output,T1 x1,T2 q1,T3 x2,T4 q2,T5 x3,T6 q3,T7 x4,
	       T8 q4,T9 x5,T10 q5,T11 x6,T12 q6,T13 x7,T14 q7,
		   T15 x8,T16 q8,T17 x9,T18 q9,T19 x10,T20 q10,
		   T21 x11,T22 q11,T23 x12,T24 q12,T25 x13,T26 q13,
		   T27 x14,T28 q14,T29 x15,T30 q15,T31 x16,T32 q16,
		   T33 x17,T34 q17,T35 x18,T36 q18,T37 x19,T38 q19,
		   T39 x20,T40 q20,T41 x21,T42 q21,T43 x22,T44 q22)
{
    Write(output,x1,q1,x2,q2,x3,q3,x4,q4,x5,q5,x6,q6,x7,q7,x8,q8,x9,q9,
		  x10,q10,x11,q11,x12,q12,x13,q13,x14,q14,x15,q15,x16,q16,x17,q17,
		  x18,q18,x19,q19,x20,q20,x21,q21);
	char *s22=Line(q22); Liter(x22,q22,s22);
	output<< " "; for(int j=0; j<q22; j++) output<< s22[j];
	delete []s22;
}

template <class T> void StrOut(T &output,Str &x,int q)
{Str s(q); s.s=Liter(x.s,q); output<< s.s;}

template <class T>
void StrOut(T &output,Str &x1,int q1,Str &x2,int q2)
{
	Str s1(q1),s2(q2);
	s1.s=Liter(x1.s,q1); s2=Liter(x2.s,q2);
	output<< s1.s<< " "<< s2.s;
}

template <class T>
void StrOut(T &output,Str &x1,int q1,
			  Str &x2,int q2,Str &x3,int q3)
{
	Str s1(q1),s2(q2),s3(q3);
	s1.s=Liter(x1.s,q1); s2=Liter(x2.s,q2); s3=Liter(x3.s,q3);
	output<< s1.s<< " "<< s2.s<< " "<< s3.s;
}

template <class T>
void StrOut(T &output,Str &x1,int q1,Str &x2,int q2,
			  Str &x3,int q3,Str &x4,int q4)
{
	Str s1(q1),s2(q2),s3(q3),s4(q4);
	s1.s=Liter(x1.s,q1); s2=Liter(x2.s,q2);
	s3=Liter(x3.s,q3); s4=Liter(x4.s,q4);
	output<< s1.s<< " "<< s2.s<< " "<< s3.s<< " "<< s4.s;
}

template <class T>
void StrOut(T &output,Str &x1,int q1,Str &x2,int q2,
			  Str &x3,int q3,Str &x4,int q4,Str &x5,int q5)
{
	Str s1(q1),s2(q2),s3(q3),s4(q4),s5(q5);
	s1.s=Liter(x1.s,q1); s2=Liter(x2.s,q2);
	s3=Liter(x3.s,q3); s4=Liter(x4.s,q4); s5=Liter(x5.s,q5);
	output<< s1.s<< " "<< s2.s<< " "<< s3.s<< " "<< s4.s<< " "<< s5.s;
}

template <class T>
void StrOut(T &output,Str &x1,int q1,Str &x2,int q2,
			  Str &x3,int q3,Str &x4,int q4,Str &x5,int q5,
			  Str &x6,int q6)
{
	Str s1(q1),s2(q2),s3(q3),s4(q4),s5(q5),s6(q6);
	s1.s=Liter(x1.s,q1); s2=Liter(x2.s,q2);
	s3=Liter(x3.s,q3); s4=Liter(x4.s,q4);
	s5=Liter(x5.s,q5); s6=Liter(x6.s,q6);
	output<< s1.s<< " "<< s2.s<< " "<< s3.s
		  << " "<< s4.s<< " "<< s5.s<< " "<< s6.s;
}

template <class T>
void StrOut(T &output,Str &x1,int q1,Str &x2,int q2,
			Str &x3,int q3,Str &x4,int q4,Str &x5,int q5,
			Str &x6,int q6,Str &x7,int q7)
{
	Str s1(q1),s2(q2),s3(q3),s4(q4),s5(q5),s6(q6),s7(q6);
	s1.s=Liter(x1.s,q1); s2=Liter(x2.s,q2);
	s3=Liter(x3.s,q3); s4=Liter(x4.s,q4);
	s5=Liter(x5.s,q5); s6=Liter(x6.s,q6); s7=Liter(x7.s,q7);
	output<< s1.s<< " "<< s2.s<< " "<< s3.s
		  << " "<< s4.s<< " "<< s5.s
		  << " "<< s6.s<< " "<< s7.s;
}

//Чтение данных из потока entry и запись их в двумерный массив
template <class TT> void ReadMatr(TT entry,Matr& A)
{int i,j,m=A.m,n=A.n; DO(i,1,m)DO(j,1,n)entry>> A(i,j);}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT> void ReadMatr(TT entry,Matr& A1,Matr& A2)
{
	int i,j,m=A1.m,n=A1.n; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT> void ReadMatr(TT entry,Matr& A1,Matr& A2,Matr& A3)
{
	int i,j,m=A1.m,n=A1.n; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
	m=A3.m; n=A3.n; DO(i,1,m)DO(j,1,n)entry>> A3(i,j);
}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT> void ReadMatr(TT entry,Matr& A1,Matr& A2,Matr& A3,Matr& A4)
{
	int i,j,m=A1.m,n=A1.n; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
	m=A3.m; n=A3.n; DO(i,1,m)DO(j,1,n)entry>> A3(i,j);
	m=A4.m; n=A4.n; DO(i,1,m)DO(j,1,n)entry>> A4(i,j);
}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT> 
void ReadMatr(TT entry,Matr& A1,Matr& A2,Matr& A3,Matr& A4,Matr& A5)
{
	int i,j,m=A1.m,n=A1.n; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
	m=A3.m; n=A3.n; DO(i,1,m)DO(j,1,n)entry>> A3(i,j);
	m=A4.m; n=A4.n; DO(i,1,m)DO(j,1,n)entry>> A4(i,j);
	m=A5.m; n=A5.n; DO(i,1,m)DO(j,1,n)entry>> A5(i,j);
}

//Чтение данных из потока entry и запись их в двумерный массив
template <class TT,class S>void ReadMatr(TT entry,S s,Matr& A)
{int i,j,m=A.m,n=A.n; entry>> s; DO(i,1,m)DO(j,1,n)entry>> A(i,j);}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT,class S>
void ReadMatr(TT entry,S s1,Matr& A1,S s2,Matr& A2)
{
	int i,j,m=A1.m,n=A1.n; entry>> s1; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; entry>> s2; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT,class S>
void ReadMatr(TT entry,S s1,Matr& A1,S s2,Matr& A2,S s3,Matr& A3)
{
	int i,j,m=A1.m,n=A1.n; entry>> s1; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; entry>> s2; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
	m=A3.m; n=A3.n; entry>> s3; DO(i,1,m)DO(j,1,n)entry>> A3(i,j);
}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT,class S>
void ReadMatr(TT entry,S s1,Matr& A1,S s2,Matr& A2,S s3,Matr& A3,
			  S s4,Matr& A4)
{
	int i,j,m=A1.m,n=A1.n; entry>> s1; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; entry>> s2; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
	m=A3.m; n=A3.n; entry>> s3; DO(i,1,m)DO(j,1,n)entry>> A3(i,j);
	m=A4.m; n=A4.n; entry>> s4; DO(i,1,m)DO(j,1,n)entry>> A4(i,j);
}

//Чтение данных из потока entry и запись их в двумерные массивы
template <class TT,class S>
void ReadMatr(TT entry,S s1,Matr& A1,S s2,Matr& A2,S s3,Matr& A3,
			  S s4,Matr& A4,S s5,Matr& A5)
{
	int i,j,m=A1.m,n=A1.n; entry>> s1; DO(i,1,m)DO(j,1,n)entry>> A1(i,j);
	m=A2.m; n=A2.n; entry>> s2; DO(i,1,m)DO(j,1,n)entry>> A2(i,j);
	m=A3.m; n=A3.n; entry>> s3; DO(i,1,m)DO(j,1,n)entry>> A3(i,j);
	m=A4.m; n=A4.n; entry>> s4; DO(i,1,m)DO(j,1,n)entry>> A4(i,j);
	m=A5.m; n=A5.n; entry>> s5; DO(i,1,m)DO(j,1,n)entry>> A5(i,j);
}

//Чтение данных из потока entry и запись их в одномерный массив
template <class TT>void ReadVect(TT entry,Vect& a)
{int i,m=a.m; DO(i,1,m)entry>> a(i);}

//Чтение данных из потока entry и запись их в одномерный массив
template <class TT,class S>void ReadVect(TT entry,S s,Vect& a)
{int i,m=a.m; entry>> s; DO(i,1,m)entry>> a(i);}

//Чтение данных из потока entry и запись их в одномерные массивы
template <class TT,class S>
void ReadVect(TT entry,S s1,Vect& a1,S s2,Vect& a2)
{
	int m=a1.m; entry>> s1; DO(i,1,m)entry>> a1(i);
	m=a2.m; entry>> s2; DO(i,1,m)entry>> a2(i);
}

//Чтение данных из потока entry и запись их в одномерные массивы
template <class TT,class S>
void ReadVect(TT entry,S s1,Vect& a1,S s2,Vect& a2,S s3,Vect& a3)
{
	int i,m=a1.m; entry>> s1; DO(i,1,m)entry>> a1(i);
	m=a2.m; entry>> s2; DO(i,1,m)entry>> a2(i);
	m=a3.m; entry>> s3; DO(i,1,m)entry>> a3(i);
}

//Чтение данных из потока entry и запись их в одномерные массивы
template <class TT,class S>
void ReadVect(TT entry,S s1,Vect& a1,S s2,Vect& a2,S s3,Vect& a3,
			  S s4,Vect& a4)
{
	int i,m=a1.m; entry>> s1; DO(i,1,m)entry>> a1(i);
	m=a2.m; entry>> s2; DO(i,1,m)entry>> a2(i);
	m=a3.m; entry>> s3; DO(i,1,m)entry>> a3(i);
	m=a4.m; entry>> s4; DO(i,1,m)entry>> a4(i);
}

//Чтение данных из потока entry и запись их в одномерные массивы
template <class TT,class S>
void ReadVect(TT entry,S s1,Vect& a1,S s2,Vect& a2,S s3,Vect& a3,
			  S s4,Vect& a4,S s5,Vect& a5)
{
	int i,m=a1.m; entry>> s1; DO(i,1,m)entry>> a1(i);
	m=a2.m; entry>> s2; DO(i,1,m)entry>> a2(i);
	m=a3.m; entry>> s3; DO(i,1,m)entry>> a3(i);
	m=a4.m; entry>> s4; DO(i,1,m)entry>> a4(i);
	m=a5.m; entry>> s5; DO(i,1,m)entry>> a5(i);
}

void OutVect(char S[30],Vect& x);
void OutMatr(char S[30],Matr& A);
void MatrixOut(char* ch,Matr& A,int O);

void OutVect(fstream& out,char S[30],Vect &x);
void OutMatr(fstream& out,char S[30],Matr &A);

template <class T> void emp(T out) {out<< "\n";}

template <class T>void emp(T out,int n){for(int i=1; i<=n; i++) out<< "\n";}

#define EC_ Ecrir(cout,

#define LK_ Look(cout,

#define VC_ VecOut(cout,

#define MT_ MatrOut(cout,

#define O_ emp();

Vect AbsOsc(Vect& x);
double Altitude(Vect& x);
double AngleSunOrb(double ud,double t,double i,double Om);
double AngleSunOrb(double ud,double t,Vect& x);
double AngleSunOrb(double ud,Vect& x,double t);
double AnodePrg(double T,double e,double om);
double Anom(double mu,double eps,double tau,double e,double a,double t);
double ArgWid(double eps,Vect& xos,double t);
void Azim(double ud,double t,Vect& x,double& sinA,double& cosA);

//Формирование номинальных орбит 72 навигационных КА
Matr Constellation(bool typeConstel);

//Формирование возмущенных орбит 72 навигационных КА
/*Matr Constellation(double& rand,bool ph,Vect& Sdr,double t);

//Формирование возмущенных орбит 72 навигационных КА
Matr Constellation(double& rand,bool ph,double t,Vect& Sdr);

//Формирование возмущенных орбит 72 навигационных КА
Matr Constellation(bool ph,Matr& dr,double t);

//Определение созвездия навигационных КА, видимых из точки q
//с упорядочиванием реперов по углам возвышения и с их выводом
Vect Constellation(double ud,double t,double FiMin,
				   Vect& q,Matr& xos,int &N,Vect& Fi);

//Определение созвездия навигационных КА, видимых из точки q
//с упорядочиванием реперов по углам возвышения
Vect Constellation(double ud,double t,double FiMin,
				   Vect& q,Matr& xos,int &N);

//Определение созвездия навигационных КА, видимых из точки q
//без упорядочивания реперов по углам возвышения и с их выводом
Vect Constellation(double ud,double t,double FiMin,
				   Matr& xos,Vect& q,int &N,Vect& Fi);

//Определение созвездия навигационных КА, видимых из точки q
//без упорядочивания реперов по углам возвышения
Vect Constellation(double ud,double t,double FiMin,
				   Matr& xos,Vect& q,int &N);

//Определение координат M навигационных КА
//на заданный момент t в ГСК
Matr Constellation(double ud,double t,Matr& xos);

//Определение координат M навигационных КА
//на заданный момент t в J2000
Matr Constellation(double t,Matr& xos);

//Определение ВС M навигационных КА
//на заданный момент t в ГСК
Matr Constellation(Matr& xos,double ud,double t);

//Определение ВС M навигационных КА
//на заданный момент t в J2000
Matr Constellation(Matr& xos,double t);*/

void Date(double ud,int &d,int &m,int &y);
Vect Date(double ud);
double Derivative(double u,Vect& xos);
double dGamma(double fip,double H);
void dGamma(double B,double H,double& dG);
Str DHMS(double ud,double t);
double dHorizon(double ud,double t,Vect& q,Vect& x);
double dHorizon(Vect& q,Vect& x);
Vect DMY(char *dat);
Str DMY(double ud);
double Eclipte(const double ud);
Vect ElemOrb(double mu,double t,Vect& x);
double FuncA(int d,int m,int y);
Vect GeodGrinv(Vect& g);
Vect GrinvGeod(Vect& x);
void GravJ2000(int n,double ud,double t,Vect& c,Vect& d,Vect& x,Vect& dg);
void GravGrinv(int n,Vect& c,Vect& d,Vect& x,Vect& dg);
Vect GravSphere(int n,Vect& c,Vect& d,Vect& s);
Matr GrinvTop(double Fi,double Lam);
Matr GrinvTop(Vect& q);
Vect HMS(char *tim);
Str HMS(double t);
double Horizon(double ud,double t,Vect& q,Vect& x);
double Horizon(Vect& q,Vect& x);
void Horizon(Vect& q,Vect& x,double& A,double& gamma);
void ISA_IGI(int d,int m,int y,double& F107,double& F81,double& Kp);
Vect Kepler(double mu,double eps,double t,Vect& x);
Vect Kepler(double mu,double eps,double t,Vect& x,double& u);
Vect Kepler(double t,Vect& xosC);
double LambdaGeogr(double s0,double t,double T,double u,Vect& xos);
double LambdaGeogr(double ud,double t,Vect& xos);
Matr MatNut(double ud);
Matr MatPrec(double ud);
Matr MatRot(int k,double u);
Matr Mat2000Gr(double ud,double t);
Matr Mat2000GrG(double ud,double t);
Matr MatGr2000G(double ud,double t);
Matr Mat2000GrV(double ud,double t);
Matr MatGr2000V(double ud,double t);
Matr Mat2000Orb(Vect& xos);
Matr Mat2000OrbG(Vect& xos);
Matr MatAbsOrb(Vect& x);
Matr MatAbsOrbG(Vect& x);
Matr MatOrbAbsG(Vect& x);
Vect ModelPrediction(double t,Vect& x);
void Moon(double JD,Vect& r);
Vect Moon(double udt);
void Nutn(const double ud,double& hpci,double& hepc);
Vect OscAbs(Vect& xos);
Vect Periodicity(Vect& Date,int Mmax,double alfaMax,double FiObject,
				 double LambdaObject,Vect& xos,Vect& dL,Matr& Nv,
				 Matr& Tv,Vect& NvTotal,Vect& TvTotal,double Ptotal);
double Pers(Vect& x);
double Pers(double a);
void PolinomM(double x,Matr& a,Vect& f);
void PolinomS(double x,Matr& a,Vect& f);
double RateFi(double u,Vect& xos);
Matr Reduct(double ud,double t);
Matr Reduct(double ud,double t,double& s);
Matr ReductG(double ud,double t);
int SelectF0(double F81);
int Silvestr(Matr& P);
int Repositive(double eps, double gamma, Matr& P);
void Decomp(fstream& out, double eps, Matr& P, Matr& S, Vect& Lam);
Vect Soleil0(double ud,double t);
Vect Soleil(double ud,double tg);
double Snul(double ud);
void SpeedMoon(double ud,double t,Vect& x,Vect& a);
void SpeedSun(double ud,double t,Vect& x,Vect& a);
void SpeedAtm(int typeAtm,double ud,double t,double b,Vect& x,Vect& dg);
double StAtmos(double H);
double Stime(double udt);
double StimeVer(double ud,double t);
void Sun(double ud,Vect& r);
Vect Sun(double udt);
Vect SunV(double ud,double t,double dt);
Vect SunV(double ud,double t);
Vect MoonV(double ud,double t,double dt);
Vect MoonV(double ud,double t);
double SunSynOrbit(double e,double a);
double TimePrg(double mu,double tau,double e,double a,double teta);
double TimeGrinv(double ud,double t,Vect& x);
double TimeLocal(double ud,double t,Vect& x);
double TraceOffset(double e,double a,double i);
void TraceOffset(double e,double a,double i,double& dL);
double Uday(int d,int m,int y);
double Ulian(int d,int m,int y);
Vect VectKinMom(double i,double Om);
void TimeRadio(double F(Vect& G,double e),Vect& xos,Vect& q,double eps,
			   double GamMin,double Lrn,double& T0,double& Tf);
void TimeRadio(double F(Vect& G,double e),double eps,double GamMin,
			   double ud,double t0,double tf,Vect& xos,Vect& q,
			   int &count,Vect& tauP,Vect& T0,Vect& Tf);
void Visible(double fiMax,double Hion,double t,Matr& XOS,Vect& x,
			 int &n1,int &n2,Vect& N1,Vect& N2);
void Visible(double fiMax,double Hion,double t,
			 Matr& XOS,Vect& x,int &n,Vect& N);
void Visible(double fiMax,double t,Matr& XOS,Vect& x,int &n,Vect& N);
//Модификация для TimeRadio. Определение нуля функции на отрезке. Дихотомия
bool BiSec(double Func(Vect& G,double t),Vect& G,double t0,
		   double tf,double epsf,double epst,double& t);

