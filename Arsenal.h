#pragma once
#include <fstream>
#include <time.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdarg.h>


#include <signal.h>

#include <float.h>
#include <complex>


using namespace std;

extern double pi, rdn, pi_4, pi_2, pi2;
extern double rE_, aE_, bE_, acm_, mu_, muM_, muS_, omE_, Apr_, vC_, ro0_;

class Vect;	//Предварительное объявление класса для
			//использования объектов типа Vect
			//в определении класса Matr

class Matr
{
public:
	static int numbMatr, numbMatrMx, NumbMatr, NumbVect;
	int m, n;
	double** A;
	Matr* H;
	Vect* h;
	Matr() {}	        //M1. Конструктор без параметров
	Matr(int p, int q);	//M2. Прототип "основного" конструктора
	~Matr();	        //M3. Прототип деструктора
	Matr* operator->() { return this; } //M4. Доступ к компонентам класса 
	double& operator()(int i, int j);    //M5. Определение значения (i, j)-го элемента
	Matr& operator()(int i1, int i2, int j1, int j2); //M6. Сечение матрицы
	Vect& operator()(int i1, int i2, int i);		  //M7. Сечение i-й строки матрицы
	Matr& operator()(int k, int n0, Vect& a); //M8. Копирование вектора в сечение k-й строки матрицы
	Vect& operator()(int Ns);				  //M9. Выделение из матрицы строки под номером Ns
	Vect& operator[](int Nc);				  //M10. Выделение из матрицы столбца под номером Nc
	Matr& operator()(Vect& a, int m0, int k); //M11. Копирование вектора в сечение k-го столбца матрицы 
	Matr& operator()(Vect& s); //M12. Выделение строк матрицы, на которые указывает вектор s
	Matr& operator[](Vect& c); //M13. Выделение столбцов матрицы, на которые указывает вектор с
	Matr& operator()(Vect& s, Vect& c); //M14. Выделение подматрицы из матрицы
	Matr& operator-(Matr& B);			//M15. Разность матриц
	Matr& operator+(Matr& B);			//M16. Сумма матриц
	Matr& operator*(Matr& B);			//M17. Перемножение матриц
	Vect& operator*(Vect& b);			//M18. Перемножение матрицы на вектор
	Matr& operator-(double q);			//M19. Вычитание скаляра из матрицы
	Matr& operator+(double q);			//M20. Сумма матрицы со скаляром
	Matr& operator*(double q);			//M21. Умножение матрицы на скаляр
	Matr& operator/(double q);			//M22. Деление матрицы на скаляр
	Matr& operator/(Matr& B);			//M23. Деление матрицы на матрицу
	Matr& operator!();					//M24. Транспонирование матрицы
	Matr& operator-();					//M25. Унарный минус
	Matr& operator~();					//M26. Обращение матрицы
	Matr& operator=(double q);			//M27. Присваивание скаляра матрицей
	Matr& operator=(Matr& B);			//M28. Присваивание матрицы матрицей
	Matr& operator()(int k, Vect& a);   //M29. Копирование вектора в k-ю строку матрицы
	Matr& operator()(Vect& a, int k);   //M30. Копирование вектора в k-й столбец матрицы
	Matr& operator()(Matr& B, int m0, int n0); //M31. Копирование матрицы в сечение матрицы
	bool TestZero();			//M32. Тестирование на ноль. Возвращает true, если матрица нулевая
	Matr& RootMat();			//M33. Извлечение корня квадратного из симметричной матрицы P>0
	Matr& Pow(double q);	    //M34. Поэлементное возведение матрицы в степень
	Matr& Dmat();				//M35. Обнуление недиагональных элементов матрицы
	Matr& Dmat(Vect& a);		//M36. Формирование диагональной матрицы по вектору
	Matr& Dmat(double q);		//M37. Формирование диагональной матрицы по скаляру
	Matr& Grevill();			//M38. Определение псевдообратной матрицы
	Matr& TrStr(int i, int j);  //M39. Перестановка двух столбцов матрицы
	Matr& TrCol(int i, int j);  //M40. Перестановка двух столбцов матрицы
	Matr& SubMatr();			//M41. Перестановка базисных строк
	Matr& SubMatr(int rang, Vect& Am); //M42. Перестановка базисных строк
	bool Negativ();        //M43. Возвращает true, если обнаруживает отрицательный диагональный элемент
	Vect& MatVecStr();				//M44. Упаковка матрицы в вектор по строкам
	Vect& MatVecCol();				//M45. Упаковка матрицы в вектор по столбцам
	Matr& VecMatStr(Vect& a);		//M46. Упаковка вектора в матрицу по строкам
	Matr& VecMatCol(Vect& a);		//M47. Упаковка вектора в матрицу по столбцам
	int Silvestr();					//M48. Возвращает true, если матрица положительно определенная
	Matr& Cut(int I, int J);		//M49. Удаляет из матрицы I-ю строку и J-й столбец
	double Max();					//M50. Определение значения максимального элемента матрицы
	double Min();					//M51. Определение значения минимального элемента матрицы
	void M(int p, int q, Matr& X);	//M52. Функция, заменяющая основной конструктор класса Matr
	Matr& M(int p, int q);			//M53. Функция, заменяющая основной конструктор класса Matr
	Matr& MatrixI(int n);		    //M54. Единичная матрица
};

//Прототипы глобальных операторов для класса Matr
Matr& operator-(double q, Matr& A); //MG1. Вычитание матрицы A из скаляра q
Matr& operator+(double q, Matr& A); //MG2. Сложение скаляра q с матрицей A
Matr& operator*(double q, Matr& A); //MG3. Умножение скаляра q с матрицей A
Matr& operator/(double q, Matr& A); //MG4. Деление скаляра q на матрицу A
Matr& operator&(Matr& A, Matr& B);  //MG5. Прямое произведение матриц

class Vect
{
public:
	static int numbVect, numbVectMx, NumbVect, NumbMatr;
	int m; //Размерность вектора
	double* a;
	Vect* h;
	Matr* H;
	Vect() {}	  //V1. Конструктор без параметров			 
	Vect(int p);  //V2. Прототип "основного" конструктора
	~Vect();      //V3. Прототип деструктора
	double& operator()(int i); //V4. Определение значения i-й координаты вектора
	double operator*(Vect& x); //V5. Скалярное произведение векторов
	double operator!();        //V6. Модуль вектора
	Vect& operator+(Vect& x);  //V7. Сложение двух векторов
	Vect& operator-(Vect& x);  //V8. Вычитание двух векторов
	Vect& operator=(double q); //V9. Присваивание скаляра вектору
	Vect& operator=(Vect& x);  //V10. Присваивание вектора вектору
	Vect& operator*(double q); //V11. Умножение вектора на скаляр
	Vect& operator/(double q); //V12. Деление вектора на скаляр
	Vect* operator->() { return this; } //V13. Доступ к компонентам класса
	Vect& operator()(int i1, int i2);   //V14. Сечение вектора
	Vect& operator()(Vect& pointer);    //V15. Выделение компонент, на которые указывает pointer
	Vect& operator*(Matr& A);   //V16. Умножение вектор-строки на матрицу
	Matr& operator^(Vect& x);   //V17. Диадное произведение
	double operator%(Vect& b);  //V18. Угол между векторами
	Vect& operator-(double q);  //V19. Вычитание скаляра из вектора
	Vect& operator+(double q);  //V20. Сложение скаляра с вектором
	Vect& operator/(Vect& b);   //V21. Деление вектора на вектор
	Vect& operator()(Vect& x, int m0);  //V22. Копирование вектора в сечение вектора
	Vect& operator()(int m0, Vect& x);  //V23. Копирование вектора в сечение вектора
	Vect& operator()(int m0, double x); //V24. Инициализация сечения вектора (от m0 до m) числом
	Vect& operator()(int m1, int m2, double x); //V25. Инициализация сечения вектора числом
	Vect& operator()(double x, int m1, int m2); //V26. Инициализация сечения вектора числом
	Vect& operator-();            //V27. Унарный минус
	Vect& Pow(double q);          //V28. Возведение в степень
	Vect& Dvec(Matr& A);          //V29. Главная диагональ квадратной матрицы
	bool TestZero();              //V30. Тестирование на ноль
	Vect& TrStr(int i, int j);    //V31. Перестановка местами i-ой и j-ой компонент
	double Min();			      //V32. Определение значения минимальной компоненты вектора
	double Max();			      //V33. Определение значения максимальной компоненты вектора
	void V(int p, Vect& x);	      //V34. Замена основного конструктора
	Vect& V(int p);			      //V35. Замена основного конструктора
	Vect& Vector_e(int n, int k); //V36. Единичный вектор
};

//Прототипы глобальных операторов для класса Vect
Vect& operator&(Vect& x, Vect& y);  //VG1. Прямое произведение векторов
Vect& operator*(double q, Vect& x); //VG2. Умножение скаляра на вектор
Vect& operator-(double q, Vect& x); //VG3. Вычитание из скаляра вектора
Vect& operator+(double q, Vect& x); //VG4. Сложение скаляра с вектором
Vect& operator/(double q, Vect& x); //VG5. Деление скаляра на вектор

class Tensor
{
public:
	static int numbTensor, numbTensorMx;
	int m1, m2, m3; double*** T;
	Tensor() {};						//Êîíñòðóêòîð áåç ïàðàìåòðîâ
	Tensor(int n1, int n2, int n3);	//Ïðîòîòèï "îñíîâíîãî" êîíñòðóêòîðà
	~Tensor();			            //Ïðîòîòèï äåñòðóêòîðà
	double& operator()(int k1, int k2, int k3);
	Vect& operator()(char* c, int k2, int k3);
	Vect& operator()(int k1, char* c, int k3);
	Vect& operator()(int k1, int k2, char* c);
	Tensor& operator()(Vect& x, int k2, int k3);
	Tensor& operator()(int k1, Vect& x, int k3);
	Tensor& operator()(int k1, int k2, Vect& x);
	Matr& operator()(int k1, char* c2, char* c3);
	Matr& operator()(char* c1, int k2, char* c3);
	Matr& operator()(char* c1, char* c2, int k3);
	Tensor& operator()(int k1, char* c2, char* c3, Matr& A);
	Tensor& operator()(char* c1, int k2, char* c3, Matr& A);
	Tensor& operator()(char* c1, char* c2, int k3, Matr& A);
	Tensor* operator->() { return this; }
	void del();
};

class Matr;

class Str
{
public:
	static int numbStr, numbStrMx, NumbStr;
	int m;
	char* s;
	Str* h;
	Str() {};	//S1. Конструктор без параметров
	Str(int p);	//S2. Прототип "основного" конструктора
	~Str();	    //S3. Прототип деструктора
	char& operator()(int i);		    //S4. Доступ к i-му символу
	Str& operator()(int i1, int i2);	//S5. Сечение строки
	Str& operator+(Str& x);			    //S6. Сцепление (конкатенация) строк
	Str operator+(char* x);				//S7. Сцепление (конкатенация) строки с переменной типа char
	Str& operator=(Str& x);				//S8. Оператор присвоения строковой переменной
	Str& operator=(char* x);			//S9. Оператор присвоения переменной типа char
	Str* operator->() { return this; }  //S10. Доступ к компонентам класса
};

//Глобальные операторы для класса Str
Str& operator+(char* x, Str& y); //SG1. Сцепление (конкатенация) строки с переменной типа char

class Ballist
{
public:
	bool Start, startMoon, startSun, startAtm;
	int Ngrav, //Порядок разложения ГПЗ в ряд по сферическим функциям
		OnMoon, OnSun, OnLux, typeAtm, //Ключи учета возмущений
		kday, F0Atm, indexAtm, nday, day, month, year;
	double ud,		//Юлианская дата
		   BalCoef, //Баллистический коэффициент
		   A_dAtm, F107Atm, F81Atm, KpAtm, ud0Atm,
		   gamma, //Коэффициент  отражения (поглащения)
		   p, //Давление солнечной радиации вблизи Земли
		   Sm, //Площадь поперечного сечения КА
		   mass; //Масса КА
	Matr* Aatm, * Batm, * Catm, * Datm, * Eatm, * Fi1Atm, * Latm, * Natm;
	Ballist(); //B1. Конструктор класса Ballist
	Ballist(int N); //B2. Конструктор класса Ballist
	~Ballist(); //B3. Деструктор класса Ballist
	void Init(); //B4. Инициализация компонент класса Ballist
	Vect& ModProgn(double t, Vect& x);         //B5. Вычисление правых частей дифференциальных уравнений
	Vect& ModProgn(Vect& U(double t, Vect& x), //B6. Вычисление правых частей дифференциальных уравнений 
		           double t, Vect& x);
	Vect& sModProgn(double t, Vect& x); //B7. Модель движения, учитывающая только коэфф. c20
	void GravJ2000(double t, Vect& c,  //B8. Определение компонент гравитационного ускорения в J2000
		           Vect& d, Vect& x, Vect& dg);
	void GravGrinv(Vect& c, Vect& d,   //B9. Определение компонент гравитационного ускорения в гринвиче
		           Vect& x, Vect& dg);
	Vect& GravSphere(Vect& c, Vect& d, Vect& s);  //B10. Определение частных производных от геопотенциала
	void SpeedMoon(double t, Vect& x, Vect& dg); //B11. Лунное ускорение
	void SpeedSun(double t, Vect& x, Vect& dg);  //B12. Солнечное ускорение
	void SpeedAtm(double t, Vect& x, Vect& dg);  //B13. Атмосферное ускорение
	bool FactorLux(double t, Vect& x);			 //B14. Светотеневая обстановка
	void FactorLux(double t, Vect& x, double& dalfa); //B15. Светотеневая обстановка
	bool FactorLux(Vect& x, double t);				  //B16. Светотеневая обстановка
	void FactorLux(Vect& x, double t, double& deld);  //B17. Светотеневая обстановка
	Vect& PressLux(double t, Vect& x);				  //B18. Давление солнечных лучей
	double Atmos(double t, Vect& x); //B19. Динамическая модель атмосферы (ГОСТ Р 25645.166-2004)
	void Moon(double udt, Vect& r);  //B20. Определение координат Луны в J2000
	void Sun(double udt, Vect& r);	 //B21. Определение координат Солнца в J2000
	void Look();					 //B22. Просмотр компонент класса Ballist
	Ballist* operator->() { return this; } //B23. Доступ к компонентам класса Ballist
};

//Прототипы глобальных функций
double* Array(int m);				  //G1. Выделение памяти под одномерный массив
double** Array(int m, int n);		  //G2. Выделение памяти под двумерный массив
double*** Array(int m, int n, int l); //G3. Выделение памяти под трехмерный массив
Vect& Vector_e(int n, int k);		  //G4. Единичный вектор
Matr& MatrixI(int n);				  //G5. Единичная матрица
char* Line(int m);					  //G6. Выделение памяти под массив символов
char** Line(int m, int ns);			  //G7. Выделение памяти под одномерный массив [m] из ns символов
char*** Line(int m, int n, int ns);   //G8. Выделение памяти под двумерный массив [m x n x l] из ns символов
char**** Line(int m, int n, int l, int ns); //G9. Выделение памяти под двумерный массив [m x n x l] из ns символов
char Lit(int x);					  //G10. Приведение к цифре (если x - не цифра, а составное число)
char* Liter(double x, int m);		  //G11. Перевод числа x в строку из m символов
char* Liter(const char* c, int m);	  //G12. Перевод символьной строки с в строку из m символов
void Liter(double x, int m, char* c); //G13. Перевод числа x в строку из m символов
void Liter(const char* c, int m, char* s); //G14. Перевод символьной строки с в строку из m символов
Vect& Abs(Vect& x);						   //G15.
Matr& Abs(Matr& A);						   //G16.
Vect& MinFunc(double Func(Vect& x), int Nit, double eps, Vect& dx, Vect& x0); //G17. Метод конфигураций
Vect& Exp(Vect& x);						   //G18. Векторная экспонента
double Sum(Matr& A);					   //G19. Сумма всех элементов матрицы
Matr& Dmat(int n, double a);			   //G20. Формирование диагональной матрицы по скаляру
Matr& Dmat(Vect& a);					   //G21. Формирование диагональной матрицы по вектору
double Sum(Vect& a);					   //G22. Сумма всех элементов вектора
int RangMatr(Matr& A, Vect& s, Vect& c, double eps); //G23. Определение ранга матрицы
int RangMatr(Matr& A, double eps);					 //G24. Определение ранга матрицы
Matr& OwnVector(double F(double lam), int Nit,		 //G25. Определение собственных значений и 
	            double eps, Matr& A, Vect& lam);	 //собственных векторов матрицы A 
Vect& NotTrivial(double eps, Matr& A); //G26. Частное нетривиальное решение
									   //линейной однородной системы уравнений
int RangMatr(double eps, Matr& B, Vect& row, Vect& col); //G27. Определение ранга матрицы
int RangMatr(double eps, Matr& A);						 //G28. Определение ранга матрицы
void ZeroFun(double F(double t), int Nit, double eps1,	 //G29. Определение вещественных нулей функции
	         double eps2, double eps3, double eta, Vect& c);
void ZeroFun(double F(double t), int Nit,				 //G30. Определение вещественных нулей функции 
	         double eps, Vect& c);
void Sort1(Vect& x);	 //G31. Сортировка массива по возрастанию 
void Sort2(Vect& x);	 //G32. Сортировка массива по убыванию
Vect& Sort3(Vect& x);	 //G33. Сортировка массива по возрастанию с указанием номеров
Vect& Sort4(Vect& x);	 //G34. Сортировка массива по убыванию с указанием номеров
Vect& Sort5(Vect& x);	 //G35. Определение номеров компонент массива по возрастанию
Vect& Sort6(Vect& x);	 //G36. Определение номеров компонент массива по убыванию
bool TestZero(Vect& x);  //G37. Тестирование на ноль. Возвращает true, если вектор нулевой
bool TestZero(Matr& A);  //G38. Тестирование на ноль. Возвращает true, если матрица нулевая
double det(Matr& A);	 //G39. Определитель матрицы
Vect& Adding(Vect& x, Vect& y);			 //G40. Прямое сложение векторов
double Forme(Vect& x, Matr& P, Vect& y); //G41. Билинейная (квадратичная при x=y) форма
double Forme(Vect& x, Matr& P);			 //G42. Квадратичная форма
void Runge(double t, double T, Vect& x,  //G43. Интегрирование системы дифференциальных уравнений 
	       Vect& F(double t, Vect& x));
void Runge(Vect& F(double t, Vect& x),	 //G44. Интегрирование системы дифференциальных уравнений 
	       double dt, double t, Vect& x);
double Runge(Vect& F(double t, Vect& x), //G45. Интегрирование системы дифференциальных уравнений
	       double eps, double t, double T, Vect& x);
void Runge(Matr& F(double t, Matr& X),       //G46. Интегрирование матричной системы  
	       double t, double dt, Matr& X);    //дифференциальных уравнений
double Runge(Matr& F(double t, Matr& X), double eps, //G47. Интегрирование матричной системы  
	       double t, double T, Matr& X);             //дифференциальных уравнений
double Comp(double a, double b, double c);   //G48. Используется в Runge
int Expon(double x);						 //G49. Используется в Runge
Vect Lagrange(Vect& T, double t);			 //G50. Множители интерполяционного полинома Лагранжа
double Lagrange(int L, double t, Vect& T, Vect& F); //G51. Интерполяция полиномами Лагранжа
int Cnm(int n, int m);						  //G52. Количество сочетаний из n по m
void Comb(bool& prim, int n, int m, Vect& c); //G53. Генератор сочетаний
Matr& Comb(int n, int m);					  //G54. Генератор сочетаний
void emp();				//G55. Пропуск строки на экране 
void emp(int n);		//G56. Пропуск n строк на экране
Matr& SubMatr(Matr& B); //G57. Выделение и перемещение подматрицы полного ранга
char* Char(Str& x);		//G58. Преобразование Str->char
void Char(Str& x, char* s); //G59. Преобразование Str->char
Vect& VecProd(Vect& x, Vect& y); //G60. Векторное произведение векторов
Vect& AbsOsc(Vect& x);//G61. Переход от вектора x в J2000 к оскулирующим элементам орбиты
double Altitude(Vect& x); //G62. Высота над поверхностью референц-эллипсоида
double AngleSunOrb(double ud,                      //G63. Угол между плоскостью орбиты и  
	               double t, double i, double Om); //направлением на Солнца из центра Земли
double AngleSunOrb(double ud, double t, Vect& x);  //G64.
double AngleSunOrb(double ud, Vect& x, double t);  //G65.
double AnodePrg(double T, double e, double om); //G66. Время перелета от восходящего узла до перигея
double Anom(double mu, double eps,				//G67. Решение уравнения Кеплера 
	        double tau, double e, double a, double t);
double ArgWid(double eps, Vect& xos, double t); //G68. Определение аргумента широты
void Azim(double ud, double t,					//G69. Определение местного азимута трассы 
	      Vect& x, double& sinA, double& cosA);
double Asin(double x, double y);    //G70. Круговой арксинус
double Arcsin(double y, double x);  //G71. Круговой арксинус
double dAsin(double x, double y);   //G72. Производная от кругового арксинуса
double dArcsin(double y, double x); //G73. Производная от кругового арксинуса
double Arctg(double y, double x);   //G74. Круговой арктангенс
double dArctg(double y, double x);  //G75. Производная от кругового арктангенса
void ISA_IGI(int d, int m, int y,   //G76. К динамической модели атмосферы
	         double& F107, double& F81, double& Kp);
double FuncA(int d, int m, int y);  //G77. К динамической модели атмосферы
int SelectF0(double F81);			//G78. К динамической модели атмосферы
Vect& Soleil(double ud, double tg); //G79. Определение координат Солнца в J2000
double TimeLocal(double ud, double t, Vect& x); //G80. Определение местного солнечного времени
double TimeGrinv(double ud, double t, Vect& x); //G81. Определение гринвичского времени
void Date(double ud, int& d, int& m, int& y);	//G82. Перевод юлианской даты в календарную
Vect& Date(double ud);							//G83. Перевод юлианской даты в календарную
Vect& DMY(char* dat); //G84. Распаковка даты из формата 00.00.00 в вектор
Str& DMY(double ud);  //G85.Перевод юлианской даты в строку вида 00/00/0000
Vect& HMS(char* tim); //G86. Распаковка времени из формата 00:00:00.0000 в вектор
Str& HMS(double t);   //G87. Перевод секунд в строку вида 00:00:00.0000
Str& DHMS(double ud, double t); //G88. Перевод юлианской даты и времени в секундах
Vect& GeodGrinv(Vect& g); //G89. Переход от геодезических координат (g=||B,L,H||) к гринвичским
Vect& GrinvGeod(Vect& x); //G90. Переход от гринвичских координат к геодезическим (g=||B,L,H||)
double dGamma(double fip, double H); //G91. Поправка к минимальному углу возвышения над горизонтом
void dGamma(double B, double H,  //G92. Поправка к минимальному углу возвышения над горизонтом
	        double& dG);
void TimeRadio(double F(Vect& G, double e),    //G93. Определение времени входа и выхода из зоны 
	           Vect& xos, Vect& q, double eps, //радиовидимости наземного пункта 
	           double GamMin, double Lrn, double& T0, double& Tf);
void TimeRadio(double F(Vect& G, double e),    //G94. Определение времени входа и выхода из зоны
	           double eps, double GamMin,      //радиовидимости наземного пункта
	           double ud, double t0, double tf, Vect& xos, Vect& q,
	           int& count, Vect& tauP, Vect& T0, Vect& Tf);
double StimeVer(double ud, double t); //G95. Вычисление истинного звёздного времени
double Pers(Vect& x); //G96. Определение сидерического периода обращения по ВС в J2000
double Pers(double a); //G97. Определение сидерического периода обращения по большой полуоси
void Nutn(const double ud,			   //G98. Вычисление долго и короткопериодических членов   
	      double& hpsi, double& heps); //нутации истинного полюса в долготе и наклоне
double Eclipte(const double ud); //G99. Вычисление среднего наклона экватора к эклиптике
double Trace(Matr& A); //G100. След матрицы
bool BiSec(double Func(Vect& G, double t), Vect& G, double t0, //G101. Модификация BiSec для 
	       double tf, double epsf, double epst, double& t);		   //использования в функции TimeRadio
Vect& ElemOrb(double mu, double t, Vect& x); //G102. Определение оскулирующих элементов орбиты
double TimePrg(double mu, double tau, //G103. Определение времени прохождения углового расстояния, 
	           double e, double a, double teta); //равного истинной аномалии
Vect& Kepler(double mu, double eps, //G104. Переход от оскулирующих элементов орбиты 
	         double t, Vect& x);    //к прямоугольным инерциальным координатам
Vect& Kepler(double mu, double eps, //G105. Переход от оскулирующих элементов орбиты 
	         double t, Vect& x, double& u); //к прямоугольным инерциальным координатам
Matr& Reduct(double ud, //G106. Вычисление матрицы перехода от J2000 
	double t);	//к J-ТЭ с учётом прецессии и нутации
Matr& ReductG(double ud, //G107. Вычисление матрицы перехода от J2000 
	          double t); //к J-ТЭ с учётом прецессии и нутации
Matr& Reduct(double ud,    //G108. Вычисление истинного звёздного времени и матрицы перехода от J2000
	         double t, double& S); //к J-ТЭ с учётом прецессии и нутации
Matr& MatRot(int k, double u); //G109. Вычисление матрицы поворота относительно одной из заданных осей
Matr& Mat2000Gr(double ud, double t); //G110. Вычисление матрицы перехода от J2000 к ГСК
Matr& Mat2000GrV(double ud, //G111. Вычисление матрицы перехода от проекций скорости в J2000 
	             double t); //к проекциям скорости в ГСК
Matr& MatGr2000V(double ud, //G112. Вычисление матрицы перехода от проекций скорости в ГСК 
	             double t); //к проекциям скорости в J2000
Matr& Mat2000GrG(double ud, //G113. Вычисление матрицы перехода от вектора состояния в J2000 
	             double t); //к вектору состояния в ГСК
Matr MatGr2000G(double ud, //G114. Вычисление матрицы перехода от вектора состояния в ГСК 
	            double t); //к вектору состояния в J2000
Matr& MatPrec(double ud); //G115. Вычисление матрицы прецессии
Matr MatNut(double ud); //G116. Вычисление матрицы нутации
double Stime(double udt); //G117. Вычисление гринв. среднего звёздного времени на заданное время суток
double Snul(double ud); //G118. Вычисление средн. звёздного времени гринв. меридиана на начало суток
double Ulian(int d, int m, int y); //G119. Перевод календарной даты в юлианскую
double Uday(int d, int m, int y);  //G120. Перевод календарной даты в юлианскую
double Polynom(Vect& c, double x); //G121. Полином степени n-1 (порядка n) от x
void PolinomM(double x, Matr& a, Vect& f); //G122.
void PolinomS(double x, Matr& a, Vect& f); //G123.
Vect& DecSpher(Vect& d); //G124. Переход от декартовых координат к сферическим координатам
Vect& SpherDec(Vect& s); //G125. Переход от сферических координат к декартовым координатам
double Legendre(int n, double x); //G126. Вычисление значений полинома Лежандра порядка n
void Legendr(int n,				  //G127. Функция Legendr вычисляет значения  
	         double x, Matr& P);  //присоединенных функций Лежандра первого рода
Matr& MatDecSpher(double fi, //G128. Матрица перехода от декартовой к сферической системе координат 
	              double lam);
Matr MatDecSpher(Vect& x); //G129. Матрица перехода от декартовой к сферической системе координат
double StAtmos(double H); //G130. Статическая атмосфера (модель ABCD-59)
Vect& VectKinMom(double i, double Om); //G131. Вектор кинетического момента нормированный
void Moon(double udt, Vect& r);        //G132. Определение координат Луны в J2000
void Sun(double udt, Vect& r);         //G133. Определение координат Солнца в J2000
Matr& Mat2000Orb(Vect& xos);           //G134. Матрица перехода из с.к. J2000 в орбитальную систему координат
Matr& Mat2000OrbG(Vect& xos);          //G135. Матрица перехода из с.к. J2000 в орбитальную с.к.
Matr& MatAbsOrb(Vect& x);              //G136. Матрица перехода из с.к. J2000 в орбитальную с.к.
Matr& MatAbsOrbG(Vect& x); //G137. Матрица перехода из с.к. J2000 в орбитальную с.к.
Matr& MatOrbAbsG(Vect& x); //G138. Матрица перехода из орбитальной с.к. в с.к. J2000
void BulSto(Vect& F(double t, Vect& x), double eps, //G139. Интегрирование системы дифференциальных ур.
	        double& dt0, double t0, double tf, Vect& x);
double Rrandom(double& x); //G140. Получение равномерно распределённого
						   //на отрезке [0, 1] псевдослучайного числа
double Nrandom(double& x); //G141. Получение нормально распределённого псевдослучайного числа
						   //с нулевым мат. ожиданием и единичным с.к.о.						   
Vect& RandVec(double& rand, Vect& x0, Matr& P); //G142. Получение реализации случайного,
											    //нормально распределенного
											    //n-мерного вектора с математическим ожиданием x
											    //и матрицей ковариаций P
double StartRand(); //G143. Инициализация датчиков ПСВ с использованием текущего времени
Matr& DifFun(Vect& Func(Vect& x), int m, Vect& x, Vect& dx); //G144. Определение матрицы
															 //частных производных A(n,m)
															 //от функции Func размерности m 
															 //по аргументу x размерности n
Matr& DifFunT(Vect& Func(Vect& x), int m, Vect& x, Vect& dx); //G145. Определение транспонированной
															  //матрицы частных производных A(m,n)
															  //от функции Func размерности m 
															  //по аргументу x размерности n
Vect& GradF(double f(Vect& x), Vect& x0, Vect& dx); //G146. Определение градиента фунции f(x)
													//конечно-разностным методом
													//в точке x0 с приращением dx
double DifFun(double Func(double x), double x, double dx); //G147. Определение производной от
														   //скалярной функции по скалярному
														   //аргументу в точке x
double dnf_dxn(double f(double x), int n, double x, double dx); //G148. Определение производной
																//n-го порядка от скалярной функции
																//скалярного аргумента
																//(рекурсивная функция)
double BlancBruit(double& rand, double D, double dt);     //G149. Белый шум на dt
Vect& BlancBruit(bool& marker, double& rand, int N, Matr& D,
				 double t, double t1, double t2); //G150. Векторный белый шум
Matr& Zero(int m, int n); //G151. Нулевая матрица
Vect& Zero(int m);        //G152 Нулевой вектор
int Silvestr(Matr& P);	  //G153. Критерий Сильвестра
int Repositive(double eps, double gamma, Matr& P); //G154. Восстановление положительной определённости
void Jacobi(Matr& B, Matr& S, Vect& lambda, double eps); //G155. Собственные значения и векторы  
void Decomp(fstream& out, double eps, Matr& P, Matr& S, Vect& Lam); //G156. Разложение матрицы
Vect& VN(int n);				 //G157				
Vect& VS(int m, double s);		 //G158
Vect& VS(double s);				 //G159
void Invert(Matr& A);			 //G160. Перезапись строк матрицы в обратном порядке
Vect& Adding(Vect& x, double y); //G161. Прямое сложение вектора и скаляра
Vect& Adding(double y, Vect& x); //G162. Прямое сложение скаляра и вектора
Vect Adding(Vect& x, Vect& y, Vect& z);    //G163. Прямое сложение трёх векторов
void Scatter(Vect& x, Vect& x1, Vect& x2); //G164. Разложение вектора x на фрагменты x1, x2
void Scatter(Vect& x, Vect& x1, Vect& x2, Vect& x3); //G165. Разложение вектора x
													 //на фрагменты x1, x2, x3
long int Fact(int x);				//G166. Факториал (рекурсивная)
Vect& Sqrt(Vect& x);				//G167. Поэлементное извлечение корня из вектора
Matr& Sqrt(Matr& A);				//G168. Поэлементное извлечение корня из матрицы
Matr Cut(int I, Matr& P);			//G169. Удаляет из матрицы I-ю строку
Matr Cut(Matr& P, int J);			//G170 . Удаляет из матрицы J-й столбец
Matr Cut(Matr& P, int I, int J);	//G171. Удаляет из матрицы I-ю строку и J-й столбец
Vect Cut(Vect& p, int k);			//G172. Удаляет из вектора p k-й элемент

void DecBit(int x, char bit[]);
char* DecBit(int x);
char* DecBase(int x, int base);
char* DecBase(double x, int base);
int NumberBase(int x, int base);

//T1.
template <class Type> Type Max(Type x1, Type x2)
{
	return x1 > x2 ? x1 : x2;
}

//T2.
template <class Type> Type Max(Type x1, Type x2, Type x3)
{
	Type x = Max(x1, x2);
	return Max(x3, x);
}

//T3.
template <class Type> Type Min(Type x1, Type x2, Type x3)
{
	Type x = Min(x1, x2);
	return Min(x3, x);
}

//T4.
template <class Type> Type Min(Type x1, Type x2)
{
	return x1 < x2 ? x1 : x2;
}

//T5. Сигнатура
template <class T> int sign(T x)
{
	int sgn;
	sgn = x < 0 ? -1 : x>0 ? 1 : 0;
	return sgn;
}

//T6.
template <class T> void emp(T& out) { out << "\n"; }

//T7.
template <class T>void emp(T& out, int n) { for (int i = 1; i <= n; i++) out << "\n"; }

//T8.
template <class TT, class T>
void Ecrir(TT& output, T x) { output << x; }

//T9.
template <class TT, class T1, class T2>
void Ecrir(TT& output, T1 x1, T2 x2) { output << x1 << " " << x2; }

//T10.
template <class TT, class T1, class T2, class T3>
void Ecrir(TT& output, T1 x1, T2 x2, T3 x3) { output << x1 << " " << x2 << " " << x3; }

//T11.
template <class TT, class T1, class T2, class T3, class T4>
void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4;
}

//T12.
template <class TT, class T1, class T2, class T3, class T4, class T5>
void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5;
}

//T13.
template <class TT, class T1, class T2, class T3, class T4, class T5, class T6>
void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6;
}

//T14.
template <class TT, class T1, class T2, class T3, class T4, class T5, class T6, class T7>
void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7)
{
	output << x1 << " " << x2 << " " << x3 << " "
		<< x4 << " " << x5 << " " << x6 << " " << x7;
}

//T15.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3,
		T4 x4, T5 x5, T6 x6, T7 x7, T8 x8)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4
		<< " " << x5 << " " << x6 << " " << x7 << " " << x8;
}

//T16.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3,
		T4 x4, T5 x5, T6 x6, T7 x7, T8 x8, T9 x9)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4
		<< " " << x5 << " " << x6 << " " << x7 << " " << x8 << " " << x9;
}

//T17.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4,
		T5 x5, T6 x6, T7 x7, T8 x8, T9 x9, T10 x10)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5
		<< " " << x6 << " " << x7 << " " << x8 << " " << x9 << " " << x10;
}

//T18.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10 << " " << x11;
}

//T19.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11, class T12>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12;
}

//T20.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13;
}

//T21.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14;
}

//T22.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15;
}

//T23.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15, T16 x16)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15
		<< " " << x16;
}

//T24.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15,
		T16 x16, T17 x17)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15
		<< " " << x16 << " " << x17;
}

//T25.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16,
	class T17, class T18>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15,
		T16 x16, T17 x17, T18 x18)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15
		<< " " << x16 << " " << x17 << " " << x18;
}

//T26.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16,
	class T17, class T18, class T19>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15,
		T16 x16, T17 x17, T18 x18, T19 x19)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15
		<< " " << x16 << " " << x17 << " " << x18 << " " << x19;
}

//T27.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16,
	class T17, class T18, class T19, class T20>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15,
		T16 x16, T17 x17, T18 x18, T19 x19, T20 x20)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15
		<< " " << x16 << " " << x17 << " " << x18 << " " << x19 << " " << x20;
}

//T28.
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16,
	class T17, class T18, class T19, class T20, class T21>
	void Ecrir(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
		T9 x9, T10 x10, T11 x11, T12 x12, T13 x13, T14 x14, T15 x15,
		T16 x16, T17 x17, T18 x18, T19 x19, T20 x20, T21 x21)
{
	output << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6
		<< " " << x7 << " " << x8 << " " << x9 << " " << x10
		<< " " << x11 << " " << x12 << " " << x13 << " " << x14 << " " << x15
		<< " " << x16 << " " << x17 << " " << x18 << " " << x19 << " " << x20
		<< " " << x21;
}

//T29.
template <class TT, class T1>
void Look(TT& output, T1 x1)
{
	int j, q = 100; char* s1 = Line(q); Liter(x1, q, s1);
	for (j = 0; j < q; j++) output << s1[j];
	delete[]s1;
}

//T30.
template <class TT, class T1>
void Look(TT& output, T1 x1, int q)
{
	int j; char* s1 = Line(q); Liter(x1, q, s1);
	for (j = 0; j < q; j++) output << s1[j];
	delete[]s1;
}

//T31.
template <class TT, class T1, class T2>
void Look(TT& output, T1 x1, T2 x2, int q)
{
	int j; char* s2 = Line(q); Liter(x2, q, s2);
	Look(output, x1, q);
	output << " "; for (j = 0; j < q; j++) output << s2[j];
	delete[]s2;
}

//T32.
template <class TT, class T1, class T2, class T3>
void Look(TT& output, T1 x1, T2 x2, T3 x3, int q)
{
	int j; char* s3 = Line(q); Liter(x3, q, s3);
	Look(output, x1, x2, q);
	output << " "; for (j = 0; j < q; j++) output << s3[j];
	delete[]s3;
}

//T33.
template <class TT, class T1, class T2, class T3, class T4>
void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, int q)
{
	int j; char* s4 = Line(q); Liter(x4, q, s4);
	Look(output, x1, x2, x3, q);
	output << " "; for (j = 0; j < q; j++) output << s4[j];
	delete[]s4;
}

//T34.
template <class TT, class T1, class T2, class T3, class T4, class T5>
void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, int q)
{
	int j; char* s5 = Line(q); Liter(x5, q, s5);
	Look(output, x1, x2, x3, x4, q);
	output << " "; for (j = 0; j < q; j++) output << s5[j];
	delete[]s5;
}

//T35.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, int q)
{
	int j; char* s6 = Line(q); Liter(x6, q, s6);
	Look(output, x1, x2, x3, x4, x5, q);
	output << " "; for (j = 0; j < q; j++) output << s6[j];
	delete[]s6;
}

//T36.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, int q)
{
	int j; char* s7 = Line(q); Liter(x7, q, s7);
	Look(output, x1, x2, x3, x4, x5, x6, q);
	output << " "; for (j = 0; j < q; j++) output << s7[j];
	delete[]s7;
}

//T37.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, int q)
{
	int j; char* s8 = Line(q); Liter(x8, q, s8);
	Look(output, x1, x2, x3, x4, x5, x6, x7, q);
	output << " "; for (j = 0; j < q; j++)output << s8[j]; delete[]s8;
}

//T38.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, int q)
{
	int j; char* s9 = Line(q); Liter(x9, q, s9);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, q);
	output << " "; for (j = 0; j < q; j++)output << s9[j]; delete[]s9;
}

//T39.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, int q)
{
	int j; char* s10 = Line(q); Liter(x10, q, s10);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, q);
	output << " "; for (j = 0; j < q; j++)output << s10[j]; delete[]s10;
}

//T40.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11, int q)
{
	int j; char* s11 = Line(q); Liter(x11, q, s11);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, q);
	output << " "; for (j = 0; j < q; j++)output << s11[j]; delete[]s11;
}

//T41.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, int q)
{
	int j; char* s12 = Line(q); Liter(x12, q, s12);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, q);
	output << " "; for (j = 0; j < q; j++)output << s12[j]; delete[]s12;
}

//T42.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, int q)
{
	int j; char* s13 = Line(q); Liter(x13, q, s13);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, q);
	output << " "; for (j = 0; j < q; j++)output << s13[j]; delete[]s13;
}

//T43.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, int q)
{
	int j; char* s14 = Line(q); Liter(x14, q, s14);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, q);
	output << " "; for (j = 0; j < q; j++)output << s14[j]; delete[]s14;
}

//T44.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, int q)
{
	int j; char* s15 = Line(q); Liter(x15, q, s15);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, q);
	output << " "; for (j = 0; j < q; j++)output << s15[j]; delete[]s15;
}

//T45.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16, int q)
{
	int j; char* s16 = Line(q); Liter(x16, q, s16);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, q);
	output << " "; for (j = 0; j < q; j++)output << s16[j]; delete[]s16;
}

//T46.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, int q)
{
	int j; char* s17 = Line(q); Liter(x17, q, s17);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, q);
	output << " "; for (j = 0; j < q; j++)output << s17[j]; delete[]s17;
}

//T47.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, int q)
{
	int j; char* s18 = Line(q); Liter(x18, q, s18);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, q);
	output << " "; for (j = 0; j < q; j++)output << s18[j]; delete[]s18;
}

//T48.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18, class T19>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, T19 x19, int q)
{
	int j; char* s19 = Line(q); Liter(x19, q, s19);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,
		x18, q);
	output << " "; for (j = 0; j < q; j++)output << s19[j]; delete[]s19;
}

//T49.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18, class T19, class T20>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, T19 x19, T20 x20, int q)
{
	int j; char* s20 = Line(q); Liter(x20, q, s20);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,
		x18, x19, q);
	output << " "; for (j = 0; j < q; j++)output << s20[j]; delete[]s20;
}

//T50.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18, class T19,
	class T20, class T21>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, T19 x19, T20 x20, T21 x21, int q)
{
	int j; char* s21 = Line(q); Liter(x21, q, s21);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,
		x18, x19, x20, q);
	output << " "; for (j = 0; j < q; j++)output << s21[j]; delete[]s21;
}

//T51.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18, class T19,
	class T20, class T21, class T22>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, T19 x19, T20 x20, T21 x21,
		T22 x22, int q)
{
	int j; char* s22 = Line(q); Liter(x22, q, s22);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,
		x18, x19, x20, x21, q);
	output << " "; for (j = 0; j < q; j++)output << s22[j]; delete[]s22;
}

//T52.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18, class T19,
	class T20, class T21, class T22, class T23>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, T19 x19, T20 x20, T21 x21,
		T22 x22, T23 x23, int q)
{
	int j; char* s23 = Line(q); Liter(x23, q, s23);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,
		x18, x19, x20, x21, x22, q);
	output << " "; for (j = 0; j < q; j++)output << s23[j]; delete[]s23;
}

//T53.
template <class TT, class T1, class T2, class T3,
	class T4, class T5, class T6, class T7,
	class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15,
	class T16, class T17, class T18, class T19,
	class T20, class T21, class T22, class T23,
	class T24>
	void Look(TT& output, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5,
		T6 x6, T7 x7, T8 x8, T9 x9, T10 x10, T11 x11,
		T12 x12, T13 x13, T14 x14, T15 x15, T16 x16,
		T17 x17, T18 x18, T19 x19, T20 x20, T21 x21,
		T22 x22, T23 x23, T24 x24, int q)
{
	int j; char* s24 = Line(q); Liter(x24, q, s24);
	Look(output, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17,
		x18, x19, x20, x21, x22, x23, q);
	output << " "; for (j = 0; j < q; j++)output << s24[j]; delete[]s24;
}

//T54. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T> void VecOut(TT& output, S s, T& x, int q)
{
	int i = 0, m = x->m, ls = strlen(s), N; char* gap = Line(ls + 2);
	if (s != "")output << "\n " << s << "="; else output << "\n";
	N = 10; //if (q <= 6)N = 20; if (q > 6 && q <= 10)N = 15; if (q > 10 && q <= 13)N = 12;
	for (i = 1; i <= m; i++)
	{
		char* v = Line(q); Liter(x(i), q, v);
		output << v << " "; if (fmod(i, N) == 0) { output << "\n" << gap; }
		delete[]v;
	}
}

//T55. Вывод значений векторов x1,x2 в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1, S s2, T& x2, int q)
{
	VecOut(output, s1, x1, q);
	VecOut(output, s2, x2, q);
}

//T56. Вывод значений векторов x1,x2,x3 в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1, S s2, T& x2, S s3, T& x3, int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q);
}

//T57. Вывод значений векторов x1,x2,x3,x4 в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1, S s2, T& x2, S s3, T& x3, S s4, T& x4, int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
}

//T58. Вывод значений векторов x1,x2,x3,x4,x5 в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1,
	S s2, T& x2,
	S s3, T& x3,
	S s4, T& x4,
	S s5, T& x5,
	int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
	VecOut(output, s5, x5, q);
}

//T59. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1,
	S s2, T& x2,
	S s3, T& x3,
	S s4, T& x4,
	S s5, T& x5,
	S s6, T& x6,
	int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
	VecOut(output, s5, x5, q); VecOut(output, s6, x6, q);
}

//T60. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1,
	S s2, T& x2,
	S s3, T& x3,
	S s4, T& x4,
	S s5, T& x5,
	S s6, T& x6,
	S s7, T& x7,
	int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
	VecOut(output, s5, x5, q); VecOut(output, s6, x6, q);
	VecOut(output, s7, x7, q);
}

//T61. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1,
	S s2, T& x2,
	S s3, T& x3,
	S s4, T& x4,
	S s5, T& x5,
	S s6, T& x6,
	S s7, T& x7,
	S s8, T& x8,
	int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
	VecOut(output, s5, x5, q); VecOut(output, s6, x6, q);
	VecOut(output, s7, x7, q); VecOut(output, s8, x8, q);
}

//T62. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1,
	S s2, T& x2,
	S s3, T& x3,
	S s4, T& x4,
	S s5, T& x5,
	S s6, T& x6,
	S s7, T& x7,
	S s8, T& x8,
	S s9, T& x9,
	int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
	VecOut(output, s5, x5, q); VecOut(output, s6, x6, q);
	VecOut(output, s7, x7, q); VecOut(output, s8, x8, q);
	VecOut(output, s9, x9, q);
}

//T63. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class S, class T>
void VecOut(TT& output, S s1, T& x1,
	S s2, T& x2,
	S s3, T& x3,
	S s4, T& x4,
	S s5, T& x5,
	S s6, T& x6,
	S s7, T& x7,
	S s8, T& x8,
	S s9, T& x9,
	S s10, T& x10,
	int q)
{
	VecOut(output, s1, x1, q); VecOut(output, s2, x2, q);
	VecOut(output, s3, x3, q); VecOut(output, s4, x4, q);
	VecOut(output, s5, x5, q); VecOut(output, s6, x6, q);
	VecOut(output, s7, x7, q); VecOut(output, s8, x8, q);
	VecOut(output, s9, x9, q); VecOut(output, s10, x10, q);
}

//T64. Вывод значений вектора x в поток output
//q>=6 - число позиций (выводимых разрядов)
template <class TT, class T>
void VectOut(TT& output, T& x, int q)
{
	int i, j, m = x->m; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* v = Line(q); if (i <= m)Liter(x(i), q, v);
		for (j = 0; j < q; j++)output << v[j]; output << " ";
		if (i < m)output << "\n";
		delete[]v;
	}
}

//T65. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2>
void VectOut(TT& output, T1& x1, T2& x2, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m;
	m = m1 > m2 ? m1 : m2; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		for (j = 0; j < q; j++)output << s1[j]; output << " ";
		for (j = 0; j < q; j++)output << s2[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2;
	}
}

//T66. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3>
void VectOut(TT& output, T1& x1, T2& x2, T3& x3, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q), * s3 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3;
	}
}

//T67. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4>
void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q), * s3 = Line(q), * s4 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
	}
}

//T68. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5>
void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m, m5 = x5->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4; m = m > m5 ? m : m5; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q), * s3 = Line(q), * s4 = Line(q), * s5 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4; delete[]s5;
	}
}

//T69. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5, class T6>
void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q), * s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6;
	}
}

//T70. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m, m7 = x7->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q), * s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q), * s7 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7;
	}
}

//T71. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m, m7 = x7->m, m8 = x8->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
	}
}

//T72. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m, m7 = x7->m, m8 = x8->m, m9 = x9->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9;
	}
}

//T73. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m,
		m5 = x5->m, m6 = x6->m, m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10;
	}
}

//T74. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m, m5 = x5->m,
		m6 = x6->m, m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m, m11 = x11->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11;
	}
}

//T75. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11, T12& x12, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m, m11 = x11->m, m12 = x12->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
	}
}

//T76. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m, m11 = x11->m, m12 = x12->m, m13 = x13->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13;
	}
}

//T77. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14;
	}
}

//T78. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15, int q)
{
	int i = 0, j = 0, m = 0,
		m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15;
	}
}

//T79. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m, m4 = x4->m,
		m5 = x5->m, m6 = x6->m, m7 = x7->m, m8 = x8->m, m9 = x9->m, m10 = x10->m,
		m11 = x11->m, m12 = x12->m, m13 = x13->m, m14 = x14->m, m15 = x15->m, m16 = x16->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
	}
}

//T80. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, T17& x17, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m,
		m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m,
		m16 = x16->m, m17 = x17->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16;
	m = m > m17 ? m : m17; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q),
			* s17 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		if (i <= m17)Liter(x17(i), q, s17); else for (j = 0; j < q; j++)s17[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j]; output << " ";
		for (j = 0; j < q; j++) output << s17[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
		delete[]s17;
	}
}

//T81. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17,
	class T18>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, T17& x17, T18& x18, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m,
		m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m,
		m16 = x16->m, m17 = x17->m, m18 = x18->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16;
	m = m > m17 ? m : m17; m = m > m18 ? m : m18; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q),
			* s17 = Line(q), * s18 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		if (i <= m17)Liter(x17(i), q, s17); else for (j = 0; j < q; j++)s17[j] = ' ';
		if (i <= m18)Liter(x18(i), q, s18); else for (j = 0; j < q; j++)s18[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j]; output << " ";
		for (j = 0; j < q; j++) output << s17[j]; output << " ";
		for (j = 0; j < q; j++) output << s18[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
		delete[]s17; delete[]s18;
	}
}

//T82. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17,
	class T18, class T19>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, T17& x17, T18& x18, T19& x19, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m,
		m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m,
		m16 = x16->m, m17 = x17->m, m18 = x18->m,
		m19 = x19->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16;
	m = m > m17 ? m : m17; m = m > m18 ? m : m18; m = m > m19 ? m : m19;
	emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q),
			* s17 = Line(q), * s18 = Line(q),
			* s19 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		if (i <= m17)Liter(x17(i), q, s17); else for (j = 0; j < q; j++)s17[j] = ' ';
		if (i <= m18)Liter(x18(i), q, s18); else for (j = 0; j < q; j++)s18[j] = ' ';
		if (i <= m19)Liter(x19(i), q, s19); else for (j = 0; j < q; j++)s19[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j]; output << " ";
		for (j = 0; j < q; j++) output << s17[j]; output << " ";
		for (j = 0; j < q; j++) output << s18[j]; output << " ";
		for (j = 0; j < q; j++) output << s19[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
		delete[]s17; delete[]s18; delete[]s19;
	}
}

//T83. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17,
	class T18, class T19, class T20>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, T17& x17, T18& x18, T19& x19,
		T20& x20, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m,
		m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m,
		m16 = x16->m, m17 = x17->m, m18 = x18->m,
		m19 = x19->m, m20 = x20->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16;
	m = m > m17 ? m : m17; m = m > m18 ? m : m18; m = m > m19 ? m : m19;
	m = m > m20 ? m : m20; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q),
			* s17 = Line(q), * s18 = Line(q),
			* s19 = Line(q), * s20 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		if (i <= m17)Liter(x17(i), q, s17); else for (j = 0; j < q; j++)s17[j] = ' ';
		if (i <= m18)Liter(x18(i), q, s18); else for (j = 0; j < q; j++)s18[j] = ' ';
		if (i <= m19)Liter(x19(i), q, s19); else for (j = 0; j < q; j++)s19[j] = ' ';
		if (i <= m20)Liter(x20(i), q, s20); else for (j = 0; j < q; j++)s20[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j]; output << " ";
		for (j = 0; j < q; j++) output << s17[j]; output << " ";
		for (j = 0; j < q; j++) output << s18[j]; output << " ";
		for (j = 0; j < q; j++) output << s19[j]; output << " ";
		for (j = 0; j < q; j++) output << s20[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
		delete[]s17; delete[]s18; delete[]s19; delete[]s20;
	}
}

//T84. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17,
	class T18, class T19, class T20, class T21>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, T17& x17, T18& x18, T19& x19,
		T20& x20, T21& x21, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m,
		m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m,
		m16 = x16->m, m17 = x17->m, m18 = x18->m,
		m19 = x19->m, m20 = x20->m, m21 = x21->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16;
	m = m > m17 ? m : m17; m = m > m18 ? m : m18; m = m > m19 ? m : m19;
	m = m > m20 ? m : m20; m = m > m21 ? m : m21; emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q),
			* s17 = Line(q), * s18 = Line(q),
			* s19 = Line(q), * s20 = Line(q),
			* s21 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		if (i <= m17)Liter(x17(i), q, s17); else for (j = 0; j < q; j++)s17[j] = ' ';
		if (i <= m18)Liter(x18(i), q, s18); else for (j = 0; j < q; j++)s18[j] = ' ';
		if (i <= m19)Liter(x19(i), q, s19); else for (j = 0; j < q; j++)s19[j] = ' ';
		if (i <= m20)Liter(x20(i), q, s20); else for (j = 0; j < q; j++)s20[j] = ' ';
		if (i <= m21)Liter(x21(i), q, s21); else for (j = 0; j < q; j++)s21[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j]; output << " ";
		for (j = 0; j < q; j++) output << s17[j]; output << " ";
		for (j = 0; j < q; j++) output << s18[j]; output << " ";
		for (j = 0; j < q; j++) output << s19[j]; output << " ";
		for (j = 0; j < q; j++) output << s20[j]; output << " ";
		for (j = 0; j < q; j++) output << s21[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
		delete[]s17; delete[]s18; delete[]s19; delete[]s20;
		delete[]s21;
	}
}

//T85. Вывод значений векторов в поток output с разрядностью q
template <class TT, class T1, class T2, class T3, class T4, class T5,
	class T6, class T7, class T8, class T9, class T10, class T11,
	class T12, class T13, class T14, class T15, class T16, class T17,
	class T18, class T19, class T20, class T21, class T22>
	void VectOut(TT& output, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6,
		T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15,
		T16& x16, T17& x17, T18& x18, T19& x19,
		T20& x20, T21& x21, T22& x22, int q)
{
	int i = 0, j = 0, m = 0, m1 = x1->m, m2 = x2->m, m3 = x3->m,
		m4 = x4->m, m5 = x5->m, m6 = x6->m,
		m7 = x7->m, m8 = x8->m, m9 = x9->m,
		m10 = x10->m, m11 = x11->m, m12 = x12->m,
		m13 = x13->m, m14 = x14->m, m15 = x15->m,
		m16 = x16->m, m17 = x17->m, m18 = x18->m,
		m19 = x19->m, m20 = x20->m, m21 = x21->m, m22 = x22->m;
	m = m1 > m2 ? m1 : m2; m = m > m3 ? m : m3; m = m > m4 ? m : m4;
	m = m > m5 ? m : m5; m = m > m6 ? m : m6; m = m > m7 ? m : m7;
	m = m > m8 ? m : m8; m = m > m9 ? m : m9; m = m > m10 ? m : m10;
	m = m > m11 ? m : m11; m = m > m12 ? m : m12; m = m > m13 ? m : m13;
	m = m > m14 ? m : m14; m = m > m15 ? m : m15; m = m > m16 ? m : m16;
	m = m > m17 ? m : m17; m = m > m18 ? m : m18; m = m > m19 ? m : m19;
	m = m > m20 ? m : m20; m = m > m21 ? m : m21; m = m > m22 ? m : m22;
	emp(output);
	for (i = 1; i <= m; i++)
	{
		char* s1 = Line(q), * s2 = Line(q),
			* s3 = Line(q), * s4 = Line(q),
			* s5 = Line(q), * s6 = Line(q),
			* s7 = Line(q), * s8 = Line(q),
			* s9 = Line(q), * s10 = Line(q),
			* s11 = Line(q), * s12 = Line(q),
			* s13 = Line(q), * s14 = Line(q),
			* s15 = Line(q), * s16 = Line(q),
			* s17 = Line(q), * s18 = Line(q),
			* s19 = Line(q), * s20 = Line(q),
			* s21 = Line(q), * s22 = Line(q);
		if (i <= m1)Liter(x1(i), q, s1); else for (j = 0; j < q; j++)s1[j] = ' ';
		if (i <= m2)Liter(x2(i), q, s2); else for (j = 0; j < q; j++)s2[j] = ' ';
		if (i <= m3)Liter(x3(i), q, s3); else for (j = 0; j < q; j++)s3[j] = ' ';
		if (i <= m4)Liter(x4(i), q, s4); else for (j = 0; j < q; j++)s4[j] = ' ';
		if (i <= m5)Liter(x5(i), q, s5); else for (j = 0; j < q; j++)s5[j] = ' ';
		if (i <= m6)Liter(x6(i), q, s6); else for (j = 0; j < q; j++)s6[j] = ' ';
		if (i <= m7)Liter(x7(i), q, s7); else for (j = 0; j < q; j++)s7[j] = ' ';
		if (i <= m8)Liter(x8(i), q, s8); else for (j = 0; j < q; j++)s8[j] = ' ';
		if (i <= m9)Liter(x9(i), q, s9); else for (j = 0; j < q; j++)s9[j] = ' ';
		if (i <= m10)Liter(x10(i), q, s10); else for (j = 0; j < q; j++)s10[j] = ' ';
		if (i <= m11)Liter(x11(i), q, s11); else for (j = 0; j < q; j++)s11[j] = ' ';
		if (i <= m12)Liter(x12(i), q, s12); else for (j = 0; j < q; j++)s12[j] = ' ';
		if (i <= m13)Liter(x13(i), q, s13); else for (j = 0; j < q; j++)s13[j] = ' ';
		if (i <= m14)Liter(x14(i), q, s14); else for (j = 0; j < q; j++)s14[j] = ' ';
		if (i <= m15)Liter(x15(i), q, s15); else for (j = 0; j < q; j++)s15[j] = ' ';
		if (i <= m16)Liter(x16(i), q, s16); else for (j = 0; j < q; j++)s16[j] = ' ';
		if (i <= m17)Liter(x17(i), q, s17); else for (j = 0; j < q; j++)s17[j] = ' ';
		if (i <= m18)Liter(x18(i), q, s18); else for (j = 0; j < q; j++)s18[j] = ' ';
		if (i <= m19)Liter(x19(i), q, s19); else for (j = 0; j < q; j++)s19[j] = ' ';
		if (i <= m20)Liter(x20(i), q, s20); else for (j = 0; j < q; j++)s20[j] = ' ';
		if (i <= m21)Liter(x21(i), q, s21); else for (j = 0; j < q; j++)s21[j] = ' ';
		if (i <= m22)Liter(x22(i), q, s22); else for (j = 0; j < q; j++)s22[j] = ' ';
		for (j = 0; j < q; j++) output << s1[j]; output << " ";
		for (j = 0; j < q; j++) output << s2[j]; output << " ";
		for (j = 0; j < q; j++) output << s3[j]; output << " ";
		for (j = 0; j < q; j++) output << s4[j]; output << " ";
		for (j = 0; j < q; j++) output << s5[j]; output << " ";
		for (j = 0; j < q; j++) output << s6[j]; output << " ";
		for (j = 0; j < q; j++) output << s7[j]; output << " ";
		for (j = 0; j < q; j++) output << s8[j]; output << " ";
		for (j = 0; j < q; j++) output << s9[j]; output << " ";
		for (j = 0; j < q; j++) output << s10[j]; output << " ";
		for (j = 0; j < q; j++) output << s11[j]; output << " ";
		for (j = 0; j < q; j++) output << s12[j]; output << " ";
		for (j = 0; j < q; j++) output << s13[j]; output << " ";
		for (j = 0; j < q; j++) output << s14[j]; output << " ";
		for (j = 0; j < q; j++) output << s15[j]; output << " ";
		for (j = 0; j < q; j++) output << s16[j]; output << " ";
		for (j = 0; j < q; j++) output << s17[j]; output << " ";
		for (j = 0; j < q; j++) output << s18[j]; output << " ";
		for (j = 0; j < q; j++) output << s19[j]; output << " ";
		for (j = 0; j < q; j++) output << s20[j]; output << " ";
		for (j = 0; j < q; j++) output << s21[j]; output << " ";
		for (j = 0; j < q; j++) output << s22[j];
		if (i < m)output << "\n";
		delete[]s1; delete[]s2; delete[]s3; delete[]s4;
		delete[]s5; delete[]s6; delete[]s7; delete[]s8;
		delete[]s9; delete[]s10; delete[]s11; delete[]s12;
		delete[]s13; delete[]s14; delete[]s15; delete[]s16;
		delete[]s17; delete[]s18; delete[]s19; delete[]s20;
		delete[]s21; delete[]s22;
	}
}

//T86. Вывод двумерного массива в поток out,
//|size|>=6 - число позиций, выделяемых под число->
//Если size<0, то имя массива указывается с каждым элементом->
//Иначе имя массива указывается один раз перед всеми его элементами,
//и после записи элементов точка с запятой не пишется
template <class TT, class S> void MatrOut(TT& out, S s, Matr& A, int size)
{
	int Size = abs(size), m = A->m, n = A->n, i, j0, j, k, kmx, r, ncol, M, mn, ij, d;
	ncol = 7; if (Size <= 6)ncol = 15; if (Size > 6 && Size <= 9)ncol = 10;
	if (Size > 9 && Size <= 12)ncol = 8; M = Min(n, ncol);
	kmx = n / M; r = n % M; if (size > 0)out << "\n " << s << ":\n"; else out << "\n";
	mn = (int)(log10(m) + 1) + (int)(log10(n) + 1);
	for (k = 1; k <= kmx; k++)
	{
		j0 = (k - 1) * M;
		for (i = 1; i <= m; i++)
		{
			for (j = j0 + 1; j <= M + j0; j++)
			{
				ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 2;
				char* a = Line(Size + d); Liter(A(i, j), Size + d, a);
				if (size > 0)out << " (" << i << "," << j << ")=" << a << " ";
				if (size < 0)out << " " << s << "(" << i << "," << j << ")=" << a << "; ";
				delete[]a;
			}
			out << "\n";
		}
		if (k < kmx || r>0)out << "\n";
	}
	if (r == 0)return;
	for (i = 1; i <= m; i++)
	{
		for (k = j; k < j + r; k++)
		{
			ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 2;
			char* a = Line(Size + d); Liter(A(i, k), Size + d, a);
			if (size > 0)out << " (" << i << "," << k << ")=" << a << " ";
			if (size < 0)out << " " << s << "(" << i << "," << k << ")=" << a << " ";
			delete[]a;
		}
		out << "\n";
	}
}

//T87. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT, class S>
void MatrOut(TT& out, S s1, Matr& A1, S s2, Matr& A2, int size)
{
	MatrOut(out, s1, A1, size);
	MatrOut(out, s2, A2, size);
}

//T88. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT, class S>
void MatrOut(TT& out, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, int size)
{
	MatrOut(out, s1, A1, size);
	MatrOut(out, s2, A2, size);
	MatrOut(out, s3, A3, size);
}

//T89. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT, class S>
void MatrOut(TT& out, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, S s4, Matr& A4, int size)
{
	MatrOut(out, s1, A1, size);
	MatrOut(out, s2, A2, size);
	MatrOut(out, s3, A3, size);
	MatrOut(out, s4, A4, size);
}

//T90. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT, class S>
void MatrOut(TT& out, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, S s4, Matr& A4, S s5, Matr& A5, int size)
{
	MatrOut(out, s1, A1, size);
	MatrOut(out, s2, A2, size);
	MatrOut(out, s3, A3, size);
	MatrOut(out, s4, A4, size);
	MatrOut(out, s5, A5, size);
}

//T91. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT, class S>
void MatrOut(TT& out, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, S s4, Matr& A4, S s5, Matr& A5,
	S s6, Matr& A6, int size)
{
	MatrOut(out, s1, A1, size);
	MatrOut(out, s2, A2, size);
	MatrOut(out, s3, A3, size);
	MatrOut(out, s4, A4, size);
	MatrOut(out, s5, A5, size);
	MatrOut(out, s6, A6, size);
}

//T92. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT, class S>
void MatrOut(TT& out, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, S s4, Matr& A4, S s5, Matr& A5,
	S s6, Matr& A6, S s7, Matr& A7, int size)
{
	MatrOut(out, s1, A1, size);
	MatrOut(out, s2, A2, size);
	MatrOut(out, s3, A3, size);
	MatrOut(out, s4, A4, size);
	MatrOut(out, s5, A5, size);
	MatrOut(out, s6, A6, size);
	MatrOut(out, s7, A7, size);
}

//T93. Вывод двумерного массива в поток out без указания имени массива,
//|size|>=6 - число позиций, выделяемых под число
template <class TT> void MatrOut(TT& out, Matr& A, int size)
{
	int Size = abs(size), m = A->m, n = A->n, i, j0, j, k,
		kmx, r, ncol = 9, M = Min(n, ncol), mn, ij, d;
	kmx = n / M; r = n % M;	mn = (int)(log10(m) + 1) + (int)(log10(n) + 1);
	for (k = 1; k <= kmx; k++)
	{
		j0 = (k - 1) * M;
		for (i = 1; i <= m; i++)
		{
			for (j = j0 + 1; j <= M + j0; j++)
			{
				ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
				char* a = Line(Size + d); Liter(A(i, j), Size + d, a);
				if (size > 0)out << " (" << i << "," << j << ")=" << a << " ";
				if (size < 0)out << " " << "(" << i << "," << j << ")=" << a << "; ";
				delete[]a;
			}
			out << "\n";
		}
		if (k < kmx || r>0)out << "\n";
	}
	if (r == 0)return;
	for (i = 1; i <= m; i++)
	{
		for (k = j; k < j + r; k++)
		{
			ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
			char* a = Line(Size + d); Liter(A(i, k), Size + d, a);
			if (size > 0)out << " (" << i << "," << k << ")=" << a << " ";
			if (size < 0)out << " " << "(" << i << "," << k << ")=" << a << " ";
			delete[]a;
		}
		out << "\n";
	}
}

//T94. Вывод двумерного массива в поток out,
//|size|>=6 - число позиций, выделяемых под число->
//Если size<0, то имя массива указывается с каждым элементом,
//и после каждого элемента пишется точка с запятой->
//Иначе имя массива указывается один раз перед всеми его элементами,
//и после записи элементов точка с запятой не пишется->
//ncol - количество столбцов, выводимых до переноса
template <class TT, class S> void MatrOut(TT& out, S s, Matr& A, int size, int ncol)
{
	int Size = abs(size), m = A->m, n = A->n, i, j0, j, k,
		kmx, r, M = Min(n, ncol), mn, ij, d;
	kmx = n / M; r = n % M; if (size > 0)out << "\n " << s << ":\n"; else out << "\n";
	mn = (int)(log10(m) + 1) + (int)(log10(n) + 1);
	for (k = 1; k <= kmx; k++)
	{
		j0 = (k - 1) * M;
		for (i = 1; i <= m; i++)
		{
			for (j = j0 + 1; j <= M + j0; j++)
			{
				ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
				char* a = Line(Size + d); Liter(A(i, j), Size + d, a);
				if (size > 0)out << " (" << i << "," << j << ")=" << a << " ";
				if (size < 0)out << " " << s << "(" << i << "," << j << ")=" << a << "; ";
				delete[]a;
			}
			out << "\n";
		}
		if (k < kmx || r>0)out << "\n";
	}
	if (r == 0)return;
	for (i = 1; i <= m; i++)
	{
		for (k = j; k < j + r; k++)
		{
			ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
			char* a = Line(Size + d); Liter(A(i, k), Size + d, a);
			if (size > 0)out << " (" << i << "," << k << ")=" << a << " ";
			if (size < 0)out << " " << s << "(" << i << "," << k << ")=" << a << "; ";
			delete[]a;
		}
		out << "\n";
	}
}

//T95. Вывод двумерного массива в поток out без индексов,
//size>=6 - число позиций, выделяемых под число
template <class TT> void MatrOut(TT& out, int size, Matr& A)
{
	int m = A->m, n = A->n, i, j0, j, k, kmx, r, ncol = 9, M = Min(n, ncol), mn, ij, d;
	kmx = n / M; r = n % M;	mn = (int)(log10(m) + 1) + (int)(log10(n) + 1);
	for (k = 1; k <= kmx; k++)
	{
		j0 = (k - 1) * M;
		for (i = 1; i <= m; i++)
		{
			for (j = j0 + 1; j <= M + j0; j++)
			{
				ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
				char* a = Line(size + d); Liter(A(i, j), size + d, a);
				out << a << " "; delete[]a;
			}
			out << "\n";
		}
		if (k < kmx || r>0)out << "\n";
	}
	if (r == 0)return;
	for (i = 1; i <= m; i++)
	{
		for (k = j; k < j + r; k++)
		{
			ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
			char* a = Line(size + d); Liter(A(i, k), size + d, a);
			out << a << " "; delete[]a;
		}
		out << "\n";
	}
}

//T96. Вывод двумерных массивов в поток out без индексов,
//size>=6 - число позиций, выделяемых под число
template <class TT>
void MatrOut(TT& out, int size, Matr& A1, Matr& A2)
{
	MatrOut(out, size, A1); emp(out);
	MatrOut(out, size, A2);
}

//T97. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT>
void MatrOut(TT& out, int size, Matr& A1, Matr& A2, Matr& A3)
{
	MatrOut(out, size, A1); emp(out);
	MatrOut(out, size, A2); emp(out);
	MatrOut(out, size, A3);
}

//T98. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT>
void MatrOut(TT& out, int size, Matr& A1,
	Matr& A2, Matr& A3, Matr& A4)
{
	MatrOut(out, size, A1); emp(out);
	MatrOut(out, size, A2); emp(out);
	MatrOut(out, size, A3); emp(out);
	MatrOut(out, size, A4);
}

//T99. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT>
void MatrOut(TT& out, int size, Matr& A1, Matr& A2,
	Matr& A3, Matr& A4, Matr& A5)
{
	MatrOut(out, size, A1); emp(out);
	MatrOut(out, size, A2); emp(out);
	MatrOut(out, size, A3); emp(out);
	MatrOut(out, size, A4); emp(out);
	MatrOut(out, size, A5);
}

//T100. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT>
void MatrOut(TT& out, int size, Matr& A1, Matr& A2,
	Matr& A3, Matr& A4, Matr& A5, Matr& A6)
{
	MatrOut(out, size, A1); emp(out);
	MatrOut(out, size, A2); emp(out);
	MatrOut(out, size, A3); emp(out);
	MatrOut(out, size, A4); emp(out);
	MatrOut(out, size, A5); emp(out);
	MatrOut(out, size), A6;
}

//T101. Вывод двумерных массивов в поток out,
//size>=6 - число позиций, выделяемых под число
template <class TT>
void MatrOut(TT& out, int size, Matr& A1, Matr& A2,
	Matr& A3, Matr& A4, Matr& A5, Matr& A6, Matr& A7)
{
	MatrOut(out, size, A1); emp(out);
	MatrOut(out, size, A2); emp(out);
	MatrOut(out, size, A3); emp(out);
	MatrOut(out, size, A4); emp(out);
	MatrOut(out, size, A5); emp(out);
	MatrOut(out, size, A6); emp(out);
	MatrOut(out, size, A7);
}

//T102. Вывод двумерного массива в поток out без индексов,
//size>=6 - число позиций, выделяемых под число->
template <class TT> void MatrOut(TT& out, int size, int ncol, Matr& A)
{
	int m = A->m, n = A->n, i, j0, j, k, kmx, r, M = Min(n, ncol), mn, ij, d;
	kmx = n / M; r = n % M; mn = (int)(log10(m) + 1) + (int)(log10(n) + 1);
	for (k = 1; k <= kmx; k++)
	{
		j0 = (k - 1) * M;
		for (i = 1; i <= m; i++)
		{
			for (j = j0 + 1; j <= M + j0; j++)
			{
				ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
				char* a = Line(size + d); Liter(A(i, j), size + d, a);
				out << a << " "; delete[]a;
			}
			out << "\n";
		}
		if (k < kmx || r>0)out << "\n";
	}
	if (r == 0)return;
	for (i = 1; i <= m; i++)
	{
		for (k = j; k < j + r; k++)
		{
			ij = (int)(log10(i) + 1) + (int)(log10(j) + 1); d = mn - ij - 3;
			char* a = Line(size + d); Liter(A(i, k), size + d, a);
			out << a << " "; delete[]a;
		}
		out << "\n";
	}
}

//T103. Чтение данных из потока entry и запись их в двумерный массив
template <class TT> void ReadMatr(TT& entry, Matr& A)
{
	int i, j, m = A->m, n = A->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A(i, j);
}

//T104. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT> void ReadMatr(TT& entry, Matr& A1, Matr& A2)
{
	int i, j, m = A1->m, n = A1->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A1(i, j);
	m = A2->m; n = A2->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
}

//T105. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT> void ReadMatr(TT& entry, Matr& A1, Matr& A2, Matr& A3)
{
	int i, j, m = A1->m, n = A1->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A1(i, j);
	m = A2->m; n = A2->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
	m = A3->m; n = A3->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A3(i, j);
}

//T106. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT>
void ReadMatr(TT& entry, Matr& A1, Matr& A2, Matr& A3, Matr& A4)
{
	int i, j, m = A1->m, n = A1->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; i++) entry >> A1(i, j);
	m = A2->m; n = A2->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
	m = A3->m; n = A3->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A3(i, j);
	m = A4->m; n = A4->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A4(i, j);
}

//T107. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT>
void ReadMatr(TT& entry, Matr& A1, Matr& A2, Matr& A3, Matr& A4, Matr& A5)
{
	int i, j, m = A1->m, n = A1->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; i++) entry >> A1(i, j);
	m = A2->m; n = A2->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
	m = A3->m; n = A3->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A3(i, j);
	m = A4->m; n = A4->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A4(i, j);
	m = A5->m; n = A5->n;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A5(i, j);
}

//T108. Чтение данных из потока entry и запись их в двумерный массив
template <class TT, class S>
void ReadMatr(TT& entry, S s, Matr& A)
{
	int i, j, m = A->m, n = A->n;
	entry >> s;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A(i, j);
}

//T109. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT, class S>
void ReadMatr(TT& entry, S s1, Matr& A1, S s2, Matr& A2)
{
	int i, j, m = A1->m, n = A1->n;
	entry >> s1;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A1(i, j);
	m = A2->m; n = A2->n;
	entry >> s2;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
}

//T110. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT, class S>
void ReadMatr(TT& entry, S s1, Matr& A1, S s2, Matr& A2, S s3, Matr& A3)
{
	int i, j, m = A1->m, n = A1->n;
	entry >> s1;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A1(i, j);
	m = A2->m; n = A2->n;
	entry >> s2;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
	m = A3->m; n = A3->n;
	entry >> s3;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A3(i, j);
}

//T111. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT, class S>
void ReadMatr(TT& entry, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, S s4, Matr& A4)
{
	int i, j, m = A1->m, n = A1->n;
	entry >> s1;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A1(i, j);
	m = A2->m; n = A2->n; entry >> s2;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
	m = A3->m; n = A3->n; entry >> s3;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A3(i, j);
	m = A4->m; n = A4->n; entry >> s4;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A4(i, j);
}

//T112. Чтение данных из потока entry и запись их в двумерные массивы
template <class TT, class S>
void ReadMatr(TT& entry, S s1, Matr& A1, S s2, Matr& A2,
	S s3, Matr& A3, S s4, Matr& A4, S s5, Matr& A5)
{
	int i, j, m = A1->m, n = A1->n;
	entry >> s1;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A1(i, j);
	m = A2->m; n = A2->n; entry >> s2;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A2(i, j);
	m = A3->m; n = A3->n; entry >> s3;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A3(i, j);
	m = A4->m; n = A4->n; entry >> s4;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A4(i, j);
	m = A5->m; n = A5->n; entry >> s5;
	for (i = 1; i <= m; i++) for (j = 1; j <= n; j++) entry >> A5(i, j);
}

//T113. Чтение данных из потока entry и запись их в одномерный массив
template <class TT>void ReadVect(TT& entry, Vect& a)
{
	int i, m = a->m;
	for (i = 1; i <= m; i++) entry >> a(i);
}

//T114. Чтение данных из потока entry и запись их в одномерный массив
template <class TT, class S>void ReadVect(TT& entry, S s, Vect& a)
{
	int i, m = a->m;
	entry >> s;
	for (i = 1; i <= m; i++) entry >> a(i);
}

//T115. Чтение данных из потока entry и запись их в одномерные массивы
template <class TT, class S>
void ReadVect(TT& entry, S s1, Vect& a1, S s2, Vect& a2)
{
	int i, m = a1->m;
	entry >> s1;
	for (i = 1; i <= m; i++) entry >> a1(i);
	m = a2->m; entry >> s2;
	for (i = 1; i <= m; i++) entry >> a2(i);
}

//T116. Чтение данных из потока entry и запись их в одномерные массивы
template <class TT, class S>
void ReadVect(TT& entry, S s1, Vect& a1, S s2, Vect& a2, S s3, Vect& a3)
{
	int i, m = a1->m; entry >> s1;
	for (i = 1; i <= m; i++) entry >> a1(i);
	m = a2->m; entry >> s2;
	for (i = 1; i <= m; i++) entry >> a2(i);
	m = a3->m; entry >> s3;
	for (i = 1; i <= m; i++) entry >> a3(i);
}

//T117. Чтение данных из потока entry и запись их в одномерные массивы
template <class TT, class S>
void ReadVect(TT& entry, S s1, Vect& a1, S s2, Vect& a2,
	S s3, Vect& a3, S s4, Vect& a4)
{
	int i, m = a1->m; entry >> s1;
	for (i = 1; i <= m; i++) entry >> a1(i);
	m = a2->m; entry >> s2;
	for (i = 1; i <= m; i++) entry >> a2(i);
	m = a3->m; entry >> s3;
	for (i = 1; i <= m; i++) entry >> a3(i);
	m = a4->m; entry >> s4;
	for (i = 1; i <= m; i++) entry >> a4(i);
}

//T118. Чтение данных из потока entry и запись их в одномерные массивы
template <class TT, class S>
void ReadVect(TT& entry, S s1, Vect& a1, S s2, Vect& a2,
	S s3, Vect& a3, S s4, Vect& a4, S s5, Vect& a5)
{
	int i, m = a1->m; entry >> s1;
	for (i = 1; i <= m; i++) entry >> a1(i);
	m = a2->m; entry >> s2;
	for (i = 1; i <= m; i++) entry >> a2(i);
	m = a3->m; entry >> s3;
	for (i = 1; i <= m; i++) entry >> a3(i);
	m = a4->m; entry >> s4;
	for (i = 1; i <= m; i++) entry >> a4(i);
	m = a5->m; entry >> s5;
	for (i = 1; i <= m; i++) entry >> a5(i);
}

//T119. 
template <class TT, class T>
void Read(TT& entry, T& x) { entry >> x; }

//T120.
template <class TT, class T1, class T2>
void Read(TT& entry, T1& x1, T2& x2)
{
	entry >> x1 >> x2;
}

//T121.
template <class TT, class T1, class T2, class T3>
void Read(TT& entry, T1& x1, T2& x2, T3& x3)
{
	entry >> x1 >> x2 >> x3;
}

//T122.
template <class TT, class T1, class T2, class T3, class T4>
void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4)
{
	entry >> x1 >> x2 >> x3 >> x4;
}

//T123.
template <class TT, class T1, class T2, class T3, class T4, class T5>
void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5;
}

//T124.
template <class TT, class T1, class T2, class T3, class T4, class T5, class T6>
void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
}

//T125.
template <class TT, class T1, class T2, class T3, class T4, class T5, class T6, class T7>
void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5, T6& x6, T7& x7)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7;
}

//T126.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8;
}

//T127.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9;
}

//T128.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9, class T10>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10;
}

//T129.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10 >> x11;
}

//T130.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11, T12& x12)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10 >> x11 >> x12;
}

//T131.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> x10 >> x11 >> x12 >> x13;
}

//T132.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14;
}

//T133.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14, class T15>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14 >> x15;
}

//T134.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14,
	class T15, class T16>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15, T16& x16)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14 >> x15 >> x16;
}

//T135.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14,
	class T15, class T16, class T17>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15, T16& x16, T17& x17)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14 >> x15 >> x16 >> x17;
}

//T136.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14,
	class T15, class T16, class T17, class T18>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15, T16& x16, T17& x17, T18& x18)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14 >> x15 >> x16 >> x17 >> x18;
}

//T137.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14,
	class T15, class T16, class T17, class T18, class T19>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15, T16& x16,
		T17& x17, T18& x18, T19& x19)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14 >> x15 >> x16 >> x17 >> x18 >> x19;
}

//T138.
template <class TT, class T1, class T2, class T3, class T4,
	class T5, class T6, class T7, class T8, class T9,
	class T10, class T11, class T12, class T13, class T14,
	class T15, class T16, class T17, class T18, class T19, class T20>
	void Read(TT& entry, T1& x1, T2& x2, T3& x3, T4& x4,
		T5& x5, T6& x6, T7& x7, T8& x8, T9& x9, T10& x10, T11& x11,
		T12& x12, T13& x13, T14& x14, T15& x15, T16& x16,
		T17& x17, T18& x18, T19& x19, T19& x20)
{
	entry >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >>
		x10 >> x11 >> x12 >> x13 >> x14 >> x15 >> x16 >> x17 >> x18 >> x19 >> x20;
}

//D1.
#define	MSTACK(m,n,C) static bool flag=1; \
					  static Matr* addr; Matr& C=*new Matr(m,n); \
					  if(flag){addr=&C; flag=0;} else{free((void*)addr); addr=&C;}

//D2.
#define	VSTACK(m,x) static bool flag=1; \
					static Vect* addr; Vect& x=*new Vect(m); \
                    if(flag){addr=&x; flag=0;} else{free((void*)addr); addr=&x;}
//D3.
#define STACKV(m) h = new Vect(m); NumbVect++; //Только для методов классов
											   //с компонентой Vect* h  
//D4.
#define STACKM(m, n) H = new Matr(m, n); NumbMatr++; //Только для методов классов
													 //с компонентой Matr* H
//D5.
#define STACKh(m, x) x->h = new Vect(m); x->NumbVect++; //Для методов и функций
														//с параметром типа Vect
//D6.
#define STACKH(m, n, A) A->H = new Matr(m, n); A->NumbMatr++; //Для методов и функций
															  //с параметром типа Matr
//D7.
#define STACKS(m) h = new Str(m); NumbStr++;

//D8.
#define STACKSh(m, s) s->h = new Str(m); s->NumbStr++;

//DD1.
#define DO(t,t1,t2) for(t=t1; t<=t2; t++)

//DD2.
#define Do(t,t1,t2,dt) for(t=t1; t<=t2; t+=dt)

//DD3.
#define dO(t,t1,t2,dt) for(t=t1; t>=t2; t+=dt)

//DD4.
#define EQ(x,x0) x>x0-1e-9 && x<x0+1e-9

//DD5.
#define END emp(3); exit(1);

//DD6.
#define EC_ Ecrir(cout,

//DD7.
#define LK_ Look(cout,

//DD8.
#define VC_ VecOut(cout,

//DD9.
#define MT_ MatrOut(cout,

//DD10.
#define RD_ Read(cin,

//DD11.
#define O_ emp();

//DD12.
#define $V(x) *(x->h)

//DD13.
#define $M(A) *(A->H)

//DD14. Для StatCharPol (модификация входной функции на случай определения СХ ее ОЛ)
#define FUNC(x,f) f=Func(x); if(type)f=f-fx0-dfdxT*(x-x0);
