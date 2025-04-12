#include "/Users/elizavetanosova/Documents/Arsenal.h"

fstream out("OutMatLib.txt", ios::out|ios::in|ios::trunc);

//Опрделение входного потока для функции ModProgn:
fstream ent("d:\\Arsenal\\Data\\CoefGravPot.txt", ios::in);
//Опрделение входных потоков для функции Atmos:
fstream multA("d:\\Arsenal\\Data\\multA.txt", ios::in);
fstream A120("d:\\Arsenal\\Data\\A120.txt", ios::in);
fstream A500("d:\\Arsenal\\Data\\A500.txt", ios::in);
fstream B120("d:\\Arsenal\\Data\\B120.txt", ios::in);
fstream B600("d:\\Arsenal\\Data\\B600.txt", ios::in);
fstream B660("d:\\Arsenal\\Data\\B660.txt", ios::in);
fstream B760("d:\\Arsenal\\Data\\B760.txt", ios::in);
fstream B800("d:\\Arsenal\\Data\\B800.txt", ios::in);
fstream B860("d:\\Arsenal\\Data\\B860.txt", ios::in);
fstream B900("d:\\Arsenal\\Data\\B900.txt", ios::in);
fstream B1000("d:\\Arsenal\\Data\\B1000.txt", ios::in);
fstream C120("d:\\Arsenal\\Data\\C120.txt", ios::in);
fstream C640("d:\\Arsenal\\Data\\C640.txt", ios::in);
fstream C700("d:\\Arsenal\\Data\\C700.txt", ios::in);
fstream C760("d:\\Arsenal\\Data\\C760.txt", ios::in);
fstream C820("d:\\Arsenal\\Data\\C820.txt", ios::in);
fstream C860("d:\\Arsenal\\Data\\C860.txt", ios::in);
fstream C920("d:\\Arsenal\\Data\\C920.txt", ios::in);
fstream C980("d:\\Arsenal\\Data\\C980.txt", ios::in);
fstream D120("d:\\Arsenal\\Data\\D120.txt", ios::in);
fstream E120("d:\\Arsenal\\Data\\E120.txt", ios::in);
fstream E600("d:\\Arsenal\\Data\\E600.txt", ios::in);
fstream E700("d:\\Arsenal\\Data\\E700.txt", ios::in);
fstream E760("d:\\Arsenal\\Data\\E760.txt", ios::in);
fstream E780("d:\\Arsenal\\Data\\E780.txt", ios::in);
fstream E800("d:\\Arsenal\\Data\\E800.txt", ios::in);
fstream E900("d:\\Arsenal\\Data\\E900.txt", ios::in);
fstream L120("d:\\Arsenal\\Data\\L120.txt", ios::in);
fstream L640("d:\\Arsenal\\Data\\L640.txt", ios::in);
fstream L660("d:\\Arsenal\\Data\\L660.txt", ios::in);
fstream L740("d:\\Arsenal\\Data\\L740.txt", ios::in);
fstream L800("d:\\Arsenal\\Data\\L800.txt", ios::in);
fstream L860("d:\\Arsenal\\Data\\L860.txt", ios::in);
fstream L900("d:\\Arsenal\\Data\\L900.txt", ios::in);
fstream N120("d:\\Arsenal\\Data\\N120.txt", ios::in);
fstream Fi1120("d:\\Arsenal\\Data\\Fi1120.txt", ios::in);
//Опрделение входных потоков для функций Moon, Sun:
fstream InMoon("d:\\Arsenal2019New\\Data\\KFMOON.txt", ios::in);
fstream InSun("d:\\Arsenal2019New\\Data\\KFSUN.txt", ios::in);

Vect e; Matr E; Str eS;

//Глобальные переменные:
double pi=3.141592653589793e0,	//pi
pi_4=0.78539816339749e0,		//pi/4
pi_2=1.57079632679489e0,		//pi/2
pi2=6.28318530717959e0,			//2*pi
rdn=57.29577951308232e0;		//180/pi

double rE_ = 6.371e6,			//Средний радиус Земли,[м]
aE_ = 6.378136e6,				//Большая полуось общеземного эллипсоида,[м]
bE_ = 6.356768264e6,			//Малая полуось общеземного эллипсоида,[м]
acm_ = 3.35280373518430e-3,		//Сжатие Земли (1/298.25784)
mu_ = 3.986004418e14,			//Гравитационный параметр Земли,[м*м*м/(с*с)]
muM_ = 4.902823e12,				//Гравитационный параметр Луны,[м*м*м/(с*с)]
muS_ = 1.32712517e20,			//Гравитационный параметр Солнца,[м*м*м/(с*с)]
omE_ = 7.292115085e-5,			//Угловая скорость вращения Земли,[рад/с]
Apr_ = 4.15196e11,				//Коэф. прецессии долготы восх. узла [м*м]
vC_ = 299792458,				//Скорость света [м/с]
ro0_ = 1.58868e-8;				//Плотность ночной атмосферы
								//на высоте 120 км [кг/(м*м*м)]

//Инициализация статических компонент классов:
int Matr::numbMatr=0, Matr::numbMatrMx=0, Matr::NumbMatr = 0, Matr::NumbVect = 0;
int Vect::numbVect=0, Vect::numbVectMx=0, Vect::NumbVect = 0, Vect::NumbMatr = 0;
int Tensor::numbTensor = 0, Tensor::numbTensorMx = 1;
int Str::numbStr = 0, Str::numbStrMx = 0, Str::NumbStr = 0;
//int VStr::numbVStr = 0, VStr::numbVStrMx = 0;
//int MStr::numbMStr = 0, MStr::numbMStrMx = 0;
//int TStr::numbTStr = 0, TStr::numbTStrMx = 0;

//					К о м п о н е н т н ы е   ф у н к ц и и  (м е т о д ы   к л а с с о в)

Matr::Matr(int p, int q) //M2. Конструктор класса Matr
{
	m = p; n = q;
	A = Array(m, n);
	numbMatr++;
	if (numbMatr > numbMatrMx)numbMatrMx = numbMatr;
}

Matr::~Matr() //M3. Деструктор класса Matr
{
	int i;
	if (numbMatr > 0)
	{
		DO(i, 0, m-1) delete A[i];
		delete[] A;
        numbMatr--;
	}
	if (NumbMatr > 0)
	{
		delete H;
		NumbMatr--;
	}
	if (NumbVect > 0)
	{
		delete h;
		NumbVect--;
	}
}

double& Matr::operator()(int i, int j) //M5. Определение (i,j)-го элемента матрицы
{
	if (i<1 || i>m || j<1 || j>n)
	{
		Ecrir(cout, "\n\n Matr::operator()(int i,int j): ");
		Ecrir(cout, "specifieded indexes overrun!");
		Ecrir(cout, "\n -----> i=", i, "   j=", j, "   m=", m, "   n=", n, "\n");
		END
	}
	return A[i - 1][j - 1];
}

Matr& Matr::operator()(int i1, int i2, int j1, int j2) //M6. Сечение матрицы
{
	if (i2<i1 || j2<j1 || i2>m || j2>n)
	{
		Ecrir(cout, "\n\n Matr::operator()(int i1,int i2,int j1,int j2): ");
		Ecrir(cout, "indexes is specified wrong!\n");
		Ecrir(cout, " -----> i1=", i1, "   i2=", i2);
		Ecrir(cout, "   j1=", j1, "   j2=", j2, "   m=", m, "   n=", n);
		END
	}
	int i, j; 
	STACKM(i2 - i1 + 1, j2 - j1 + 1)
	DO(i, i1, i2)DO(j, j1, j2)(*H)(i - i1 + 1, j - j1 + 1) = A[i - 1][j - 1];
	return *H;
}

Vect& Matr::operator()(int i1, int i2, int i) //M7. Сечение i-й строки матрицы
{
	int j;
	STACKV(i2 - i1 + 1)
	DO(j, i1, i2)(*h)(j-i1+1) = A[i - 1][j - 1];
	return *h;
}

Matr& Matr::operator()(int k, int n0, Vect& a) //M8. Копирование вектора в сечение k-й строки матрицы
{
	int i, nmx = min(a->m, n - n0 + 1);
	DO(i, n0, n0 + nmx - 1)A[k - 1][i - 1] = a(i - n0 + 1);
	return *this;
}

Vect& Matr::operator()(int Ns) //M9. Выделение из матрицы строки под номером Ns
{
	if (Ns<0 || Ns>m)
	{
		Ecrir(cout, "\n\n Matr::operator()(int Ns): error in index!\n");
		Ecrir(cout, " Ns=", Ns, "  m=", m);
		END
	}
	int i;
	STACKV(n)
	DO(i, 1, n)(*h)(i) = A[Ns - 1][i - 1];
	return *h;
}

Vect& Matr::operator[](int Nc) //M10. Выделение из матрицы столбца под номером Nc
{
	if (Nc<0 || Nc>n)
	{
		Ecrir(cout, "\n\n Matr::operator[]: specified a number of line overruns!\n");
		Ecrir(cout, " Nc=", Nc, "  n=", n);
		END
	}
	int i;
	STACKV(m)
	DO(i, 1, m)(*h)(i) = A[i - 1][Nc - 1];
	return *h;
}

Matr& Matr::operator()(Vect& a, int m0, int k) //M11. Копирование вектора в сечение k-го столбца матрицы
{
	int i, mmx = min(a->m, m - m0 + 1);
	DO(i, m0, m0 + mmx - 1)A[i - 1][k - 1] = a(i - m0 + 1);
	return *this;
}

Matr& Matr::operator()(Vect& s) //M12. Выделение строк матрицы, на которые указывает вектор s
{
	int ms = s->m, i, j, k;
	STACKM(ms, n)
	DO(i, 1, ms)DO(j, 1, n) 
	{ 
		k = int(s(i));
		(*H)(i, j) = A[k - 1][j - 1]; 
	} 
	return *H;
}

Matr& Matr::operator[](Vect& c) //M13. Выделение столбцов матрицы, на которые указывает вектор с
{
	int nc = c.m, i, j, k;
	STACKM(m, nc)
	DO(i, 1, m)DO(j, 1, nc) 
	{ 
		k = int(c(j)); 
		(*H)(i, j) = A[i - 1][k - 1]; 
	} 
	return *H;
}

Matr& Matr::operator()(Vect& s, Vect& c) //M14. Выделение подматрицы из матрицы
{										 //по векторам-указателям s,c
	int ms = s->m, nc = c->m, i, j;
	Matr B(m, n), C(ms, n);
	STACKM(ms, nc)
	DO(i, 1, m)DO(j, 1, n)B(i, j) = A[i - 1][j - 1];
	C = B(s);
	*H = C[c];
	return *H;
}

//Размерность результата совпадает с размерностью первого операнда
Matr& Matr::operator-(Matr& B) //M15. Разность матриц
{
	int i, j, mMin = min(m, B->m), nMin = min(n, B->n);
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)
		if (i <= mMin && j <= nMin)(*H)(i, j) = A[i - 1][j - 1] - B(i, j);
		else (*H)(i, j) = A[i - 1][j - 1];
	return *H;
}

//Размерность результата совпадает с размерностью первого операнда
Matr& Matr::operator+(Matr& B) //M16. Сумма матриц
{
	int i, j, mMin = min(m, B->m), nMin = min(n, B->n);
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)
		if (i <= mMin && j <= nMin)(*H)(i, j) = A[i - 1][j - 1] + B(i, j);
		else (*H)(i, j) = A[i - 1][j - 1];
	return *H;
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

Vect& Matr::operator*(Vect& b) //M18. Перемножение матрицы на вектор
{
	if (n != b->m)
	{
		Ecrir(cout, "\n\n Matr::operator*(Vect& b): dimensionality of ");
		Ecrir(cout, "multiplicands are not coordinated!\n");
		Ecrir(cout, " -----> A->n=", n, "   b->m=", b->m); END
	}
	int i, j;
	STACKV(m)
	DO(i, 1, m)DO(j, 1, n)(*h)(i) += A[i - 1][j - 1] * b(j);
	return *h;
}

Matr& Matr::operator-(double q) //M19. Вычитание скаляра из матрицы
{
	int i, j;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = A[i - 1][j - 1] - q;
	return *H;
}

Matr& Matr::operator+(double q) //M20. Сумма матрицы со скаляром
{
	int i, j;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = A[i - 1][j - 1] + q;
	return *H;
}

Matr& Matr::operator*(double q) //M21. Умножение матрицы на скаляр
{
	int i, j;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = A[i - 1][j - 1] * q;
	return *H;
}

Matr& Matr::operator/(double q) //M22. Деление матрицы на скаляр
{
	if (q == 0)
	{
		Ecrir(cout,"\n	Error in Matr::operator/(double q), q=",q,"\n");
		END
	}
	int i, j;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = A[i - 1][j - 1] / q;
	return *H;
}

Matr& Matr::operator/(Matr& B) //M23. Деление матрицы на матрицу
{
	int i, j;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)
		if (B(i, j) != 0)(*H)(i, j) = A[i - 1][j - 1] / B(i, j);
	return *H;
}

Matr& Matr::operator!() //M24. Транспонирование матрицы
{
	int i, j;
	STACKM(n, m)
	DO(i, 1, n)DO(j, 1, m)(*H)(i, j) = A[j - 1][i - 1];
	return *H;
}

Matr& Matr::operator-() //M25. Унарный минус
{
	int i, j;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = -A[i - 1][j - 1];
	return *H;
}

Matr& Matr::operator~() //M26. Обращение матрицы
{
	int F, G, O, P, Q, R, T, U, BB, BD, j, i, m, k, l;
	double W, BA, BJ, BN, BQ, BE, X, Y, Z; Matr S(n, 2);
	Matr C(n,n);
	STACKM(n, n)
	DO(i, 1, n)DO(j, 1, n)C(i, j) = A[i - 1][j - 1];
	if (n == 1)goto L18; F = 1; G = 1;
L1: U = G; W = 0;
	DO(j, G, n)
	{
		BB = F; BA = 0;
		DO(i, F, n)
		{
			BJ = fabs(C(i, j));
			if (BA >= BJ)continue;
			BA = BJ; BB = i;
		}
		if (W > BA)continue;
		W = BA; T = BB; U = j;
	}
	DO(i, 1, n) { BN = C(i, U); C(i, U) = C(i, F); C(i, F) = BN; }
	DO(j, 1, n) { BQ = C(T, j); C(T, j) = C(F, j); C(F, j) = BQ; }
	S(F, 2) = U; S(F, 1) = T; j = G; X = 1 / C(F, F);
	DO(i, 1, n)C(i, F) = -(C(i, F) * X);
	C(F, F) = X; O = 0;
L7: O++;
	if (O == F)goto L7;
	if (n < O)goto L10;
	BE = C(F, O); BD = F - 1;
	if (BD == 0)goto L30;
	DO(P, 1, BD)C(P, O) = C(P, O) + C(P, F) * BE;
L30:
	BD = F + 1;
	if (BD > n)goto L7;
	DO(P, BD, n)C(P, O) = C(P, O) + C(P, F) * BE;
	goto L7;
L10:
	m = j + 1; if (m > n)goto L32;
	DO(O, m, n)C(F, O) = C(F, O) * C(F, G);
L32:
	m = j - 1; if (m == 0)goto L31;
	DO(O, 1, m)C(F, O) = C(F, O) * C(F, G);
L31:
	F++; G++; if (F > n)goto L13; if (G <= n)goto L1;
L13:
	j = 1;
L14:
	Q = sign(1 - n); k = n - Q;
L15:
	k += int(Q); i = int(S(k, 2)); Y = C(i, j); C(i, j) = C(k, j); C(k, j) = Y;
	if (k != 1)goto L15; j++; if (j <= n)goto L14; i = 1;
L16:
	R = sign(1 - n); l = n - R;
L17:
	l = l + R; j = int(S(l, 1)); Z = C(i, j); C(i, j) = C(i, l); C(i, l) = Z;
	if (l != 1)goto L17; i++; if (i <= n)goto L16; 
	*H = C;  return *H;
L18:
	C(1, 1) = 1 / C(1, 1);
	*H = C;  return *H;
}

Matr& Matr::operator=(double q) //M27. Присваивание скаляра матрицей
{
	int i, j;
	DO(i, 1, m)DO(j, 1, n)A[i - 1][j - 1] = q;
	return *this;
}

Matr& Matr::operator=(Matr& B) //M28. Присваивание матрицы матрицей
{
	int i, j, mB = B->m, nB = B->n, mm = min(m, mB), nn = min(n, nB);
	DO(i, 1, mm)DO(j, 1, nn)A[i - 1][j - 1] = B(i, j);
	return *this;
}

Matr& Matr::operator()(int k, Vect& a) //M29. Копирование вектора в k-ю строку матрицы
{
	int i, nmx = min(a->m, n);
	DO(i, 1, nmx)A[k - 1][i - 1] = a(i);
	return *this;
}

Matr& Matr::operator()(Vect& a, int k) //M30. Копирование вектора в k-й столбец матрицы
{
	int i, mmx = min(a->m, m);
	DO(i, 1, mmx)A[i - 1][k - 1] = a(i);
	return *this;
}

//M31. Копирование матрицы в сечение матрицы.
//B - копируемая матрица;
//m0,n0 - начало сечения матрицы-приемника 
Matr& Matr::operator()(Matr& B, int m0, int n0)
{
	int i, j, mmx = min(B->m + m0 - 1, m), nmx = min(B->n + n0 - 1, n);
	DO(i, m0, mmx)DO(j, n0, nmx)A[i - 1][j - 1] = B(i - m0 + 1, j - n0 + 1);
	return *this;
}

//M32. Тестирование на ноль. Возвращает true, если матрица нулевая
bool Matr::TestZero()
{
	int i, j;
	DO(i, 1, m)DO(j, 1, n)if (fabs(A[i - 1][j - 1])>1e-15)return false;
	return true;
}

//M33. Извлечение корня квадратного из симметричной матрицы P>0,
//															   T
//т.е. определение такой нижней треугольной матрицы C, что P=CC
Matr& Matr::RootMat()
{
	int i, j, k; double s;
	Matr C(n, n);
	STACKM(n, n)
	if (TestZero()) { *H = 0; return *H; }
	DO(i, 1, n)C(i, 1) = A[i - 1][0] / sqrt(A[0][0]);
	DO(i, 2, n)
		DO(j, 2, i)
			if (j != i)
			{
				s = 0;
				DO(k, 1, j - 1)s += C(i, k) * C(j, k);
				if (C(j, j) == 0)C(j, j) = 1e-50;
				C(i, j) = (A[i - 1][j - 1] - s) / C(j, j);
			}
			else
			{
				s = 0;
				DO(k, 1, i - 1)s += C(i, k) * C(i, k);
				C(i, i) = sqrt(A[i - 1][i - 1] - s);
			}
	*H = C;
	return *H;
}

Matr& Matr::Pow(double q) //M34. Поэлементное возведение матрицы в степень
{
	int i, j;
	Matr B(m, n);
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)
	{
		if (A[i - 1][j - 1] == 0 && q<0) 
		{ 
			B(i, j) = 0; 
			continue; 
		}
		B(i, j) = pow(A[i - 1][j - 1], q);
	}
	*H = B;
	return *H;
}

Matr& Matr::Dmat() //M35. Обнуление недиагональных элементов матрицы
{
	int i, j;
	DO(i, 1, m - 1)DO(j, i + 1, m) 
	{ 
		A[i - 1][j - 1] = 0; 
		A[j - 1][i - 1] = 0; 
	}
	return *this;
}

Matr& Matr::Dmat(Vect& a) //M36. Формирование диагональной матрицы по вектору
{
	int i;
	DO(i, 1, m)A[i - 1][i - 1] = a(i); 
	return *this;
}

Matr& Matr::Dmat(double q) //M37. Формирование диагональной матрицы по скаляру
{
	int i;
	DO(i, 1, m)A[i - 1][i - 1] = q; 
	return *this;
}

Matr& Matr::Grevill() //M38. Определение псевдообратной матрицы
{
	int i, j, k; double s;
	Matr Apin(n, m);
	STACKM(n, m)
	if (m >= n)
	{
		Vect b(m), c(m), a(m); DO(i, 1, m)a(i) = A[i - 1][0];
		s = a * a; DO(i, 1, m)Apin(1, i) = a(i) / s;
		DO(k, 2, n)
		{
			Vect d(k - 1); Matr Asec(m, k - 1), ApinSec(k - 1, m);
			DO(i, 1, k - 1)DO(j, 1, m)ApinSec(i, j) = Apin(i, j);
			DO(i, 1, m)DO(j, 1, k - 1)Asec(i, j) = A[i - 1][j - 1];
			DO(i, 1, m)a(i) = A[i - 1][k - 1];
			d = ApinSec * a; c = a - Asec * d; s = c * c;
			if (Sum(Abs(c)) != 0)b = c / s;
			else { b = d * ApinSec; s = 1 + d * d; b = b / s; }
			DO(i, 1, k - 1)DO(j, 1, m)Apin(i, j) = ApinSec(i, j) - d(i) * b(j);
			DO(i, 1, m)Apin(k, i) = b(i);
		}
	}
	else
	{
		Vect b(n), c(n), a(n); DO(i, 1, n)a(i) = A[0][i - 1];
		s = a * a; DO(i, 1, n)Apin(i, 1) = a(i) / s;
		DO(k, 2, m)
		{
			Vect d(k - 1); Matr Asec(k - 1, n), ApinSec(n, k - 1);
			DO(i, 1, n)DO(j, 1, k - 1)ApinSec(i, j) = Apin(i, j);
			DO(i, 1, k - 1)DO(j, 1, n)Asec(i, j) = A[i - 1][j - 1];
			DO(i, 1, n)a(i) = A[k - 1][i - 1];
			d = a * ApinSec; c = a - d * Asec; s = c * c;
			if (Sum(Abs(c)) != 0)b = c / s;
			else { b = ApinSec * d; s = 1 + d * d; b = b / s; }
			DO(i, 1, n)DO(j, 1, k - 1)
				Apin(i, j) = ApinSec(i, j) - b(i)*d(j);
			DO(i, 1, n)Apin(i, k) = b(i);
		}
	}
	*H = Apin;
	return *H;
}

//M39. Перестановка двух строк матрицы в соответствии
//с указателями номеров переставляемых строк i,j;
//TrStr(i,j)=TrStr(j,i)
Matr& Matr::TrStr(int i, int j)
{
	if (i <= 0 || j <= 0 || i > m || j > m)
	{
	 	Ecrir(cout, "\n\n Matr::TrStr: Error in index!");
		Ecrir(cout, "\n m=", m, " n=", n, " i=", i, " j=", j, "\n"); END
	}
	if (i == j)return *this;
	int k, l;
	Matr B(m, n); Vect b(n);
	DO(k, 1, m)DO(l, 1, n)B(k, l) = A[k - 1][l - 1];
	b = B(i); B(i, B(j)); B(j, b);
	DO(k, 1, m)DO(l, 1, n)A[k - 1][l - 1] = B(k, l);
	return *this; 
}

//M40. Перестановка двух столбцов матрицы в соответствии
//с указателями номеров переставляемых столбцов i,j;
//TrCol(i,j)=TrCol(j,i)
Matr& Matr::TrCol(int i, int j)
{
	if (i <= 0 || j <= 0 || i > n || j > n)
	{
		Ecrir(cout, "\n\n Matr::TrCol: Error in index!");
		Ecrir(cout, "\n m=", m, " n=", n, " i=", i, " j=", j, "\n");
		END
	}
	if (i == j)return *this;
	int k, l;
	Matr B(m, n); Vect b(m);
	DO(k, 1, m)DO(l, 1, n)B(k, l) = A[k - 1][l - 1];
	b = B[i]; B(B[j], i); B(b, j);
	DO(k, 1, m)DO(l, 1, n)A[k - 1][l - 1] = B(k, l);
	return *this;
}

//M41. Функция SubMatr преобразует матрицу полного ранга
//A(m, n) (m > n). В результате преобразования первые
//m строк образуют базис
Matr& Matr::SubMatr()
{
	int i, j, rang, k = 1, si;
	Vect s(m), c(n); 
	Matr B(m, n);
	DO(i, 1, m)DO(j, 1, n)B(i, j) = A[i - 1][j - 1];
	rang = RangMatr(B, s, c, 1e-12);
	DO(i, 1, rang)
	{ 
		si = int(s(k));
		if (i != si) B->TrStr(i, si);
		else { k++; continue; }
	}
	DO(i, 1, m)DO(j, 1, n)A[i - 1][j - 1] = B(i, j);
	return *this;
}

//M42. Функция SubMatr преобразует матрицу полного ранга
//A(m, n) (m > n). В результате преобразования первые
//m строк образуют базис.
//s(m)	- указатель базисных строк (вход)
Matr& Matr::SubMatr(int rang, Vect& s)
{
	int i, j, k = 1, si;
	Matr B(m, n);
	DO(i, 1, m)DO(j, 1, n)B(i, j) = A[i - 1][j - 1];
	DO(i, 1, rang)
	{
		si = int(s(k));
		if (i != si) B->TrStr(i, si);
		else { k++; continue; }
	}
	DO(i, 1, m)DO(j, 1, n)A[i - 1][j - 1] = B(i, j);
	return *this;
}

//M43. Возвращает true, если обнаруживает отрицательный диагональный элемент
bool Matr::Negativ()
{
	int i; 
	DO(i, 1, m)if (A[i - 1][i - 1]<0)return true; 
	return false;
}

Vect& Matr::MatVecStr() //M44. Упаковка матрицы в вектор по строкам
{
	int i, j;
	STACKV(m*n)
	DO(i, 1, m)DO(j, 1, n)(*h)((i - 1)*n + j) = A[i - 1][j - 1];
	return *h;
}

Vect& Matr::MatVecCol() //M45. Упаковка матрицы в вектор по столбцам
{
	int i, j; 
	STACKV(m*n)
	DO(i, 1, n)DO(j, 1, m)(*h)((i - 1)*m + j) = A[j - 1][i - 1];
	return *h;
}

//M46. Упаковка вектора, размерности mn в матрицу
//размерности (m,n) по строкам
Matr& Matr::VecMatStr(Vect& a) 
{	
	int i, j; 
	DO(i, 1, m)DO(j, 1, n)A[i - 1][j - 1] = a((i - 1)*n + j); 
	return *this;
}

//M47. Упаковка вектора, размерности mn в матрицу
//размерности (m,n) по столбцам
Matr& Matr::VecMatCol(Vect& a)
{
	int i, j;
	DO(i, 1, n)DO(j, 1, m)A[j - 1][i - 1] = a((i - 1)*m + j); 
	return *this;
}

//M48. Критерий Сильвестра
//Возвращает 1, если матрица положительно определенная;
//возвращает -1, если матрица отрицательно определенная;
//возвращает 0, если квадратичная форма знакопеременная
//или полупределённая
int Matr::Silvestr() 
{
	int i, j, k; double d0, d; Matr P(m, m);
	DO(i, 1, m)DO(j, 1, m)P(i, j) = A[i - 1][j - 1];
    static_cast<void>(d0 = P(1, 1)), d;
	if (d0 > 0)
	{
		DO(k, 2, m)
		{
			P->m = k; P->n = k; //Это вполне корректно (Денис Звягинцев)
			d = det(P);
			if (d0 * d > 0)d0 = d;
			else return 0;
		}
		return 1;
	}
	if (d0 < 0)
	{
		DO(k, 2, m)
		{
			P->m = k; P->n = k;
			d = det(P);
			if (d0 * d < 0)d0 = d;
			else return 0;
		}
		return -1;
	}
	return 0;
}

Matr& Matr::Cut(int I, int J) //M49. Удаляет из матрицы I-ю строку и J-й столбец
{
	int i, j, k, l, m1 = m - 1, n1 = n - 1; MSTACK(m1, n1, P)
		DO(i, 1, m1)
	{
		if (i<I)k = i - 1; else k = i;
		DO(j, 1, n1)
		{
			if (j<J)l = j - 1; else l = j;
			P(i, j) = A[k][l];
		}
	}
	return P;
}

double Matr::Max() //M50. Определение значения максимального элемента матрицы
{
	int i, j; double a = A[0][0];
	DO(i, 0, m - 1)DO(j, 0, n - 1)if (A[i][j]>a)a = A[i][j]; 
	return a;
}

double Matr::Min() //M51. Определение значения минимального элемента матрицы
{
	int i, j; double a = A[0][0];
	DO(i, 0, m - 1)DO(j, 0, n - 1)if (A[i][j]<a)a = A[i][j]; 
	return a;
}

//M52. Функция, заменяющая основной конструктор класса Matr,
//для объекта, созданного конструктором без параметров 
void Matr::M(int p, int q, Matr& X)
{
	int i, j; m = p; n = q; 
	A = Array(m, n);
	DO(i, 1, m)DO(j, 1, n)X(i, j) = A[i - 1][j - 1];
	numbMatr++; if (numbMatr>numbMatrMx)numbMatrMx = numbMatr;
}

//M53. Функция, заменяющая основной конструктор класса Matr,
//для объекта, созданного конструктором без параметров 
Matr& Matr::M(int p, int q)
{
	int i, j; m = p; n = q;
	STACKM(m, n) A = Array(m, n);
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = A[i - 1][j - 1];
	numbMatr++; if (numbMatr>numbMatrMx)numbMatrMx = numbMatr;
	return *H;
}

Matr& Matr::MatrixI(int n) //M54. Единичная матрица
{
	int i;
	STACKM(n, n)
		DO(i, 1, n) (*H)(i, i) = 1;
	return *H;
}

//MG1. Вычитание матрицы A из скаляра q.
//Глобальный оператор для класса Matr
Matr& operator-(double q, Matr& A)
{
	int i, j, m = A->m, n = A->n;
	A->H = new Matr(m, n); A->NumbMatr++;
	DO(i, 1, m)DO(j, 1, n)(*(A->H))(i, j) = q - A(i, j);
	return *(A->H);
}

//MG2. Сложение скаляра q с матрицей A.
//Глобальный оператор для класса Matr
Matr& operator+(double q, Matr& A)
{
	int i, j, m = A->m, n = A->n;
	A->H = new Matr(m, n); A->NumbMatr++;
	DO(i, 1, m)DO(j, 1, n)(*(A->H))(i, j) = q + A(i, j);
	return *(A->H);
}

//MG3. Умножение скаляра q с матрицей A.
//Глобальный оператор для класса Matr
Matr& operator*(double q, Matr& A)
{
	int i, j, m = A->m, n = A->n;
	A->H = new Matr(m, n); A->NumbMatr++;
	DO(i, 1, m)DO(j, 1, n)(*(A->H))(i, j) = q * A(i, j);
	return *(A->H);
}

//MG4. Деление скаляра q на матрицу A.
//Глобальный оператор для класса Matr
Matr& operator/(double q, Matr& A)
{
	int i, j, m = A->m, n = A->n;
	A->H = new Matr(m, n); A->NumbMatr++;
	DO(i, 1, m)DO(j, 1, n)(*(A->H))(i, j) = q / A(i, j);
	return *(A->H);
}

//MG5. Прямое произведение матриц.
//Глобальный оператор для класса Matr
Matr& operator&(Matr& A, Matr& B)
{
	int i, j, m = A->m, n = A->n;
	A->H = new Matr(m, n); A->NumbMatr++;
	DO(i, 1, m)DO(j, 1, n)(*(A->H))(i, j) = A(i, j) * B(i, j);
	return *(A->H);
}

Vect::Vect(int p) //V2. Конструктор класса Vect
{
	m = p;
	a=Array(m); numbVect++;
	if(numbVect>numbVectMx)numbVectMx=numbVect;
	//Ecrir(cout, "\n Constructor Vect(): numbVect=", numbVect,
	//	" numbVectMx=", numbVectMx, " NumbVect=", NumbVect, " NumbMatr=", NumbMatr);
}

Vect::~Vect() //V3. Деструктор класса Vect
{
	if (numbVect > 0)
	{
		delete[] a;
		numbVect--;
	}
	/*if (NumbVect > 0)
	{
		delete h;
		NumbVect--;
	}
	if (NumbMatr > 0)
	{
		delete H;
		NumbMatr--;
	}*/
	//Ecrir(cout, "\n Destructor ~Vect(): numbVect=", numbVect,
	//	" numbVectMx=", numbVectMx, " NumbVect=", NumbVect, " NumbMatr=", NumbMatr);
}

//V4. Определение значения i-й координаты вектора.
//i-я координата может стоять,
//как справа от знака =, так и слева,
//т.е. может изменяться.
//Поэтому возвращаемое значение имеет тип double&,
//а не просто double 
double& Vect::operator()(int i)
{
	if (i > 0 && i <= m)return a[i - 1];
	else
	{
		Ecrir(cout, "\n\n Vect::operator(): error in index!");
		Ecrir(cout, "\n m=", m, "  i=", i,
			  "  numbVect=", numbVect, "  numbVectMx=", numbVectMx);
		END
	}
	return a[i - 1];
}

double Vect::operator*(Vect& x) //V5. Скалярное произведение векторов
{
	int i, m_ = min(m, x->m); double s = 0;
	DO(i, 1, m_)s += a[i - 1] * x(i);
	return s;
}

double Vect::operator!() //V6. Модуль вектора
{
	int i; double s = 0;
	DO(i, 1, m)s += a[i - 1] * a[i - 1];
	return sqrt(s);
}

/*//V7. Сложение двух векторов.
//Размерность результата совпадает с размерностью первого (левого) операнда
Vect& Vect::operator+(Vect& x)
{
	int i, m_ = min(m, x->m);
	STACKV(m)
	DO(i, 1, m)
	if (i <= m_)(*h)(i) = a[i - 1] + x(i);
	else (*h)(i) = a[i - 1];
	return *h;
}*/

//V7. Сложение двух векторов.
//Размерность результата совпадает с размерностью первого (левого) операнда
Vect& Vect::operator+(Vect& x)
{
	int i, m_ = min(m, x->m);
	VSTACK(m, y)
	DO(i, 1, m)
	if (i <= m_)y(i) = a[i - 1] + x(i);
	else y(i) = a[i - 1];
	return y;
}

//V8. Вычитание двух векторов.
//Размерность результата совпадает с размерностью первого (левого) операнда
Vect& Vect::operator-(Vect& x)
{
	int i, m_= min(m, x->m);
	STACKV(m)
	DO(i, 1, m)
	if (i <= m_)(*h)(i) = a[i - 1] - x(i);
	else (*h)(i) = a[i - 1];
	return *h;
}

Vect& Vect::operator=(double q) //V9. Присваивание скаляра вектору
{
	int i;
	DO(i, 1, m) a[i - 1] = q;
	return *this;
}

Vect& Vect::operator=(Vect& x) //V10. Присваивание вектора вектору
{
	int i, m_ = min(m, x->m);
	DO(i, 1, m_) a[i - 1] = x(i);
	return *this;
}

Vect& Vect::operator*(double q) //V11. Умножение вектора на скаляр
{
	int i;
	STACKV(m)
	DO(i, 1, m)(*h)(i) = a[i - 1] * q;
	return *h;
}

Vect& Vect::operator/(double q) //V12. Деление вектора на скаляр
{
	int i;
	STACKV(m)
	DO(i, 1, m)
		if(q!=0)(*h)(i) = a[i - 1] / q;
		else
		{
			Ecrir(cout,"\n		Error! Vect::operator/(double q):  q=",q,"\n");
			END
		}
	return *h;
}

Vect& Vect::operator()(int i1, int i2) //V14. Сечение вектора
{
	int i;
	STACKV(i2 - i1 + 1)
		if (i1 > 0 && i1 <= m && i2 > 0 && i2 <= m && i1 <= i2)
			DO(i, i1, i2)(*h)(i - i1 + 1) = a[i - 1];
		else
		{
			cout << "\n\n Vect::operator()(int i1,int i2): error in index!";
			Ecrir(cout, "\n m=", m, "  i1=", i1, "  i2=", i2,
				"  numbVect=", numbVect, "  numbVectMx=", numbVectMx);
			END
		}
	return *h;
}

Vect& Vect::operator()(Vect& pointer) //V15. Выделение компонент, на которые указывает pointer
{
	int mp = pointer->m, i, k;
	STACKV(mp)
		DO(i, 1, mp)
	{
		k = pointer(i);
		if (k > 0 && k <= m)(*h)(i) = a[k - 1];
	}
	return *h;
}

Vect& Vect::operator*(Matr& A) //V16. Умножение вектор-строки на матрицу
{
	if (m != A->m)
	{
		Ecrir(cout, "\n\n Vect::operator*(Matr& A): dimensionality of ");
		Ecrir(cout, "multiplicands are not coordinated!\n");
		Ecrir(cout, " -----> a->m=", m, "   A->m=", A->m); END
	}
	int i, j, m = A->m, n = A->n;
	STACKV(n)
	DO(i, 1, n)
	DO(j, 1, m) (*h)(i) = (*h)(i) + a[j - 1] * A(j, i);
	return *h;
}

Matr& Vect::operator^(Vect& x) //V17. Диадное произведение
{
	int i, j, n = x->m;
	STACKM(m, n)
	DO(i, 1, m)DO(j, 1, n)(*H)(i, j) = a[i - 1] * x(j);
	return *H;
}

double Vect::operator%(Vect& b) //V18. Угол между векторами
{
	int i; double aa = 0, bb = 0, ab = 0;
	DO(i, 1, m)
	{ 
		aa += a[i - 1] * a[i - 1]; bb += b(i) * b(i); 
		ab += a[i - 1] * b(i); 
	}
	return acos(ab / sqrt(aa*bb));
}

Vect& Vect::operator-(double q) //V19. Вычитание скаляра из вектора
{
	int i;
	STACKV(m)
	DO(i, 1, m)(*h)(i) = a[i - 1] - q;
	return *h;
}

/*Vect& Vect::operator+(double q) //V20. Сложение скаляра с вектором
{
	int i;
	STACKV(m)
	DO(i, 1, m)(*h)(i) = a[i - 1] + q;
	return *h;
}*/

Vect& Vect::operator+(double q) //V20. Сложение скаляра с вектором
{
	int i;
	VSTACK(m, x)
	DO(i, 1, m)x(i) = a[i - 1] + q;
	return x;
}

Vect& Vect::operator/(Vect& b) //V21. Деление вектора на вектор
{
	int i;
	STACKV(m)
	DO(i, 1, m) if (b(i) != 0)(*h)(i) = a[i - 1] / b(i);
	return *h;
}

Vect& Vect::operator()(Vect& x, int m0) //V22. Копирование вектора в сечение вектора

{
	int i, mmx = min(x->m + m0 - 1, m);
	DO(i, m0, mmx)a[i - 1] = x(i - m0 + 1);
	return *this;
}

//От V22 отличается только порядком следования аргументов
Vect& Vect::operator()(int m0, Vect& x) //V23. Копирование вектора в сечение вектора
{
	int i, mmx = min(x->m + m0 - 1, m);
	DO(i, m0, mmx)a[i - 1] = x(i - m0 + 1);
	return *this;
}

Vect& Vect::operator()(int m0, double x) //V24. Инициализация сечения вектора (от m0 до m) числом
{
	int i;
	DO(i, m0, m)a[i - 1] = x;
	return *this;
}

Vect& Vect::operator()(int m1, int m2, double x) //V25. Инициализация сечения вектора числом
{
	int i;
	DO(i, m1, m2)a[i - 1] = x;
	return *this;
}

Vect& Vect::operator()(double x, int m1, int m2) //V26. Инициализация сечения вектора числом
{
	int i;
	DO(i, m1, m2)a[i - 1] = x;
	return *this;
}

Vect& Vect::operator-() //V27. Унарный минус
{
	int i;
	STACKV(m)
	DO(i, 1, m)(*h)(i) = -a[i - 1];
	return *h;
}

Vect& Vect::Pow(double q) //V28. Возведение в степень
{
	int i;
	STACKV(m)
	DO(i, 1, m)(*h)(i) = pow(a[i - 1], q);
	return *h;
}

Vect& Vect::Dvec(Matr& A) //V29. Главная диагональ квадратной матрицы
{
	int i,m = A->m; 
	STACKV(m)
	DO(i, 1, m)(*h)(i) = A(i, i);
	return *h;
}

//Возвращает true, если вектор нулевой
bool Vect::TestZero() //V30. Тестирование на ноль
{ 
	int i; 
	DO(i, 1, m)if (a[i - 1] != 0)return false; 
	return true; 
}

Vect& Vect::TrStr(int i, int j) //V31. Перестановка местами i-ой и j-ой компонент 
{
	double b;
	b = a[i - 1]; a[i - 1] = a[j - 1];
	a[j - 1] = b;
	return *this;
}

double Vect::Min() //V32. Определение значения минимальной компоненты вектора
{
	int i; double b = a[0];
	DO(i, 1, m - 1) if(a[i]<b)b = a[i];
	return b;
}

double Vect::Max() //V33. Определение значения максимальной компоненты вектора
{
	int i; double b = a[0];
	DO(i, 1, m - 1) if(a[i]>b)b = a[i];
	return b;
}

void Vect::V(int p, Vect& x) //V34. Замена основного конструктора
{
	int i; m = p;
	a = Array(m);
	DO(i, 1, m)x(i) = a[i - 1];
	numbVect++;
	if (numbVect > numbVectMx)numbVectMx = numbVect;
}

Vect& Vect::V(int p) //V35. Замена основного конструктора
{
	int i; m = p;
	VSTACK(m, x)
	a = Array(m); DO(i, 1, m)x(i) = a[i - 1];
	numbVect++; if (numbVect > numbVectMx)numbVectMx = numbVect;
	return x;
}

Vect& Vect::Vector_e(int n, int k) //V36. Единичный вектор
{
	STACKV(n)
		* h = 0; (*h)(k) = 1;
	return *h;
}

//VG1. Прямое произведение векторов.
//Глобальный оператор для класса Vect
Vect& operator&(Vect& x, Vect& y)
{
	int i, m = x->m;
	//x->h = new Vect(m); x->NumbVect++;
	STACKh(m, x) //Память отводится не x, а x->h  
	DO(i, 1, m)(*(x->h))(i) = x(i) * y(i);
	return *(x->h);
}

//VG2. Умножение скаляра на вектор.
//Глобальный оператор для класса Vect
Vect& operator*(double q, Vect& x)
{
	int i, m=x->m;
	x->h = new Vect(m); x->NumbVect++;
	DO(i, 1, m)(*(x->h))(i) = q * x(i);
	return *(x->h);
}

//VG3. Вычитание из скаляра вектора.
//Глобальный оператор для класса Vect
Vect& operator-(double q, Vect& x)
{
	int i, m= x->m;
	x->h = new Vect(m); x->NumbVect++;
	DO(i, 1, m) (*(x->h))(i) = q - x(i);
	return *(x->h);
}

//VG4. Сложение скаляра с вектором.
//Глобальный оператор для класса Vect
Vect& operator+(double q, Vect& x)
{
	int i, m= x->m;
	x->h = new Vect(m); x->NumbVect++;
	DO(i, 1, m)(*(x->h))(i) = q + x(i);
	return *(x->h);
}

//VG5. Деление скаляра на вектор.
//Глобальный оператор для класса Vect
Vect& operator/(double q, Vect& x)
{
	int i, m= x->m;
	x->h = new Vect(m); x->NumbVect++;
	DO(i, 1, m)(*(x->h))(i) = q / x(i);
	return *(x->h);
}

//Конструктор класса Tensor
Tensor::Tensor(int n1, int n2, int n3)
{
	m1 = n1; m2 = n2; m3 = n3;
	T = Array(n1, n2, n3); numbTensor++;
	numbTensorMx++;
}

//Деструктор класса Tensor
Tensor::~Tensor()
{
	for (int i = 0; i < m1; i++)
		for (int j = 0; j < m2; j++)delete T[i][j];
	delete[]T;	numbTensor--;
}

double& Tensor::operator()(int k1, int k2, int k3)
{
	if (k1 > 0 && k1 <= m1 && k2 > 0 && k2 <= m2 && k3 > 0 && k3 <= m3)
		return T[k1 - 1][k2 - 1][k3 - 1];
	else { cout << "\n\nTensor(k1,k2,k3): error in index!"; END }
	return T[k1 - 1][k2 - 1][k3 - 1];
}

//Сечение тензора по второму и третьему измерению
//(множество элементов, для которых 2-й и 3-й индексы фиксированы) 
Vect& Tensor::operator()(char* c, int k2, int k3)
{
	VSTACK(m1, v)
	if (k2 > 0 && k2 <= m2 && k3 > 0 && k3 <= m3)
	{
		for (int i = 0; i < m1; i++)v(i + 1) = T[i][k2 - 1][k3 - 1];
		return v;
	}
	else { cout << "\n\n Tensor(char *c,int k2,int k3): error in index!"; END }
	return v;
}

//Копирование вектора в сечение по второму и третьему измерению 
Tensor& Tensor::operator()(Vect& x, int k2, int k3)
{
	if (k2 > 0 && k2 <= m2 && k3 > 0 && k3 <= m3)
	{
		for (int i = 0; i < m1; i++)T[i][k2 - 1][k3 - 1] = x(i + 1);
		return *this;
	}
	else { cout << "\n\n Tensor(Vect& x,int k2,int k3): error in index!"; END }
	return *this;
}

//Сечение тензора по первому и третьему измерению 
Vect& Tensor::operator()(int k1, char* c, int k3)
{
	VSTACK(m2, v)
	if (k1 > 0 && k1 <= m1 && k3 > 0 && k3 <= m3)
	{
		for (int i = 0; i < m2; i++)v(i + 1) = T[k1 - 1][i][k3 - 1];
		return v;
	}
	else { cout << "\n\n Tensor(int k1,char *c,int k3): error in index!"; END }
	return v;
}

//Копирование вектора в сечение по первому и третьему измерению 
Tensor& Tensor::operator()(int k1, Vect& x, int k3)
{
	if (k1 > 0 && k1 <= m1 && k3 > 0 && k3 <= m3)
	{
		for (int i = 0; i < m2; i++)T[k1 - 1][i][k3 - 1] = x(i + 1);
		return *this;
	}
	else { cout << "\n\n Tensor(int k1,Vect& x,int k3): error in index!"; END }
	return *this;
}

//Сечение тензора по первому и второму измерению 
Vect& Tensor::operator()(int k1, int k2, char* c)
{
	VSTACK(m3, v)
	if (k1 > 0 && k1 <= m1 && k2 > 0 && k2 <= m2)
	{
		for (int i = 0; i < m3; i++)v(i + 1) = T[k1 - 1][k2 - 1][i];
		return v;
	}
	else { cout << "\n\n Tensor(int k1,int k2,char *c): error in index!"; END }
	return v;
}

//Копирование вектора в сечение по первому и второму измерению 
Tensor& Tensor::operator()(int k1, int k2, Vect& x)
{
	if (k1 > 0 && k1 <= m1 && k2 > 0 && k2 <= m2)
	{
		for (int i = 0; i < m3; i++)T[k1 - 1][k2 - 1][i] = x(i + 1);
		return *this;
	}
	else { cout << "\n\n Tensor(int k1,int k2,Vect& x): error in index!"; END }
	return *this;
}

//Сечение тензора по первому измерению
Matr& Tensor::operator()(int k1, char* c2, char* c3)
{
	int i = 0, j = 0; MSTACK(m2, m3, M)
		if (k1 > 0 && k1 <= m1)
		{
			for (i = 0; i < m2; i++) for (j = 0; j < m3; j++)M(i + 1, j + 1) = T[k1 - 1][i][j];
			return M;
		}
		else { cout << "\n\n Tensor(int k1,char *c2,char *c3): error in index!"; END }
	return M;
}

//Копирование матрицы в сечение тензора по первому измерению
Tensor& Tensor::operator()(int k1, char* c2, char* c3, Matr& A)
{
	int i = 0, j = 0;
	if (k1 > 0 && k1 <= m1)
	{
		for (i = 0; i < m2; i++) for (j = 0; j < m3; j++)T[k1 - 1][i][j] = A(i + 1, j + 1);
		return *this;
	}
	else
	{
		cout << "\n\n Tensor(int k1,char *c2,char *c3,Matr& A): error in index!";
		END
	}
	return *this;
}

//Сечение тензора по второму измерению
Matr& Tensor::operator()(char* c1, int k2, char* c3)
{
	int i = 0, j = 0; MSTACK(m1, m3, M)
		if (k2 > 0 && k2 <= m2)
		{
			for (i = 0; i < m1; i++) for (j = 0; j < m3; j++)M(i + 1, j + 1) = T[i][k2 - 1][j];
			return M;
		}
		else { cout << "\n\n Tensor(char *c1,int k2,char *c3): error in index!"; END }
	return M;
}

//Копирование матрицы в сечение тензора по второму измерению
Tensor& Tensor::operator()(char* c1, int k2, char* c3, Matr& A)
{
	int i = 0, j = 0;
	if (k2 > 0 && k2 <= m2)
	{
		for (i = 0; i < m1; i++) for (j = 0; j < m3; j++)T[i][k2 - 1][j] = A(i + 1, j + 1);
		return *this;
	}
	else
	{
		cout << "\n\n Tensor(char *c1,int k2,char *c3,Matr& A): error in index!";
		END
	}
	return *this;
}

//Сечение тензора по третьему измерению
Matr& Tensor::operator()(char* c1, char* c2, int k3)
{
	int i = 0, j = 0; MSTACK(m1, m2, M)
		if (k3 > 0 && k3 <= m3)
		{
			for (i = 0; i < m1; i++) for (j = 0; j < m2; j++)
				M(i + 1, j + 1) = T[i][j][k3 - 1];
			return M;
		}
		else { cout << "\n\n Tensor(char *c1,char *c2,int k3): error in index!"; END }
	return M;
}

//Копирование матрицы в сечение тензора по третьему измерению
Tensor& Tensor::operator()(char* c1, char* c2, int k3, Matr& A)
{
	int i = 0, j = 0; MSTACK(m1, m2, M)
		if (k3 > 0 && k3 <= m3)
		{
			for (i = 0; i < m1; i++) for (j = 0; j < m2; j++)T[i][j][k3 - 1] = A(i + 1, j + 1);
			return *this;
		}
		else
		{
			cout << "\n\n Tensor(char *c1,char *c2,int k3,Matr& A): error in index!";
			END
		}
	return *this;
}

void Tensor::del()
{
	for (int i = 0; i < m1; i++)for (int j = 0; j < m2; j++)delete T[i][j];
	delete[]T;
}

//S2. Конструктор класса Str
Str::Str(int p)
{
	m = p; 
	s = Line(p); 
	numbStr++;
	if(numbStr > numbStrMx)
		numbStrMx = numbStr + 1;
}

//S3. Деструктор класса Str
Str::~Str() 
{ 
	if(numbStr > 0)
	{
		delete[]s;
		numbStr--;
	}
}

char& Str::operator()(int i) //S4. Доступ к i-му символу
{
	if (i > 0 && i <= m) return s[i - 1];
	else
	{
		cout << "Str: error in index!" << "\n";
		exit(1);
	}
	return s[i - 1];
}

//S5. Сечение строки
Str& Str::operator()(int i1, int i2)
{
	int i;
	STACKS(i2 - i1 + 1)
	if (i1 > 0 && i1 <= m && i2 > 0 && i2 <= m && i1 <= i2)
	{
		DO(i, i1, i2) (*h)(i - i1 + 1) = s[i - 1];
		return *h;
	}
	else
	{
		cout << "Str (Section): error in index!" << "\n";
		END
	}
}

//S6. Сцепление (конкатенация) строк
Str& Str::operator+(Str& x)
{
	int i = 0, n = m + x->m;
	STACKS(n)
	DO(i, 1, n)
		if(i <= m) (*h)(i) = s[i - 1]; 
		else (*h)(i) = x(i - m);
	return *h;
}

//S7. Сцепление (конкатенация) строки с переменной типа char
Str Str::operator+(char* x)
{
	int i = 0, n = m + strlen(x);
	STACKS(n)
	DO(i, 1, n)
		if (i <= m) (*h)(i) = s[i - 1]; 
		else (*h)(i) = x[i - m - 1];
	return *h;
}

//S8. Оператор присвоения строковой переменной
Str& Str::operator=(Str& x)
{
	int i = 0, n = x->m, mn = min(m, n);
	DO(i, 1, mn)s[i - 1] = x(i);
	return *this;
}

//S9. Оператор присвоения переменной типа char
Str& Str::operator=(char* x)
{
	int i = 0, n = strlen(x), mn = min(m, n);
	DO(i, 0, mn-1)s[i] = x[i];
	return *this;
}

//SG1. Сцепление (конкатенация) строки с переменной типа char.
//Глобальный оператор для класса Str
Str& operator+(char* x, Str& y)
{
	int i = 0, mx = strlen(x), my = y->m, n = mx + my;
	STACKSh(n, y)
	DO(i, 1, n)
		if (i <= mx) (*(y->h))(i) = x[i - 1]; 
		else (*(y->h))(i) = y(i - mx);
	return *(y->h);
}

//B1. Конструктор класса Ballist
Ballist::Ballist()
{
	Init();
	Aatm = new Matr(7, 7); Batm = new Matr(5, 7);
	Catm = new Matr(5, 7); Datm = new Matr(5, 7);
	Eatm = new Matr(9, 7); Fi1Atm = new Matr(1, 7);
	Latm = new Matr(5, 7); Natm = new Matr(3, 7);
}

//B2. Конструктор класса Ballist
Ballist::Ballist(int N)
{
	Ngrav = N; Init();
	Aatm = new Matr(7, 7); Batm = new Matr(5, 7);
	Catm = new Matr(5, 7); Datm = new Matr(5, 7);
	Eatm = new Matr(9, 7); Fi1Atm = new Matr(1, 7);
	Latm = new Matr(5, 7); Natm = new Matr(3, 7);
}

//B3. Деструктор класса Ballist
Ballist::~Ballist() 
{ 
	delete Aatm, Batm, Catm, Datm, Eatm, Fi1Atm, Latm, Natm; 
}

//B4. Инициализация компонент класса Ballist
void Ballist::Init()
{
	Start = 1; startMoon = 1; startSun = 1; startAtm = 1;
	kday = 0; OnMoon = 1; OnSun = 1; OnLux = 1; typeAtm = 2;
	gamma = 1.22; p = 4.65e-6;
}

//B5. Вычисление правых частей дифференциальных уравнений,            
//соответствующих модели движения орбитального объекта в J2000	
//(коэффициенты разложения геопотенциала - ненормированные).	      
//Входные параметры:
//t		- независимая переменная (время);                              
//x	    - фазовый вектор объекта в J2000.							   
//Возвращаемое значение:                                                    
//f		- вектор правых частей системы д.у. (производных)              
//Примечание: переменная t есть время в секундах равное сумме			   
//            текущего гринвичского времени и целого числа суток		
//            (в секундах), прошедших с начала интегрирования
Vect& Ballist::ModProgn(double t, Vect& x)
{
	int n = x->m; double tg, r3;
	Vect dg(3), aM(3), aS(3), dgAtm(3), dgLux(3), f(n);
	Matr CoefGrav(700, 2); static Vect cProgn(700), dProgn(700);
	if (Start) //Только при первом обращении к ModProgn
	{
		char s[50];
		ReadMatr(ent, s, CoefGrav); //Считывание коэффициентов разложения ГПЗ
								    //из потока ent, связанного с файлом
								    //...//Data//CoefGravPot.txt
		cProgn = CoefGrav[1]; dProgn = CoefGrav[2];
		ent.seekg(0); Start = 0;
	}
	if (t >= 86400 * (kday + 1.0))ud++; if (t < 86400.0 * kday)ud--; //Учет скачков юлианской даты
																     //из-за колебаний времени в
																     //алгоритмах интегрирования
																     //с переменным шагом
	kday = t / 86400; tg = t - kday * 86400.0;
	r3 = !x(1, 3); r3 = r3 * r3 * r3;
	f(1) = x(4); f(2) = x(5); f(3) = x(6);
	f(4) = -mu_ * x(1) / r3; f(5) = -mu_ * x(2) / r3; f(6) = -mu_ * x(3) / r3;
	if (Ngrav > 1)
	{
		GravJ2000(tg, cProgn, dProgn, x, dg);
		f(4) += dg(1); f(5) += dg(2); f(6) += dg(3);
	}
	if (OnMoon)
	{
		SpeedMoon(tg, x(1, 3), aM);
		f(4) += aM(1); f(5) += aM(2); f(6) += aM(3);
	}
	if (OnSun)
	{
		SpeedSun(tg, x(1, 3), aS);
		f(4) += aS(1); f(5) += aS(2); f(6) += aS(3);
	}
	if (typeAtm > 0)
	{
		SpeedAtm(tg, x, dgAtm);
		f(4) += dgAtm(1); f(5) += dgAtm(2); f(6) += dgAtm(3);
	}
	if (OnLux)
	{
		dgLux = PressLux(tg, x);
		f(4) += dgLux(1); f(5) += dgLux(2); f(6) += dgLux(3);
	}
	STACKh(n,x)
	*(x->h) = f;
	return *(x->h);
}

//B6. Вычисление правых частей дифференциальных уравнений,            
//соответствующих модели движения орбитального объекта в J2000	
//(коэффициенты разложения геопотенциала - ненормированные).	      
//Входные параметры:
//U(t,x)	- управляющее воздействие; 
//t			- независимая переменная (время);                              
//x			- фазовый вектор объекта в J2000.							   
//Возвращаемое значение:                                                    
//f			- вектор правых частей системы д.у. (производных)              
//Примечание: переменная t есть время в секундах равное сумме			   
//            текущего гринвичского времени и целого числа суток		
//            (в секундах), прошедших с начала интегрирования
Vect& Ballist::ModProgn(Vect& U(double t, Vect& x), double t, Vect& x)
{
	int n = x->m; double tg, r3;
	Vect u(3), dg(3), aM(3), aS(3), dgAtm(3), dgLux(3), f(n);
		Matr CoefGrav(700, 2); static Vect cProgn(700), dProgn(700);
	if (Start)
	{
		char s[50];
		ReadMatr(ent, s, CoefGrav);
		cProgn = CoefGrav[1]; dProgn = CoefGrav[2];
		ent.seekg(0); Start = 0;
	}
	if (t >= 86400 * (kday + 1))ud++; if (t < 86400 * kday)ud--; //Учет скачков юлианской даты
																 //из-за колебаний времени в
																 //алгоритмах интегрирования
																 //с переменным шагом
	kday = t / 86400;
	tg = t - kday * 86400.0;
	r3 = !x(1, 3); r3 = r3 * r3 * r3;
	f(1) = x(4); f(2) = x(5); f(3) = x(6);
	f(4) = -mu_ * x(1) / r3; f(5) = -mu_ * x(2) / r3; f(6) = -mu_ * x(3) / r3;
	u = U(t, x); f(4) += u(1); f(5) += u(2); f(6) += u(3);
	if (Ngrav > 1)
	{
		GravJ2000(tg, cProgn, dProgn, x, dg);
		f(4) += dg(1); f(5) += dg(2); f(6) += dg(3);
	}
	if (OnMoon)
	{
		SpeedMoon(tg, x(1, 3), aM);
		f(4) += aM(1); f(5) += aM(2); f(6) += aM(3);
	}
	if (OnSun)
	{
		SpeedSun(tg, x(1, 3), aS);
		f(4) += aS(1); f(5) += aS(2); f(6) += aS(3);
	}
	if (typeAtm > 0)
	{
		SpeedAtm(tg, x, dgAtm);
		f(4) += dgAtm(1); f(5) += dgAtm(2); f(6) += dgAtm(3);
	}
	if (OnLux)
	{
		dgLux = PressLux(tg, x);
		f(4) += dgLux(1); f(5) += dgLux(2); f(6) += dgLux(3);
	}
	STACKh(n, x)
	*(x->h) = f;
	return *(x->h);
}

//B7. Модель движения, учитывающая только коэфф. c20
//(нецентральность грав. поля за счет сжатия)
Vect& Ballist::sModProgn(double t, Vect& x)
{
	double c20 = -0.0010826257, tg, sf, cf, r = !x(1, 3),
		   r3 = r * r * r, coef = c20 * mu_ / r / r * (aE_ / r) * (aE_ / r);
	Vect xg(3), dg(3); Matr A(3, 3); VSTACK(6, f)
	if (t >= 86400 * (kday + 1.0))ud++; if (t < 86400.0 * kday)ud--;
	kday = t / 86400; tg = t - kday * 86400.0; A = Mat2000Gr(ud, t);
	xg = A * x(1, 3); sf = xg(3) / r; cf = !xg(1, 2) / r;
	dg(1) = -1.5 * coef * (3 * sf * sf - 1); dg(2) = 3 * coef * sf * cf;
	dg = !MatDecSpher(xg) * dg; dg = !A * dg;
	f(1) = x(4); f(2) = x(5); f(3) = x(6);
	f(4) = -mu_ * x(1) / r3 + dg(1);
	f(5) = -mu_ * x(2) / r3 + dg(2);
	f(6) = -mu_ * x(3) / r3 + dg(3);
	return f;
}

//B8. Определение компонент ускорения, обусловленного аномалиями 
//гравитационного поля, в J2000	
//(коэффициенты гравитационного потенциала - ненормированные) 
//Входные параметры:                                                        
//t                 - гринвичское время в с;                       
//c((n+1)*(n+2)/2),                                                        
//d((n+1)*(n+2)/2)	- коэффициенты разложения геопотенциала        
//                    в ряд по сферическим функциям;               
//x(6)              - ВС в J2000									 
//Выходные параметры:                                                       
//dg                - проекции вектора ускорения на оси J2000		 
void Ballist::GravJ2000(double t, Vect& c, Vect& d, Vect& x, Vect& dg)
{
	Vect xg(3), axg(3); 
	Matr A(3, 3); 
	A = Mat2000Gr(ud, t);
	xg = A * x(1, 3);
	GravGrinv(c, d, xg, axg); 
	dg = !A * axg;
}

//B9. Определение компонент ускорения, обусловленного аномалиями     
//гравитационного поля, в гринвичской системе координат          
//(коэффициенты гравитационного потенциала - ненормированные)    
//Входные параметры:                                                       
//c((n+1)*(n+2)/2),d((n+1)*(n+2)/2)	- коэффициенты разложения геопотенциала   
//									  в ряд по сферическим функциям;
//x(3)								- гринвичские координаты                       
//Выходные параметры:                                                      
//dg								- проекции вектора ускорения на оси            
//									  гринвичской системы координат              
void Ballist::GravGrinv(Vect& c, Vect& d, Vect& x, Vect& dg)
{
	double R, r, cF, tF; Vect s(3), as(3);
	R = !x(1, 2); r = !x; cF = R / r; tF = x(3) / R;
	s(1) = r; s(2) = atan(tF); s(3) = atan2(x(2), x(1)); //Сферические координаты
	as = GravSphere(c, d, s); dg = !MatDecSpher(s(2), s(3)) * as;
}

//B10. Определение частных производных от геопотенциала	   
//по сферическим координатам					
//(коэффициенты геопотенциала - ненормированные)      
//Входные параметры:                                                         
//c((n+1)*(n+2)/2),	d((n+1)*(n+2)/2) - коэффициенты разложения геопотенциала
//           						   в ряд по сферическим функциям;
//s									 - сферические координаты
//									   (радиус, широта, долгота);                    
//Выходные параметры:                                                        
//dg(3)								 - производные от геопотенциала по радиусу,
//									   широте и долготе соответственно
Vect& Ballist::GravSphere(Vect& c, Vect& d, Vect& s)
{
	int i, j, k = 0;
	double R, F, L, sf = sin(s(2)), cf = cos(s(2)), tf = tan(s(2)),
		sL = sin(s(3)), cL = cos(s(3)), mur, S, h, jL, sjL, cjL;
	Matr P(Ngrav + 2, Ngrav + 2); VSTACK(3, dg)
	Legendr(Ngrav + 1, sf, P); mur = mu_ / s(1) / s(1);
	//Pij(x) - присоединенные функции Лежандра 1-го рода
	DO(i, 2, Ngrav)
	{
		R = 0; F = 0; L = 0;
		DO(j, 0, i)
		{
			k++; jL = j * s(3);
			sjL = sin(jL); cjL = cos(jL);
			S = c(k) * cjL + d(k) * sjL;
			R += S * P(i + 1, j + 1);
			F += S * (P(i + 1, j + 2) - j * tf * P(i + 1, j + 1));
			L += j * (d(k) * cjL - c(k) * sjL) * P(i + 1, j + 1);
		}
		h = pow(aE_ / s(1), i);
		dg(1) += (i + 1.0) * h * R; dg(2) += h * F; dg(3) += h * L;
	}
	dg(1) = -mur * dg(1); dg(2) = mur * dg(2); dg(3) = mur * dg(3) / cf;
	return dg;
}

//B11. Определение компонент ускорения,
//обусловленного притяжением Луны.      
//Входные параметры:                                                     
//t		- гринвичское время;                                       
//x		- координаты орбитального объекта в J2000				     
//Выходные параметры:                                                    
//dg	- проекции вектора ускорения на оси J2000 				
void Ballist::SpeedMoon(double t, Vect& x, Vect& dg)
{
	double udt, r, d; Vect q(3), qx(3); udt = ud + t / 86400;
	Moon(udt, q); qx = q - x; d = !qx; r = !q; r = r * r * r; d = d * d * d;
	dg = muM_ * (qx / d - q / r);
}

//B12. Определение компонент ускорения,
//обусловленного притяжением Солнца     
//Входные  :                                                                 
//t		- гринвичское время в с;                                      
//x		- координаты орбитального объекта в J2000					  
//Выходные параметры:                                                        
//dg    - проекции вектора ускорения на оси J2000					  
void Ballist::SpeedSun(double t, Vect& x, Vect& dg)
{
	double udt, r, d; Vect q(3), qx(3); udt = ud + t / 86400;
	Sun(udt, q);	qx = q - x; d = !qx; r = !q; r = r * r * r; d = d * d * d;
	dg = muS_ * (qx / d - q / r);
}

//B13. Учет атмосферного торможения в J2000.
//t		- гринвичское время в с;
//x		- фазовый вектор в J2000;
//dg	- вектор возмущающего ускорения в проекциях на оси J2000
void Ballist::SpeedAtm(double t, Vect& x, Vect& dg)
{
	//double ro,khog=0.99371887105960277;
	double ro, khog = 1;
	Vect v(3), w(3), xg(3);
	//v - скорость КА отн. атмосферы;
	//w - скорость движ. атмосферы относительно J2000;
	//khog - коэфф. захвата атмосферы
	w(1) = -omE_ * x(2); w(2) = omE_ * x(1); v = x(4, 6) - khog * w;
	if (typeAtm == 1) { double h = !x(1, 3) - rE_; ro = StAtmos(h); }
	if (typeAtm == 2) { xg = Mat2000Gr(ud, t) * x(1, 3); ro = Atmos(t, xg); }
	dg = -BalCoef * ro * !v * v;
}

//B14. Светотеневая обстановка (Lux=1 - свет, Lux=0 - в тени).
//Точечный источник света.
//t	- гринвичское время в с; x - фазовый вектор в J2000
bool Ballist::FactorLux(double t, Vect& x)
{
	bool Lux; double L, L2, dL, alfa0, alfa, ca, R2 = aE_ * aE_;
	//Vect r(3),rS(3); r=x(1,3); rS=Soleil(ud,t); 
	Vect r(3), rS(3); r = x(1, 3); Sun(ud + t / 86400, rS);
	L2 = rS * rS; L = sqrt(L2); dL = (L2 - R2) / L;
	ca = (rS + r) * rS / !(rS + r) / !rS;
	alfa = acos(ca); alfa0 = asin(aE_ / L);
	Lux = !(rS + r) * ca<dL || alfa>alfa0;
	return Lux;
}

//B15. Светотеневые параметры. Точечный источник света.
//Вход: t - гринвичское время в с; x - фазовый вектор в J2000
//Выход: dalfa=alfa-alfa0,
//		 где alfa0 - угловой размер Земли (от Солнца);
//		 alfa - угловое положение КА (от Солнца).
void Ballist::FactorLux(double t, Vect& x, double& dalfa)
{
	double alfa0, alfa, ca, L, L2, dL, R2 = aE_ * aE_;
	//Vect r(3),rS(3); r=x(1,3); rS=Soleil(ud,t); 
	Vect r(3), rS(3); r = x(1, 3); Sun(ud + t / 86400, rS);
	L2 = rS * rS; L = sqrt(L2); dL = (L2 - R2) / L;
	ca = (rS - r) * rS / !(rS - r) / !rS;
	alfa = acos(ca); alfa0 = asin(aE_ / L);
	if (!(rS - r) * ca < dL) { dalfa = 1; return; }
	dalfa = alfa - alfa0;
}

//B16. Светотеневая обстановка (Lux=1 - свет, Lux=0 - в тени).
//Бесконечно удаленный источник света.
//t	- гринвичское время в с; x - фазовый вектор в J2000
bool Ballist::FactorLux(Vect& x, double t)
{
	bool Lux; double L, d, a;
	//Vect r(3),rS(3); r=x(1,3); rS=Soleil(ud,t); 
	Vect r(3), rS(3); r = x(1, 3); Sun(ud + t / 86400, rS);
	L = !rS; a = rS * r / L; if (a > 0)return 1;
	d = sqrt(r * r - a * a); Lux = d > aE_;
	return Lux;
}

//B17. Светотеневые параметры. Бесконечно удаленный источник света.
//Вход: t - гринвичское время в с; x - фазовый вектор в J2000
//Выход: deld - разность между длиной нормали от КА к вектору
//				Земля - Солнце и радиусом Земли
void Ballist::FactorLux(Vect& x, double t, double& deld)
{
	double L, d, a; Vect r(3), rS(3);
	//r=x(1,3); rS=Soleil(ud,t); 
	r = x(1, 3); Sun(ud + t / 86400, rS);
	L = !rS; a = rS * r / L; if (a > 0) { deld = 1e9; return; }
	d = sqrt(r * r - a * a); deld = d - aE_;
}

//B18. Ускорение в J2000, обусловленное световым давлением.
//t		- гринвичское время в с;
//x		- фазовый вектор в J2000
Vect& Ballist::PressLux(double t, Vect& x)
{
	Vect r(3), rS(3), r_rS(3); VSTACK(3, W)
	//r=x(1,3); rS=Soleil(ud,t); r_rS=r-rS; 
	r = x(1, 3); Sun(ud + t / 86400, rS); r_rS = r - rS;
	W = FactorLux(t, x) * gamma * p * Sm / mass * r_rS / !r_rS;
	return W;
}

//B19. Динамическая модель атмосферы (ГОСТ Р 25645.166-2004).
//Входные параметры: ud - юлианская дата;
//t - гринвичское время в с; x(3) - гринвичские координаты
double Ballist::Atmos(double t, Vect& x)
{
	double ro, ron, h, K0, K1, K2, K3, K4, cosfi, cosfi2, fi = 0, bet, S, rx;
	Vect rS(3), sS(3);
	Matr A(7, 7), B(5, 7), C(5, 7), D(5, 7), E(9, 7), Fi1(1, 7), L(5, 7), N(3, 7);
	rx = !x; h = Altitude(x); if (rx > 1e6)h /= 1000; if (startAtm)ud0Atm = ud;
	if (startAtm || ud != ud0Atm)
	{
		Vect AA(9); char s[50]; ReadVect(multA, s, AA);
		if (h >= 120 && h < 500)ReadMatr(A120, s, *Aatm);
		if (h >= 500)ReadMatr(A500, s, *Aatm);
		if (h >= 120 && h < 600)ReadMatr(B120, s, *Batm);
		if (h >= 600 && h < 660)ReadMatr(B600, s, *Batm);
		if (h >= 660 && h < 760)ReadMatr(B660, s, *Batm);
		if (h >= 760 && h < 800)ReadMatr(B760, s, *Batm);
		if (h >= 800 && h < 860)ReadMatr(B800, s, *Batm);
		if (h >= 860 && h < 900)ReadMatr(B860, s, *Batm);
		if (h >= 900 && h < 1000)ReadMatr(B900, s, *Batm);
		if (h >= 1000)ReadMatr(B1000, s, *Batm);
		if (h >= 120 && h < 640)ReadMatr(C120, s, *Catm);
		if (h >= 640 && h < 700)ReadMatr(C640, s, *Catm);
		if (h >= 700 && h < 760)ReadMatr(C700, s, *Catm);
		if (h >= 760 && h < 820)ReadMatr(C760, s, *Catm);
		if (h >= 820 && h < 860)ReadMatr(C820, s, *Catm);
		if (h >= 860 && h < 920)ReadMatr(C860, s, *Catm);
		if (h >= 920 && h < 980)ReadMatr(C920, s, *Catm);
		if (h >= 980)ReadMatr(C980, s, *Catm);
		if (h >= 120)ReadMatr(D120, s, *Datm);
		if (h >= 120 && h < 600)ReadMatr(E120, s, *Eatm);
		if (h >= 600 && h < 700)ReadMatr(E600, s, *Eatm);
		if (h >= 700 && h < 760)ReadMatr(E700, s, *Eatm);
		if (h >= 760 && h < 780)ReadMatr(E760, s, *Eatm);
		if (h >= 780 && h < 800)ReadMatr(E780, s, *Eatm);
		if (h >= 800 && h < 900)ReadMatr(E900, s, *Eatm);
		if (h >= 900)ReadMatr(E900, s, *Eatm);
		if (h >= 120)ReadMatr(Fi1120, s, *Fi1Atm);
		if (h >= 120 && h < 640)ReadMatr(L120, s, *Latm);
		if (h >= 640 && h < 660)ReadMatr(L640, s, *Latm);
		if (h >= 660 && h < 740)ReadMatr(L660, s, *Latm);
		if (h >= 740 && h < 800)ReadMatr(L740, s, *Latm);
		if (h >= 800 && h < 860)ReadMatr(L800, s, *Latm);
		if (h >= 860 && h < 900)ReadMatr(L860, s, *Latm);
		if (h >= 900)ReadMatr(L900, s, *Latm);
		if (h >= 120)ReadMatr(N120, s, *Natm);
		Date(ud, day, month, year); nday = ud - Ulian(1, 1, year);
		A_dAtm = Polynom(AA, nday);
		F107Atm = 200; F81Atm = 200; KpAtm = 2.6667; //Только для сверки (отд.116)
		//ISA_IGI(day,month,year,F107Atm,F81Atm,KpAtm);
		F0Atm = SelectF0(F81Atm);
		switch (F0Atm)
		{
			case 75:  indexAtm = 1; break; case 100: indexAtm = 2; break;
			case 125: indexAtm = 3; break; case 150: indexAtm = 4; break;
			case 175: indexAtm = 5; break; case 200: indexAtm = 6; break;
			case 250: indexAtm = 7;
		}
		multA.seekg(0);
		A120.seekg(0); A500.seekg(0);
		B120.seekg(0); B600.seekg(0); B660.seekg(0);
		B760.seekg(0); B800.seekg(0); B860.seekg(0);
		B900.seekg(0); B1000.seekg(0);
		C120.seekg(0); C640.seekg(0); C700.seekg(0);
		C760.seekg(0); C820.seekg(0); C860.seekg(0);
		C920.seekg(0); C980.seekg(0);
		D120.seekg(0);
		E120.seekg(0); E600.seekg(0); E700.seekg(0);
		E760.seekg(0); E780.seekg(0); E900.seekg(0);
		Fi1120.seekg(0);
		L120.seekg(0); L640.seekg(0); L660.seekg(0);
		L740.seekg(0); L800.seekg(0); L860.seekg(0);
		L900.seekg(0); N120.seekg(0);
	}
	A = *Aatm; B = *Batm; C = *Catm; D = *Datm;
	E = *Eatm; Fi1 = *Fi1Atm; L = *Latm; N = *Natm;
	Sun(ud + t / 86400, rS); sS = DecSpher(rS);
	S = StimeVer(ud, t); bet = sS(3) - S + Fi1(1, indexAtm);
	cosfi = (x(3) * sin(sS(2)) + cos(sS(2)) * (x(1) * cos(bet) + x(2) * sin(bet))) / rx;
	cosfi2 = sqrt(0.5 * (1 + cosfi));
	ron = ro0_ * exp(Polynom(A[indexAtm], h));
	K0 = 1 + Polynom(Fi1[indexAtm], h) * (F81Atm - F0Atm) / F0Atm;
	K1 = Polynom(C[indexAtm], h) * pow(cosfi2, Polynom(N[indexAtm], h));
	K2 = Polynom(D[indexAtm], h) * A_dAtm;
	K3 = Polynom(B[indexAtm], h) * (F107Atm - F81Atm) / (F81Atm + fabs(F107Atm - F81Atm));
	K4 = Polynom(E[indexAtm](1, 5), h) * Polynom(E[indexAtm](6, 9), KpAtm);
	ro = ron * K0 * (1 + K1 + K2 + K3 + K4);
	startAtm = 0; ud0Atm = ud;
	return ro;
}

//B20. Определение координат Луны в J2000
void Ballist::Moon(double udt, Vect& r)
{
	int i, j, NInt, NStr;
	double DTJ0, D_DTJ, DTJ, NI;
	static Matr C(1449, 10); Matr KF(10, 3);
	DTJ0 = 2449600; D_DTJ = udt - DTJ0;
	NInt = floor(D_DTJ / 20); NStr = NInt * 3;
	NI = NInt; DTJ = D_DTJ - NI * 20;
	if (startMoon)
	{
		ReadMatr(InMoon, C);
		startMoon = 0;
	}
	DO(i, 1, 10)DO(j, 1, 3)KF(i, j) = C(NStr + j)(i);
	if (muM_ > 1e9)PolinomM(DTJ, KF * 1000, r);
	else PolinomM(DTJ, KF, r);
}

//B20. Определение координат Солнца в J2000
void Ballist::Sun(double udt, Vect& r)
{
	int NInt, NStr, i, j;
	double DTJ0, D_DTJ, DTJ, NI;
	static Matr C(291, 11); Matr KF(11, 3);
	DTJ0 = 2449600; D_DTJ = udt - DTJ0;
	NInt = floor(D_DTJ / 100); NStr = NInt * 3;
	NI = NInt; DTJ = D_DTJ - NI * 100;
	if (startSun)
	{
		ReadMatr(InSun, C);
		startSun = 0;
	}
	DO(i, 1, 11)DO(j, 1, 3)KF(i, j) = C(NStr + j)(i);
	if (muS_ > 1e15)PolinomS(DTJ, KF * 1000, r);
	else PolinomS(DTJ, KF, r);
}

//B22. Просмотр компонент класса Ballist
void Ballist::Look()
{
	Ecrir(cout, "\n Ballist::Look():\n------------------",
		  "\n ud=", ud, "Ngrav=", Ngrav, "OnMoon=", OnMoon,
		  "OnSun=", OnSun, "OnLux=", OnLux,
		  "\n typeAtm=", typeAtm, "BalCoef=", BalCoef,
		  "mass=", mass, "Sm=", Sm);
}

//										ГЛОБАЛЬНЫЕ ФУНКЦИИ

//G1. Выделение памяти под одномерный массив
double *Array(int m)
{
	double *a = new double[m];
	if (m<0 || a == 0)
	{
		Ecrir(cout, "\n Array1: Mistake of separation to memories!");
		Ecrir(cout, "\n m=", m);
		return a;
	}
	int i;
	//Инициализация:
	DO(i, 0, m-1) a[i] = 0;
	return a;
}

//G2. Выделение памяти под двумерный массив
double **Array(int m, int n)
{
	int i, j; double **A = 0;
	A = new double*[m];
	DO(i, 0, m-1)A[i] = new double[n];
	if (m<0 || n<0 || A == 0)
	{
		Ecrir(cout, "\n Array2: Mistake of separation to memories!");
		Ecrir(cout, "\n m=", m, "    n=", n);
		return A;
	}
	//Инициализация:
	DO(i, 0, m-1) DO(j, 0, n-1) A[i][j] = 0;
	return A;
}

//G3. Выделение памяти под трехмерный массив
double ***Array(int m, int n, int l)
{
	int i = 0, j = 0, k = 0;
	double ***A = 0;
	A = new double**[m];
	DO(i, 0, m-1) A[i] = new double*[n];
	DO(i, 0, m - 1) DO(j, 0, n - 1) A[i][j] = new double[l];
	if (m<0 || n<0 || l<0 || A == 0)
	{
		cout << "	Array3: Mistake of separation to memories!/n";
		return A;
	}
	//Инициализация:
	DO(i, 0, m - 1) DO(j, 0, n - 1) DO(k, 0, l-1) A[i][j][k] = 0;
	return A;
}

Vect& Vector_e(int n, int k) //G4. Единичный вектор
{
	e->h = new Vect(n);	e->NumbVect++; //Поскольку e имеет тип Vect, 
	*(e->h) = 0; (*(e->h))(k) = 1;	   //то e->h относится к классу Vect, а не Matr	
	return *(e->h);
}

Matr& MatrixI(int n) //G5. Единичная матрица
{
	int i;
	E->H = new Matr(n, n);		//Поскольку E имеет тип Matr, то E->H относится к классу Matr,
	*(E->H) = 0; E->NumbMatr++; //а не Vect
	DO(i, 1, n)(*(E->H))(i, i) = 1;
	return *(E->H);
}

//G6. Выделение памяти под массив символов
char *Line(int m)
{
	int i;
	char *a = new char[m + 1];
	//Инициализация:
	DO(i, 0, m - 1) a[i] = ' '; a[m] = '\0';
	return a;
}

//G7. Выделение памяти под одномерный массив [m] из ns символов 
char **Line(int m, int ns)
{
	int i, j; char **s = 0;
	s = new char*[m];
	DO(i, 0, m - 1) s[i] = new char[ns + 1];
	//Инициализация:
	DO(i, 0, m - 1)
	{
		DO(j, 0, ns-1) s[i][j] = ' ';
		s[i][ns] = '\0';
	}
	return s;
}

//G8. Выделение памяти под двумерный массив [m x n] из ns символов
char ***Line(int m, int n, int ns)
{
	int i, j, k; char ***s = 0;
	s = new char**[m];
	DO(i, 0, m - 1) s[i] = new char*[n];
	DO(i, 0, m - 1) DO(j, 0, n - 1) s[i][j] = new char[ns + 1];
	//Инициализация:
	DO(i, 0, m - 1) DO(j, 0, n - 1)
	{
		DO(k, 0, ns-1) s[i][j][k] = ' ';
		s[i][j][ns] = '\0';
	}
	return s;
}

//G9. Выделение памяти под двумерный массив [m x n x l] из ns символов
char ****Line(int m, int n, int l, int ns)
{
	int i, j, k, p; char ****s = 0;
	s = new char***[m];
	DO(i, 0, m - 1) s[i] = new char**[n];
	DO(i, 0, m - 1) DO(j, 0, n - 1) s[i][j] = new char*[l];
	DO(i, 0, m - 1) DO(j, 0, n - 1) DO(k, 0, l - 1)	s[i][j][k] = new char[ns + 1];
	//Инициализация:
	DO(i, 0, m - 1) DO(j, 0, n - 1) DO(k, 0, l - 1)
	{
		DO(p, 0, ns-1) s[i][j][k][p] = ' ';
		s[i][j][k][ns] = '\0';
	}
	return s;
}

//G10. Приведение к цифре (если x - не цифра, а составное число)
char Lit(int x)
//x - цифра от 0 до 9
{
	int k, L;
	if (x != 0) k = (int)log10(x); if (x == 0) k = 1;
	L = (int)(x / pow(10, k)); //Первая цифра числа x
	char s[11] = "0123456789";
	return s[L];
}

//G11. Перевод числа x в строку из m символов
char *Liter(double x, int m)
{
	char *c = Line(m);
	if (x == 0) { c[m / 2] = '0'; return c; }
	int j0 = 0, j = 0, jj = 0, jf = 0, k = 0, L = 0, L1 = 0, iz = 0, ka = 0, kk = 0;
	double z=0, y=0; 
	z = fabs(x);
	if (z != 0)k = (int)log10(z);
	if (z<1)k = k - 1;
	if (z == 0)k = 1; if (k<-75)k = -75;
	L = (int)(z / pow(10, k)); //Первая цифра числа x
	if (abs(k)<9)
	{
		iz = (int)z; j0 = (m - k) / 2;
		if (z - iz == 0)
		{
			if (x<0)c[j0 - 1] = '-';
			y = z / pow(10, k + 1); L1 = 0;
			DO(j, j0, j0 + k)
			{
				y = (y - L1) * 10; L1 = (int)(y + 0.000001);
				c[j] = Lit(L1);
			}
			return c;
		}
	}
	if (x<0)c[0] = '-';
	if (z<pow(10, 5 - m) || z>pow(10, m - 2))
	{
		c[1] = Lit(L); c[2] = '.'; jf = m - 7;
		y = (z - L * pow(10, k)) / pow(10, k); L1 = 0;
		for (j = 0; j<jf; j++) { y = (y - L1) * 10; L1 = (int)y; c[3 + j] = Lit(L1); }
		j = 3 + jf; c[j] = 'e'; j++; c[j] = '+'; if (k<0)c[j] = '-'; ka = abs(k);
		kk = ka / 10; c[j + 1] = Lit(kk); kk = ka - kk * 10; c[j + 2] = Lit(kk);
	}
	if (z >= pow(10, 5 - m) && z <= 1)
	{
		c[1] = '0'; c[2] = '.';
		jf = m - 3; y = z; L1 = 0;
		DO(j, 0, jf-1) { y = (y - L1) * 10; L1 = (int)y; c[3 + j] = Lit(L1); }
	}
	if (z>1 && z <= pow(10, m - 2))
	{
		c[1] = Lit(L); y = (z - L * pow(10, k)) / pow(10, k); L1 = 0;
		DO(j, 0, k-1) { y = (y - L1) * 10; L1 = (int)y; c[2 + j] = Lit(L1); }
		jj = 3 + k; c[jj - 1] = '.'; jf = m - 3 - k;
		DO(j, 0, jf-1) { y = (y - L1) * 10; L1 = (int)y; c[jj + j] = Lit(L1); }
	}
	return c;
}

//G12. Перевод символьной строки с в строку из m символов
char* Liter(const char* c, int m)
{
	char* s = Line(m); int j0, j, L = strlen(c);
	if (L >= m) { DO(j, 0, m - 1)s[j] = c[j]; return s; }
	j0 = (int)((m - L) / 2) + 1;
	DO(j, j0, j < j0 + L - 1)s[j] = c[j - j0];
	return s;
}

//G13. Перевод числа x в строку из m символов
void Liter(double x, int m, char *c)
{
	//if(x==0){c[m/2]='0'; return;}
	int j0 = 0, j = 0, jj = 0, jf = 0, k = 0, L = 0, L1 = 0, iz = 0, ka = 0, kk = 0;
	double z = 0, y = 0;
	z = fabs(x);
	if (z != 0)k = (int)log10(z); if (z<1)k = k - 1;
	if (z == 0)k = 1; if (k<-75)k = -75;
	L = (int)(z / pow(10, k)); //Первая цифра числа x
	if (abs(k)<9)
	{
		iz = (int)z; j0 = (m - k) / 2;
		if (z - iz == 0)
		{
			if (x<0)c[j0 - 1] = '-';
			y = z / pow(10, k + 1); L1 = 0;
			for (j = j0; j<j0 + k + 1; j++)
			{
				y = (y - L1) * 10; L1 = (int)(y + 0.000001); c[j] = Lit(L1);
			}
			return;
		}
	}
	if (x<0) c[0] = '-';
	if (z<pow(10, 5 - m) || z>pow(10, m - 2))
	{
		c[1] = Lit(L); c[2] = '.'; jf = m - 7;
		y = (z - L * pow(10, k)) / pow(10, k); L1 = 0;
		DO(j, 0, jf-1) { y = (y - L1) * 10;
		L1 = (int)y; c[3 + j] = Lit(L1); }
		j = 3 + jf; c[j] = 'e'; j++; c[j] = '+';
		if (k<0)c[j] = '-'; ka = abs(k); kk = ka / 10;
		c[j + 1] = Lit(kk); kk = ka - kk * 10; c[j + 2] = Lit(kk);
	}
	if (z >= pow(10, 5 - m) && z <= 1)
	{
		c[1] = '0'; c[2] = '.';
		jf = m - 3; y = z; L1 = 0;
		DO(j, 0, jf-1) { y = (y - L1) * 10; L1 = (int)y; c[3 + j] = Lit(L1); }
	}
	if (z>1 && z <= pow(10, m - 2))
	{
		c[1] = Lit(L); y = (z - L * pow(10, k)) / pow(10, k); L1 = 0;
		DO(j, 0, k-1) { y = (y - L1) * 10; L1 = (int)y; c[2 + j] = Lit(L1); }
		jj = 3 + k; c[jj - 1] = '.'; jf = m - 3 - k;
		DO(j, 0, jf-1) { y = (y - L1) * 10; L1 = (int)y; c[jj + j] = Lit(L1); }
	}
}

//G14. Перевод символьной строки с в строку из m символов
void Liter(const char *c, int m, char *s)
{
	int j0, j, L = strlen(c);
	if (L >= m) { DO(j, 0, m-1)s[j] = c[j]; return; }
	j0 = (int)((m - L) / 2) + 1; 
	DO(j, j0, j0 + L -1)s[j] = c[j - j0];
}

Vect& Abs(Vect& x) //G15. 
{
	int m = x->m, i; 
	x->h = new Vect(m); x->NumbVect++;
	DO(i, 1, m)(*(x->h))(i) = fabs(x(i));
	return *(x->h);
}

Matr& Abs(Matr& A) //G16.
{
	int m = A->m, n = A->n, i, j; 
	A->H = new Matr(m, n); A->NumbMatr++;
	DO(i, 1, m)DO(j, 1, n) (*(A->H))(i, j) = fabs(A(i, j));
	return *(A->H);
}

//G17. Метод конфигураций	
//dx - начальный шаг по x;
//x0 - начальное приближение оптимального значения,
//	   возвращаемое значение - argmin Func(x)
Vect& MinFunc(double Func(Vect& x), int Nit, double eps, Vect& dx, Vect& x0)
{
	int i, eval, n = x0->m; bool On;
	double Fx, F, Fz, teta;
	Vect s(n), z(n), x(n);
	x0->h = new Vect(n); x0->NumbVect++;
	s = dx; x = x0;	Fx = Func(x0); eval = 1;
Lab1: F = Fx; z = x;
	DO(i, 1, n)
	{
		z(i) += s(i); Fz = Func(z); eval++;
		if (Fz < F)F = Fz;
		else
		{
			s(i) = -s(i); z(i) += 2 * s(i);
			if (eval < Nit)eval++; else { *(x0->h) = x; return *(x0->h); }
			Fz = Func(z); if (Fz < F)F = Fz; else z(i) -= s(i);
		}
	}
	if (F < Fx)
	{
	Lab2:
		DO(i, 1, n)
		{
			On = z(i) > x(i) && s(i) < 0 || z(i) <= x(i) && s(i) >= 0;
			if (On)s(i) = -s(i);
			teta = x(i); x(i) = z(i); z(i) = 2 * z(i) - teta;
		}
		Fx = F;
		if (eval < Nit)eval++; else { *(x0->h) = x; return *(x0->h); }
		Fz = Func(z); F = Fz;
		DO(i, 1, n)
		{
			z(i) += s(i); Fz = Func(z); eval++;
			if (Fz < F)F = Fz;
			else
			{
				s(i) = -s(i); z(i) += 2 * s(i);
				if (eval < Nit)eval++; else { *(x0->h) = x; return *(x0->h); }
				Fz = Func(z);
				if (Fz < F)F = Fz; else z(i) -= s(i);
			}
		}
		if (F >= Fx)goto Lab1;
		DO(i, 1, n)if (fabs(z(i) - x(i)) > 0.5 * fabs(s(i)))goto Lab2;
	}
	if (eval > Nit) { *(x0->h) = x; return *(x0->h); }
	DO(i, 1, n)if (fabs(s(i)) >= eps * dx(i))s(i) *= 0.25;
	goto Lab1;
}

Vect& Exp(Vect& x) //G18. Векторная экспонента
{
	int i, n = x->m;
	x->h = new Vect(n);  x->NumbVect++;
	DO(i, 1, n)(*(x->h))(i) = exp(x(i));
	return *(x->h);
}

double Sum(Matr& A) //G19. Сумма всех элементов матрицы
{
	int i, j; double S = 0;
	DO(i, 1, A->m)DO(j, 1, A->n) S += A(i, j);
	return S;
}

Matr& Dmat(int n, double a) //G20. Формирование диагональной матрицы по скаляру 
{
	int i;
	E->H = new Matr(n, n); E->NumbMatr++;
	DO(i, 1, n)(*(E->H))(i, i) = a;
	return *(E->H);
}

Matr& Dmat(Vect& a) //G21. Формирование диагональной матрицы по вектору
{
	int i, n = a->m;
	a->H = new Matr(n, n); a->NumbMatr++;
	DO(i,1,n)(*(a->H))(i,i)=a(i);
	return *(a->H);
}

double Sum(Vect& a) //G22. Сумма всех элементов вектора
{
	int i;  double s = 0;
	DO(i, 1, a->m)s += a(i);
	return s;
}

//G23. Определение ранга матрицы, номеров 
//линейно независимых строк s и линейно 
//независимых столбцов c
int RangMatr(Matr& A, Vect& s, Vect& c, double eps)
{
	bool flag = 1;
	int i, j, k, mn, mx, m1, m2, m = A->m, n = A->n, rang = 0;
	double detM, d0 = 0;
	if (TestZero(A))return 0;
	mn = min(m, n); mx = max(m, n);
	dO(k, mn, 1, -1)
	{
		m1 = Cnm(mn, k); m2 = Cnm(mx, k);
		Matr C1(m1, k), C2(m2, k), M(k, k);
		C1 = Comb(mn, k); C2 = Comb(mx, k);
		DO(i, 1, m1)DO(j, 1, m2)
		{
			if (m <= n)M = A(C1(i), C2(j));
			else M = A(C2(j), C1(i));
			detM = fabs(det(M));
			if (detM > eps)
			{
				if (flag) { rang = k; flag = 0; }
				if (detM > d0 && k == rang)
				{
					d0 = detM;
					//s = C1(i); c = C2(j);
					if (m <= n)s = C1(i); else c = C1(i);
					if (n > m)c = C2(j); else s = C2(j);
				}
			}
		}
		if (k == rang)break;
	}
	return rang;
}

//G24. Определение ранга матрицы
int RangMatr(Matr& A, double eps)
{
	int m = A->m, n = A->n;
	Vect s(m), c(n);
	return RangMatr(A, s, c, eps);
}

//G25. Определение собственных значений и собственных
//векторов матрицы A.
//Входные параметры:
//F		- функция, определяющая характеристическое
//		  уравнение det(A-lam*MatrixI(n))=0;
//Nit	- число итераций для решения уравнения F(lam) = 0;
//eps	- наименьшая величина базисного минора,
//		  при которой он считается отличным от нуля
//		  (для определения ранга матрицы в функции NotTrivial)
//A(n,n)- заданная матрица;
//Выходные параметры:
//lam	- вектор собственных значений.
//Возвращаемое значение:
//V(n,n)- матрица, i-й столбец которой - собствнный
//		  вектор, соответствующий lam(i), i=1,...,n
Matr& OwnVector(double F(double lam), int Nit,
	double eps, Matr& A, Vect& lam)
{
	int i, n = A->n; //Определение размерности матрицы A(n,n)
	Matr V(n, n);
	Vect v(n); Matr B(n, n);
	ZeroFun(F, Nit, eps, lam); //Определение всех корней
							   //(собственных значений lam(n))
							   //характеристического уравнения
	DO(i, 1, n)
	{
		B = A - lam(i) * MatrixI(n); //Матрица однородной системы уравнений
									 //[A - lam(i) * MatrixI(n)]*v=0
									 //для определения i-го собственного вектора
		v = NotTrivial(eps, B); //Частное нетривиальное решение однородной системы уравнений
								//[A - lam(i) * MatrixI(n)]*v=0 относительно вектора v
								//(в качестве свободной переменной выбрана v(n) - последняя 
								//компонента вектора v, частное значение которой принято
								//равной единице)
		V(v, i); //Вставка вектора v в i-й столбец матрицы V
	}
	A->H = new Matr(n, n); A->NumbMatr++;
	*(A->H) = V;
	return *(A->H);
}

//G26. Частное нетривиальное решение линейной однородной системы уравнений
Vect& NotTrivial(double eps, Matr& A)
{
	int n = A->n, r, nr;
	Vect x(n);
	r = RangMatr(A, eps);
	nr = n - r; Vect x1(r), x2(nr);
	x2 = Vector_e(nr, nr);
	e->h = new Vect(n); e->NumbVect++; //Только после x2 = Vector_e(nr, nr); 
									   //Иначе размерность *(e->h) станет равной nr 
	if (r == n && n == A->m) { *(e->h) = 0; return *(e->h); }
	x1 = -(~A(1, r, 1, r) * A(1, r, r + 1, n) * x2);
	x = Adding(x1, x2); *(e->h) = x;
	return *(e->h);
}

//G27. Определение ранга матрицы B(m,n),
//номеров линейно независимых строк row
//и столбцов col, 1e-15<=eps<=1e-9 (рекомендуется)
int RangMatr(double eps, Matr& B, Vect& row, Vect& col)
{
	if (B->TestZero())return 0;
	int i, j, l, ii, jj, kk, ll, mm, r, c, nm, ncol, rang, m = B->m, n = B->n;
	double piv, hold, sav, tol; Vect A(m * n);
	DO(i, 1, n)DO(j, 1, m)A((i - 1) * m + j) = B(j, i);
	if (m <= 0 || n <= 0) { rang = -1; goto L3; }
	rang = 0; piv = 0; jj = 0;
	DO(j, 1, n)
	{
		col(j) = j;
		DO(i, 1, m)
		{
			jj++; hold = A(jj);
			if (fabs(piv) - fabs(hold) < 0) { piv = hold; r = i; c = j; }
		}
	}
	DO(i, 1, m)row(i) = i; tol = fabs(eps * piv); nm = n * m;
	Do(ncol, m, nm, m)
	{
		if (fabs(piv) - tol <= 0)break; rang++; jj = r - rang;
		if (jj > 0)
		{
			Do(j, rang, nm, m) { i = j + jj; sav = A(j); A(j) = A(i); A(i) = sav; }
			jj = int(row(r)); row(r) = row(rang); row(rang) = jj;
		}
		jj = (c - rang) * m;
		if (jj > 0)
		{
			kk = ncol;
			DO(j, 1, m) { i = kk + jj; sav = A(kk); A(kk) = A(i); kk--; A(i) = sav; }
			jj = int(col(c)); col(c) = col(rang); col(rang) = jj;
		}
		kk = rang + 1; mm = rang - m; ll = ncol + mm; if (mm >= 0)goto L25;
		jj = ll; sav = piv; piv = 0;
		DO(j, kk, m)
		{
			jj++; hold = A(jj) / sav; A(jj) = hold; l = j - rang;
			if (rang >= n)break; ii = jj;
			DO(i, kk, n)
			{
				ii += m; mm = ii - l; A(ii) = A(ii) - hold * A(mm);
				if (fabs(A(ii)) - fabs(piv) >= 0) { piv = A(ii); r = j; c = i; }
			}
		}
	}
	if (rang < 1)goto L3; if (rang == 1)goto L25; r = ll;
	DO(j, 2, rang)
	{
		ii = j - 1; r -= m; jj = ll;
		DO(i, kk, m)
		{
			hold = 0; jj++; mm = jj; c = r;
			DO(l, 1, ii) { hold = hold + A(mm) * A(c); c--; }
			mm -= m;
			A(mm) = A(mm) - hold;
		}
	}
L25:
	if (n <= rang)goto L3; r = ll; kk = ll + m;
	DO(j, 1, rang)
	{
		Do(i, kk, nm, m)
		{
			jj = r; ll = i; hold = 0; ii = j;
		L27:
			ii--;
			if (ii > 0) { hold -= A(jj) * A(ll); jj -= m; ll--; goto L27; }
			A(ll) = (hold - A(ll)) / A(jj);
		}
	}
	r--; DO(i, rang + 1, m)row(i) = 0; DO(i, rang + 1, n)col(i) = 0;
L3:
	if (rang > 0)
	{
		Vect rc(rang); rc = row(1, rang); Sort1(rc); row = rc;
		rc = col(1, rang); Sort1(rc); col = rc; return rang;
	}
	return rang;
}

//G28. Определение ранга RangMatr матрицы A(m,n)
//1e-15<=eps<=1e-9 (рекомендуется)
int RangMatr(double eps, Matr& A)
{
	int m = A->m, n = A->n;
	Vect Am(m), An(n);
	return RangMatr(eps, A, Am, An);
}

//G29. Определение нулей функции скалярного аргумента.
//Nit - максимально допустимое число итераций;
//eps1 - относительный критерий сходимости итераций;
//eps2 - абсолютный критерий сходимости по значениям функции;
//eps3, eta - используются для кратных корней так, что если
//            |r-c(i)|<eps3, то r заменяется на r+eta;
void ZeroFun(double F(double t), int Nit, double eps1,
	double eps2, double eps3, double eta, Vect& c)
{
	int i, j, k, n = c.m;
	double p, p1, p2, x0, x1, x2, r, f, f1, d, e, h, u, v, w, t;
	DO(k, 1, n)
	{
		j = 0; p = 0.9 * c(k); p1 = 1.1 * c(k); p2 = c(k);
		if (c(k) == 0) { p = -1; p1 = 1; p2 = 0; }
		r = p; goto reg;
	s1:		r = p1; x0 = f1; goto reg;
	s2:		r = p2; x1 = f1; goto reg;
	s3:		x2 = f1; d = -0.5; h = c(k) == 0 ? -1 : -0.1 * c(k);
	iter:	e = d + 1; t = x0 * d * d - x1 * e * e + x2 * (d + e);
		u = t * t - 4 * x2 * d * e * (x0 * d - x1 * e + x2); u = u <= 0 ? 0 : sqrt(u);
		v = t + u; w = t - u; u = fabs(v) <= fabs(w) ? w : v; if (u == 0)u = 1;
		t = -2 * x2 * e / u; h = t * h; r += h; if (fabs(h / r) < eps1)goto calc; else goto reg;
	s4:		if (fabs(f1) < fabs(x2 * 10)) { x0 = x1; x1 = x2; x2 = f1; d = t; goto iter; }
		t *= 0.5; h *= 0.5; r -= h;
	reg:	j++;
	calc:	f = F(r); f1 = f; if (j >= Nit)goto fin;
		DO(i, 2, k)
		{
			e = r - c(i - 1);
			if (fabs(e) >= eps3)f1 /= e; else { r += eta; goto calc; }
		}
		if (fabs(f) < eps2 && fabs(f1) < eps2)goto fin;
		else
		{
			if (j == 1)goto s1; if (j == 2)goto s2;
			if (j == 3)goto s3; if (j > 3)goto s4;
		}
	fin:	c(k) = r;
	}
}

//G30. Определение нулей функции скалярного аргумента.
//Упрощённое обращение к ZeroFun.
//Входные параметры:
//F(t)	- функция, нули которой определяются; 
//Nit	- максимально допустимое число итераций;
//eps	- обобщенная точность сходимости итераций;
//c		- начальное приближение к вектору нулей функции F(t)
//Выходные параметры:
//c		- вектор, для i-й компоненты которого F[c(i)]=0.  
void ZeroFun(double F(double t), int Nit, double eps, Vect& c)
{
	ZeroFun(F, Nit, eps, eps, eps, eps, c);
}

void Sort1(Vect& x) //G31. Сортировка массива по возрастанию
{
	double t; int m = x->m, i, j, k;
	DO(i, 2, m)
	{
		k = i; t = x(k);
		dO(j, i - 1, 1, -1)
			if (x(j) <= t)break;
			else { k = j; x(j + 1) = x(j); }
		x(k) = t;
	}
}

void Sort2(Vect& x) //G32. Сортировка массива по убыванию
{
	double t; int m = x->m, i, j, k;
	DO(i, 2, m)
	{
		k = i; t = x(k);
		dO(j, i - 1, 1, -1)
			if (x(j) >= t)break;
			else { k = j; x(j + 1) = x(j); }
		x(k) = t;
	}
}

Vect& Sort3(Vect& x) //G33. Сортировка массива по возрастанию с указанием номеров
{
	double t; int m = x->m, i, j, k;
	Vect x0(m), N(m);
	x0 = x;
	DO(i, 2, m)
	{
		k = i; t = x(k);
		dO(j, i - 1, 1, -1)
			if (x(j) <= t)break;
			else { k = j; x(j + 1) = x(j); }
		x(k) = t;
	}
	k = 0;
	DO(i, 1, m)DO(j, 1, m)
		if (x(i) == x0(j)) { k++; N(k) = j; break; }
	STACKh(m, x)
	*(x->h) = N;
	return *(x->h);
}

Vect& Sort4(Vect& x) //G34. Сортировка массива по убыванию с указанием номеров
{
	double t; int m = x->m, i, j, k;
	Vect x0(m), N(m);
	x0 = x;
	DO(i, 2, m)
	{
		k = i; t = x(k);
		dO(j, i - 1, 1, -1)
			if (x(j) >= t)break;
			else { k = j; x(j + 1) = x(j); }
		x(k) = t;
	}
	k = 0;
	DO(i, 1, m)DO(j, 1, m)
		if (x(i) == x0(j)) { k++; N(k) = j; break; }
	STACKh(m, x)
	*(x->h) = N;
	return *(x->h);
}

//G35. Определение номеров компонент массива по возрастанию.
//Исходный вектор сохраняется
Vect& Sort5(Vect& x)
{
	double t; int m = x->m, i, j, k;
	Vect x0(m), N(m);
	x0 = x;
	DO(i, 2, m)
	{
		k = i; t = x(k);
		dO(j, i - 1, 1, -1)
			if (x(j) <= t)break;
			else { k = j; x(j + 1) = x(j); }
		x(k) = t;
	}
	k = 0;
	DO(i, 1, m)DO(j, 1, m)
		if (x(i) == x0(j)) { k++; N(k) = j; break; }
	x = x0;
	STACKh(m, x)
	*(x->h) = N;
	return *(x->h);
}

//G36. Определение номеров компонент массива по убыванию
//Исходный вектор сохраняется
Vect& Sort6(Vect& x)
{
	double t; int m = x.m, i, j, k;
	Vect x0(m), N(m);
	x0 = x;
	DO(i, 2, m)
	{
		k = i; t = x(k);
		dO(j, i - 1, 1, -1)
			if (x(j) >= t)break;
			else { k = j; x(j + 1) = x(j); }
		x(k) = t;
	}
	k = 0;
	DO(i, 1, m)DO(j, 1, m)
		if (x(i) == x0(j)) { k++; N(k) = j; break; }
	x = x0;
	STACKh(m, x)
	*(x->h) = N;
	return *(x->h);
}

bool TestZero(Vect& x) //G37. Тестирование на ноль. Возвращает true, если вектор нулевой 
{
	int i;
	DO(i, 1, x.m)if (x(i) != 0)return false;
	return true;
}

bool TestZero(Matr& A) //G38. Тестирование на ноль. Возвращает true, если матрица нулевая
{
	int i, j; DO(i, 1, A->m)DO(j, 1, A->n)if (A(i, j) != 0)return false;
	return true;
}

double det(Matr& A) //G39. Определитель матрицы
{
	int i, j, k, n = A.n;
	double mx, d = 1, t;
	Matr P(n, n);
	DO(i, 1, n)DO(j, 1, n)P(i, j) = A(i, j);
	DO(k, 1, n)
	{
		mx = 0; DO(i, k, n) { t = P(i, k); if (fabs(t) > fabs(mx)) { mx = t; j = i; } }
		if (mx == 0) { d = 0; return d; }
		if (j != k) { d = -d; DO(i, k, n) { t = P(j, i); P(j, i) = P(k, i); P(k, i) = t; } }
		DO(i, k + 1, n) { t = P(i, k) / mx; DO(j, k + 1, n)P(i, j) -= t * P(k, j); }
		d *= P(k, k);
	}
	return d;
}

Vect& Adding(Vect& x, Vect& y) //G40. Прямое сложение векторов
{
	int mx = x->m, my = y->m, n = mx + my;
	VSTACK(n, xy)
	xy(x, 1); xy(y, mx + 1);
	return xy;
}

//													  T
//G41. Билинейная (квадратичная при x=y) форма Forme=x Py 
double Forme(Vect& x, Matr& P, Vect& y)
{
	int i, j, m = P.m, n = P.n; double s = 0;
	DO(i, 1, m)DO(j, 1, n)s += x(i) * P(i, j) * y(j);
	return s;
}

//						         T -1      T
//G42. Квадратичная форма Forme=x P  x (P=P ) 
double Forme(Vect& x, Matr& P)
{
	int i, j, k, n = x.m; Matr Q(n + 1, n + 1); double mu, a;
	DO(i, 1, n) { Q(i, n + 1) = x(i); Q(n + 1, i) = x(i); DO(j, 1, n)Q(i, j) = P(i, j); }
	Q(n + 1, n + 1) = 0;
	DO(i, 1, n)
		if (Q(i, i) == 0)
		{
			DO(j, i + 1, n)
				if (Q(j, i) != 0)DO(k, 1, n + 1) { a = Q(i, k); Q(i, k) = Q(j, k); Q(j, k) = a; }
		}
		else
			DO(j, i + 1, n + 1)
		{
			mu = Q(j, i) / Q(i, i); DO(k, i + 1, n + 1)Q(j, k) = Q(j, k) - mu * Q(i, k);
		}
	return -Q(n + 1, n + 1);
}

//G43. Интегрирование системы дифференциальных уравнений
//методом Рунге-Кутта с постоянным шагом
//Входные параметры:                                                      
//t			- начальное значение независимой переменной;	                              
//T			- конечное значение независимой переменной;	                              
//x=x(t)	- начальный вектор интегрируемых параметров	                      
//F			- вектор-функция правых частей системы уравнений f(t,x)
//Выходные параметры:	                                                  
//x=x(T)	- конечный вектор интегрируемых параметров             	          
void Runge(double t, double T, Vect& x, Vect& F(double t, Vect& x))
{
	int j, k, n = x->m;
	double t0 = t, dt = T - t;
	Vect x0(n), z(n), w(n), a(5);
	a(1) = 0.5 * dt; a(2) = a(1); a(5) = a(1); a(3) = dt; a(4) = dt; x0 = x; w = x;
	DO(j, 1, 4)
	{
		z = F(t, w); t = t0 + a(j);
		DO(k, 1, n) { x(k) += a(j + 1) * z(k) / 3; w(k) = x0(k) + a(j) * z(k); }
	}
}

//G44. Интегрирование системы дифференциальных уравнений
//методом Рунге-Кутта с постоянным шагом
//Входные параметры                                                      
//F			- вектор-функция правых частей системы уравнений f(t,x)
//dt		- шаг интегрирования;	                              
//t			- начальное значение независимой переменной;	                              
//x=x(t)	- начальный вектор интегрируемых параметров	                      
//Выходные параметры:	                                                  
//x=x(t+dt)	- конечный вектор интегрируемых параметров             	          
void Runge(Vect& F(double t, Vect& x), double t, double dt, Vect& x)
{
	int j, k, n = x->m;
	double t0 = t;
	Vect x0(n), z(n), w(n), a(5);
	a(1) = 0.5 * dt; a(2) = a(1); a(5) = a(1); a(3) = dt; a(4) = dt; x0 = x; w = x;
	DO(j, 1, 4)
	{
		z = F(t, w); t = t0 + a(j);
		DO(k, 1, n) { x(k) += a(j + 1) * z(k) / 3; w(k) = x0(k) + a(j) * z(k); }
	}
}

//G45. Интегрирование системы дифференциальных уравнений
//методом Рунге-Кутта с автоматическим выбором шага
//F		- вектор-функция правых частей системы уравнений f(t,x);
//eps	- точность интегрирования;	                              
//t		- начальное значение независимой переменной;	                              
//T		- конечное значение независимой переменной;	                              
//x		- начальный вектор интегрируемых параметров x(t),
//		  на выходе x - конечный вектор интегрируемых
//		  параметров x(T) 
double Runge(Vect& F(double t, Vect& x), double eps, double t, double T, Vect& x)
{
	int k, n = x->m; static int ss;
	static bool prim = 1; bool out;
	double t1, t2, t3, h; static double hs = T - t;
	Vect x1(n), x2(n), x3(n);
	if (prim) { h = T - t; ss = 0; }
	else h = hs; out = 0;
lab1:
	if ((t + 2.01 * h - T > 0) == (h > 0)) { hs = h; out = 1; h = 0.5 * (T - t); }
	t1 = t + 2 * h; x1 = x; Runge(F, t, 2 * h, x1);
lab2:
	t2 = t + h; x2 = x; Runge(F, t, h, x2);
	t3 = t2 + h; x3 = x2; Runge(F, t2, h, x3);
	DO(k, 1, n)if (Comp(x1(k), x3(k), eps) > eps)goto lab3;
	t = t3; if (out)goto fin; x = x3; if (ss == 5) { ss = 0; h = 2 * h; }
	ss++; goto lab1;
lab3:
	h *= 0.5; out = 0; t1 = t2; x1 = x2; goto lab2;
fin:
	prim = 0; x = x3;
	return hs;
}

//G46. Интегрирование матричной системы дифференциальных
//уравнений методом Рунге-Кутта с постоянным шагом
//Входные параметры                                                      
//F			- матричная функция правых частей
//			  системы уравнений f(t,X)
//dt		- шаг интегрирования;	                              
//t			- начальное значение независимой переменной;	                              
//X=X(t)	- начальная матрица интегрируемых параметров	                      
//Выходные параметры:	                                                  
//X=X(t+dt)	- конечная матрица интегрируемых параметров             	          
void Runge(Matr& F(double t, Matr& X), double t, double dt, Matr& X)
{
	int j, k, n = X->m, nn = n * n;
	double t0 = t;
	Vect x(nn), x0(nn), z(nn), w(nn), a(5);
	Matr Z(n, n), W(n, n);
	x = X->MatVecStr(); //Упаковка матрицы X в вектор x
	a(1) = 0.5 * dt; a(2) = a(1); a(5) = a(1); a(3) = dt; a(4) = dt; x0 = x; w = x;
	DO(j, 1, 4)
	{
		W->VecMatStr(w); //Упаковка вектора w в матрицу W
		Z = F(t, W); z = Z->MatVecStr(); t = t0 + a(j);
		DO(k, 1, nn) { x(k) += a(j + 1) * z(k) / 3; w(k) = x0(k) + a(j) * z(k); }
	}
	X->VecMatStr(x);
}

//G47. Интегрирование матричной системы дифференциальных уравнений
//методом Рунге-Кутта с автоматическим выбором шага
//Входные параметры:	                                                  
//F			- матричная функция правых частей системы уравнений f(t,X);
//eps		- точность интегрирования;	                              
//t			- начальное значение независимой переменной;	                              
//T			- конечное значение независимой переменной;	                              
//X			- начальная матрица интегрируемых параметров X(t),
//Выходные параметры:	                                                  
//X=X(T)	- конечная матрица интегрируемых параметров             	          
double Runge(Matr& F(double t, Matr& X), double eps, double t, double T, Matr& X)
{
	int k, n = X->m, nn = n * n;
	static int ss; static bool prim = 1; bool out;
	double t1, t2, t3, h; static double hs;
	Vect x(nn), x1(nn), x2(nn), x3(nn);
	Matr X1(n, n), X2(n, n), X3(n, n);
	x = X->MatVecStr();
	if (prim) { h = T - t; ss = 0; }
	else h = hs; out = 0;
lab1:
	if ((t + 2.01 * h - T > 0) == (h > 0)) { hs = h; out = 1; h = 0.5 * (T - t); }
	t1 = t + 2 * h; x1 = x; X1->VecMatStr(x1);
	Runge(F, t, 2 * h, X1); x1 = X1->MatVecStr();
lab2:
	t2 = t + h; x2 = x; X2->VecMatStr(x2);
	Runge(F, t, h, X2); x2 = X2->MatVecStr();
	t3 = t2 + h; x3 = x2; X3->VecMatStr(x3);
	Runge(F, t2, h, X3); x3 = X3->MatVecStr();
	DO(k, 1, nn)if (Comp(x1(k), x3(k), eps) > eps)goto lab3;
	t = t3; if (out)goto fin; x = x3; if (ss == 5) { ss = 0; h = 2 * h; }
	ss++; goto lab1;
lab3:
	h *= 0.5; out = 0; t1 = t2; x1 = x2; goto lab2;
fin:
	prim = 0; x = x3; X->VecMatStr(x3);
	return hs;
}

//G48. Определяет абсолютное значение разности мантисс a и b после
//того как порядки этих величин выравнены до наибольшего
//порядка параметров a,b и c
double Comp(double a, double b, double c)
{
	int ae, be, ce; ae = Expon(a); be = Expon(b); ce = Expon(c);
	if (ae < be)ae = be; if (ae < ce)ae = ce;
	return fabs(a - b) / pow(10, ae);
}

//G49. Эта функция вычисляет показатель десятичного порядка
//нормализованного числа x с плавающей точкой
//int Expon(double x){return x==0? -300: floor(0.4342944819*log(fabs(x))+1);}
int Expon(double x)
{
	return x == 0 ? DBL_MIN_10_EXP : floor(0.4342944819 * log(fabs(x)) + 1);
}

//G50. Множители интерполяционного полинома
//Лагранжа при узловых значениях функции f().
//t - точка (неузловая в общем случае), в которой требуется
//вычислить функцию, T - вектор значений аргумента размерности n.
//					T
//Примечание: f(t)=c [f(T(1)),f(T(2)),...,f(T(n))]
Vect Lagrange(Vect& T, double t)
{
	int i, j, n = T->m;
	Vect c(n);
	DO(i, 1, n)
	{
		c(i) = 1;
		DO(j, 1, n)if (i != j)c(i) *= (t - T(j)) / (T(i) - T(j));
	}
	STACKh(n, T)
	*(T->h) = c;
	return *(T->h);
}

//G51. Интерполяция полиномами Лагранжа степени L
//для равноотстоящих узловых значений аргумента
//Вход:
//L - наибольшая степень "старшего" из полиномов Лагранжа;
//t - значение аргумента, для которого
//    определяется значение функции;
//T - массив значений аргумента;
//F - массив значений функции f()
//Возвращаемое значение: f(t)
double Lagrange(int L, double t, Vect& T, Vect& F)
{
	int n = T.m, m = (t - T(1)) / (T(2) - T(1)) + 2, dm1 = L / 2, dm2 = (L - 1) / 2;
	Vect p(L), q(L), c(L);
	if (m < L - 1) { p = T; q = F; }
	else if (m > n - L + 1) { p = T(n - L + 1, n); q = F(n - L + 1, n); }
	else { p = T(m - dm1, m + dm2); q = F(m - dm1, m + dm2); }
	return Lagrange(p, t) * q;
}

int Cnm(int n, int m) //G52. Количество сочетаний из n по m
{
	int m0, i, j, k, N; m0 = m;
	if (m <= n - m)k = m; else k = n - m;
	m = n - k; if (k == 0)j = 1; else j = m + 1;
	DO(i, 2, k)j = (double)(m + i) / i * j;
	N = j; m = m0;
	return N;
}

//G53. Генератор сочетаний.
//										    m
//Определяет все сочетания из n по m после C  обращений к Comb.
//										    n
//Перед первым обращением следует задать prim=true, 
//В процессе промежуточных обращений prim будет получать значение false.
//После получения всех сочетаний prim снова примет значение true
void Comb(bool& prim, int n, int m, Vect& c)
{
	int i, j;
	if (prim) { prim = 0; DO(j, 1, m)c(j) = j; return; }
	if (c(m) < n) { c(m)++; return; }
	dO(j, m - 1, 1, -1)
		if (c(j) < n - m + j)
		{
			c(j)++; DO(i, j + 1, m)c(i) = c(j) + i - j;
			return;
		}
	prim = 1;
}

//G54. Генератор сочетаний.
//Определяет все сочетания из n по m.
//Каждое сочетание занимает одну из строк
//						  m
//возвращаемой матрицы A(C ,m)
//						  n
Matr& Comb(int n, int m)
{
	bool prim = 1; int k, N = Cnm(n, m);
	Vect c(m);
	Matr C(N, m);
	//E->H = new Matr(N, m); E->NumbMatr++;
	STACKH(N, m, E)
		DO(k, 1, N) { Comb(prim, n, m, c); C(k, c); }
	*(E->H) = C;
	return *(E->H);
}
void emp() //G55. Пропуск строки на экране 
{
	cout << "\n";
}

void emp(int n) //G56. Пропуск n строк на экране
{
	int i;
	DO(i, 1, n)cout << "\n";
}

//G57. Функция SubMatr возвращает выделенную из заданной
//матрицы квадратную подматрицу полного ранга.
//В заданной матрице меняется расположение строк
//и столбцов так, чтобы в левом верхнем
//углу разместилась возвращаемая матрица
Matr& SubMatr(Matr& B)
{
	int i, j, m = B->m, n = B->n, rang, k = 1, si, ci;
	Vect s(m), c(n);
	rang = RangMatr(B, s, c, 1e-12);
	if (rang < m)
		DO(i, 1, rang)
	{
		si = int(s(k));
		if (i != si) B->TrStr(i, si);
		else { k++; continue; }
	}
	k = 1;
	if (rang < n)
		DO(i, 1, rang)
	{
		ci = int(c(k));
		if (i != ci) B->TrCol(i, ci);
		else { k++; continue; }
	}
	return B;
}

char* Char(Str& x) //G58. Преобразование Str->char
{
	int i = 0, m = x->m;
	x->s = Line(m + 1); x->numbStr++;
	DO(i, 1, m) x->s[i - 1] = x(i);
	return  x->s;
}

void Char(Str& x, char* s) //G59. Преобразование Str->char 
{
	int i = 0, m = x->m;
	DO(i, 1, m)s[i - 1] = x(i);
}

Vect& VecProd(Vect& x, Vect& y) //G60. Векторное произведение векторов
{
	Vect z(3);
	z(1) = x(2) * y(3) - x(3) * y(2);
	z(2) = x(3) * y(1) - x(1) * y(3);
	z(3) = x(1) * y(2) - x(2) * y(1);
	STACKh(3, x)
	*(x->h) = z;
	return z;
}

//G61. Переход от вектора x в J2000 к оскулирующим
//элементам орбиты xos=[e,om,a,i,Om,u]
//(0<= i <=pi, 0<= om <=2*pi, 0<= Om <=2*pi, 0<= u <=2*pi) 
Vect& AbsOsc(Vect& x)
{
	double r, v, rr, vv, rv, mc, c2, cT, sOm, cOm, si, ci;
	Vect c(3), X(3), V(3), xos(6);
	X = x; V = x(4, 6); c = VecProd(X, V); r = !x(1, 3); v = !x(4, 6);
	c2 = c * c; mc = sqrt(c2); rr = X * X; vv = V * V;
	xos(4) = acos(c(3) / mc); si = sin(xos(4)); ci = cos(xos(4));
	xos(5) = atan2(c(1), -c(2)); sOm = sin(xos(5)); cOm = cos(xos(5));
	xos(6) = atan2(x(3) / si, x(1) * cOm + x(2) * sOm);
	xos(3) = c2 / mu_; xos(1) = 1 - xos(3) / mu_ * (2 * mu_ / r - vv);
	if (xos(1) >= 0)xos(1) = sqrt(xos(1)); else xos(1) = 0;
	rv = X * V; cT = sqrt(rr * vv - rv * rv) / r / v;
	if (xos(1) > 1e-6)xos(2) = xos(6) - atan2(rv * mc / mu_, c2 / mu_ - r);
	xos(3) /= (1 - xos(1) * xos(1));
	if (xos(2) < 0)xos(2) += pi2; if (xos(5) < 0)xos(5) += pi2;
	STACKh(6, x)
	*(x->h) = xos;
	return *(x->h);
}

//G62. Высота над поверхностью референц-эллипсоида.
//x(3) - гринвичские координаты
double Altitude(Vect& x)
{
	double acm_ = 0.003350153712621, e, fi, r, Re, del, eps, H, sf, cf, cf2, a;
	e = 2 * acm_ - acm_ * acm_; r = !x; sf = x(3) / r; 
	fi = asin(sf); cf = cos(fi); cf2 = cf * cf;
	a = 1 - e * cf2; Re = aE_ * sqrt((1 - e) / a); 
	del = e * sf * cf / a; eps = Re * del / r;
	H = (r - Re) * (1 - 0.5 * eps * del);
	return H;
}

//G63. Угол между плоскостью орбиты и направлением на Солнца из центра Земли.
//t - время по гринвичу;
//i - наклонение; Om - долгота восходящего узла в J2000
double AngleSunOrb(double ud, double t, double i, double Om)
{
	double udt = ud + t / 86400, alf; 
	Vect rS(3); 
	Sun(udt, rS);
	alf = VectKinMom(i, Om) % rS;
	alf = alf < pi_2 ? pi_2 - alf : alf - pi_2;
	return alf;
}

//G64. Угол между плоскостью орбиты и направлением на Солнца из центра Земли.
//t - время по гринвичу; x - ВС в J2000
double AngleSunOrb(double ud, double t, Vect& x)
{
	double udt = ud + t / 86400, alf; 
	Vect rS(3); 
	Sun(udt, rS);
	alf = VecProd(x(1, 3), x(4, 6)) % rS;
	alf = alf < pi_2 ? pi_2 - alf : alf - pi_2;
	return alf;
}

//G65. Угол между плоскостью орбиты и направлением на Солнца из центра Земли.
//t - время по гринвичу; x - ВС в ГСК
double AngleSunOrb(double ud, Vect& x, double t)
{
	double udt = ud + t / 86400, alf; Vect rS(3); 
	Sun(udt, rS);
	rS = Mat2000Gr(ud, t) * rS; alf = VecProd(x(1, 3), x(4, 6)) % rS;
	alf = alf < pi_2 ? pi_2 - alf : alf - pi_2;
	return alf;
}

//G66. Время перелета от восходящего узла до перигея 
double AnodePrg(double T, double e, double om)
{
	double q, dt, com = cos(om); q = (sqrt(1 - e * e)) * sin(om) / (1 + e * com);
	dt = T / pi2 * (Asin(q, sign(e + com))) - e * q; return dt;
}

//G67. Решение уравнения Кеплера
//(определение истинной аномалии по заданному времени)
double Anom(double mu, double eps, double tau, double e, double a, double t)
{
	double T, E0, E, cE, c; T = sqrt(mu / a) / a * (t - tau);
	if (e == 0)return atan2(sin(T), cos(T)); E0 = T + e * sin(T);
Lab:
	E = T + e * sin(E0); if (fabs(E - E0) > eps) { E0 = E; goto Lab; }
	cE = cos(E); c = 1 - e * cE; return atan2(sin(E) / c * sqrt(1 - e * e), (cE - e) / c);
}

//G68. Определение аргумента широты
double ArgWid(double eps, Vect& xos, double t)
{
	double u; u = xos(4) + Anom(mu_, eps, xos(1), xos(2), xos(3), t); return u;
}

//G69. Определение местного азимута трассы.
//Входные параметры:
//ud	- юлианскач дата;
//t		- гринвичское времяб с; 
//x(6)	- фазовый вектор КА в J2000.
//Выходные параметры:
//sinA	- синус местного азимута;
//sinA	- косинус местного азимута;
void Azim(double ud, double t, Vect& x, double& sinA, double& cosA)
{
	Vect b(3), c(3), rg(3); 
	Matr B(3, 3); 
	B = Mat2000Gr(ud, t);
	rg = B * x(1, 3); b(1) = -rg(2) / rg(1); 
	b(2) = 1; c = B * VecProd(x(1, 3), x(4, 6));
	cosA = b * c / !b / !c; sinA = sqrt(1 - cosA * cosA);
}

//G70. Круговой арксинус.
//Определяет угол в пределах от 0 радиан до 2*pi радиан.
//х - синус угла; у - косинус угла или сигнатура косинуса угла  
double Asin(double x, double y)
{
	if (x >= 0 && y >= 0)return asin(x);
	if (x >= 0 && y < 0)return pi - asin(x);
	if (x < 0 && y <= 0)return pi - asin(x); //"-" потому что asin(x) - нечетная
	if (x < 0 && y>0)return pi2 + asin(x); //"+" потому что asin(x) - нечетная
	return 0;
}

//G71. Круговой арксинус.
//Определяет угол в пределах от 0 радиан до 2*pi радиан.
//y - противолежащий, x - прилежащий катеты
double Arcsin(double y, double x)
{
	if (x == 0 && y == 0)
	{
		Ecrir(cout, "\n dArcsin---> error: x=", x, " y=", y, "\n"); END
	}
	double as = asin(y / sqrt(x * x + y * y));
	if (x >= 0 && y >= 0)return as;
	if (x <= 0 && y >= 0)return pi - as;
	if (x <= 0 && y <= 0)return pi - as;  //"-" потому что asin(x) - нечетная
	if (x >= 0 && y <= 0)return pi2 + as; //"+" потому что asin(x) - нечетная
	return 0;
}

//G72. Производная от кругового арксинуса.
//х - синус угла; у - косинус угла или сигнатура косинуса угла  
double dAsin(double x, double y)
{
	if (fabs(x) == 1) { Ecrir(cout, "\n dAsin---> error: x=", x, "\n"); END }
	double dA = 1 / sqrt(1 - x * x);
	if (x >= 0 && y >= 0 || x < 0 && y>0)return dA;
	if (x >= 0 && y < 0 || x < 0 && y <= 0)return -dA;
	return 0;
}

//G73. Производная от кругового арксинуса.
//y - противолежащий, x - прилежащий катеты
double dArcsin(double y, double x)
{
	if (x == 0) { Ecrir(cout, "\n dArcsin---> error: x=", x, " y=", y, "\n"); END }
	double dA = sqrt(x * x + y * y) / fabs(x);
	if (x >= 0 && y >= 0 || x >= 0 && y <= 0)return dA;
	if (x <= 0 && y >= 0 || x <= 0 && y <= 0)return -dA;
	return 0;
}

//G74. Круговой арктангенс.
//y - противолежащий, x - прилежащий катеты.
//Определяет угол в пределах от 0 радиан до 2*pi радиан
double Arctg(double y, double x)
{
	if (x == 0) { Ecrir(cout, "\n Arctg---> error: x=", x, "\n"); END }
	double at = atan(y / x);
	if (x >= 0 && y >= 0)return at;
	if (x <= 0 && y >= 0)return pi + at;
	if (x <= 0 && y <= 0)return pi + at;
	if (x >= 0 && y <= 0)return pi2 + at;
	return 0;
}

//G75. Производная от кругового арктангенса.
//y - противолежащий, x - прилежащий катеты.
double dArctg(double y, double x) { double x2 = x * x; return x2 / (x2 + y * y); }

//G76. К динамической модели атмосферы.
//Среднесуточный F107 и средневзвешенный F81 индексы солнечной активности.
//Индекс геомагнитной возмущенности Kp. Считывание из файла dateAtm.txt
void ISA_IGI(int d, int m, int y, double& F107, double& F81, double& Kp)
{
	int y0, m0, d0; fstream dateAtm("d:\\Arsenal\\Data\\dateAtm.txt", ios::in);
	do Read(dateAtm, y0, m0, d0, F107, F81, Kp); while (!(y0 == y && m0 == m && d0 == d));
}

//G77. К динамической модели атмосферы.
//Множитель A(n) (n - число суток от начала года)
double FuncA(int d, int m, int y)
{
	int n = Ulian(d, m, y) - Ulian(1, 1, y) + 1, s;
	double A, A0 = -2.53418e-2, A1 = -2.44075e-3, A2 = 3.08389e-6,
		A3 = 2.90115e-6, A4 = -4.99606e-8, A5 = 3.36327e-10,
		A6 = -1.0966e-12, A7 = 1.73227e-15, A8 = -1.06271e-18;
	s = n; A = A0 + A1 * s; s *= n; A += A2 * s; s *= n; A += A3 * s; s *= n; A += A4 * s;
	s *= n; A += A5 * s; s *= n; A += A6 * s; s *= n; A += A7 * s;  s *= n; A += A8 * s;
	return A;
}

//G78. К динамической модели атмосферы.
//Определение фиксированного ИСА F0 по средневзвешенному F81
int SelectF0(double F81)
{
	int i, j = 1, F0[7] = { 75,100,125,150,175,200,250 };
	double dF0 = fabs(F81 - F0[0]), dF;
	DO(i, 2, 7)
	{
		dF = fabs(F81 - F0[i - 1]);
		if (dF < dF0) { j = i; dF0 = dF; }
	}
	return F0[j - 1];
}

Vect& Soleil(double ud, double tg) //G79. Определение координат Солнца в J2000
{
	int ns; double d, t, t2, e, r, v, h, rr;
	Vect xs(3);
	ns = int(tg / 86400); ud = ud + ns; tg = tg - ns * 86400;
	d = ud + tg / 86400 - 2415020; t = d / 36525; t2 = t * t;
	e = 0.1675104e-1 - (0.418e-4 + 0.126e-6 * t) * t;
	r = 6.25658378411 + 1.72019697677e-2 * d - 2.61799387799e-6 * t2;
	v = 4.90822940869 + 8.21498553644e-7 * d + 7.90634151156e-6 * t2;
	r = r + 2 * e * sin(r) + 1.25 * e * e * sin(2 * r);
	h = 4.09319747446e-1 -
		(2.27110968916e-4 + (2.86233997327e-8 - 8.77900613756e-9 * t) * t) * t +
		4.46513400244e-5 * cos(4.52360151485 -
		(33.7571462465 - (3.62640633471e-5 + 3.87850944887e-8 * t) * t) * t);
	rr = 149600034.408 * (1 - e * e) / (1 + e * cos(r)); v = v + r; r = sin(v);
	xs(1) = rr * cos(v); xs(2) = rr * r * cos(h); xs(3) = rr * r * sin(h);
	xs = xs * 1000; xs = !Reduct(ud, tg) * xs;
	STACKh(3, ::e)
	*(::e->h) = xs;
	return *(::e->h);
}

//G80. Определение местного солнечного времени [сек]
//по заданной юлианской дате ud, гривичскому
//времени t [сек] и положению в J2000 x(3) [м]
double TimeLocal(double ud, double t, Vect& x)
{
	double L, tL; Vect xg(3); xg = Mat2000Gr(ud, t) * x;
	L = Arcsin(xg(2), xg(1)); //Географическая долгота
	tL = t + 43200 / pi * L; return tL;
}

//G81. Определение гринвичского времени [сек]
//по заданной юлианской дате ud, местному
//солнечному времени t [сек]
//и положению в J2000 x(3) [м]
double TimeGrinv(double ud, double t, Vect& x)
{
	double L, tG; Vect xg(3); xg = Mat2000Gr(ud, t) * x;
	L = Arcsin(xg(2), xg(1)); //Географическая долгота
	tG = t - 43200 / pi * L; return tG;
}

//G82. Перевод юлианской даты в календарную дату.                 
//Входные параметры: ud - юлианская дата.
//Выходные параметры: d - день; m - месяц; y - год
void Date(double ud, int& d, int& m, int& y)
{
	double j = ud - 1721119 + 0.5;
	y = floor((4 * j - 1) / 146097); d = floor((4 * j - 1 - 146097 * y) / 4);
	j = floor((4 * d + 3) / 1461); d = floor((4 * d + 7 - 1461 * j) / 4);
	m = floor((5 * d - 3) / 153); d = floor((5 * d + 2 - 153 * m) / 5); y = 100 * y + j;
	if (m < 10)m += 3; else { m -= 9; y++; }
}

//G83. Перевод юлианской даты в календарную дату в виде вектора CD=[d,m,y]
Vect& Date(double ud)
{
	int d, m, y; double j = ud - 1721119 + 0.5; 
	Vect CD(3);
	y = floor((4 * j - 1) / 146097); d = floor((4 * j - 1 - 146097 * y) / 4);
	j = floor((4 * d + 3) / 1461); d = floor((4 * d + 7 - 1461 * j) / 4);
	m = floor((5 * d - 3) / 153); d = floor((5 * d + 2 - 153 * m) / 5); y = 100 * y + j;
	if (m < 10)m += 3; else { m -= 9; y++; }
	CD(1) = d; CD(2) = m; CD(3) = y;
	STACKh(3, e)
	*(e->h) = CD;
	return *(e->h);
}

//G84. Распаковка даты из формата 00.00.00 в вектор
Vect& DMY(char* dat)
{
	char sd[3], sm[3], sy[5];
	sd[0] = dat[0]; sd[1] = dat[1];
	Vect v(3);
	v(1) = atoi(sd);
	sm[0] = dat[3]; sm[1] = dat[4];
	v(2) = atoi(sm);
	sy[0] = dat[6]; sy[1] = dat[7]; sy[2] = dat[8];
	sy[3] = dat[9];
	v(3) = atoi(sy);
	STACKh(3, e)
	*(e->h) = v;
	return *(e->h);
}

//G85. Перевод юлианской даты в строку вида 00/00/0000
//(день/месяц/год)
Str& DMY(double ud)
{
	int d, m, y; Str s1(2), s2(2), s3(4);
	Str S(10);
	Date(ud, d, m, y);	s1 = Liter(d, 2); 
	s2 = Liter(m, 2); s3 = Liter(y, 4);
	S = s1 + (char*)"/" + s2 + (char*)"/" + s3;
	*(eS->h) = S;
	return *(eS->h);
}

//G86. Распаковка времени из формата 00:00:00.0000 в вектор
Vect& HMS(char* tim)
{
	char sh[3], sm[3], ss[8], ** end = NULL;
	Vect v(3);
	sh[0] = tim[0]; sh[1] = tim[1]; 
	v(1) = atoi(sh);
	sm[0] = tim[3]; sm[1] = tim[4]; 
	v(2) = atoi(sm);
	ss[0] = tim[6]; ss[1] = tim[7]; ss[2] = tim[8]; ss[3] = tim[9];
	ss[4] = tim[10]; ss[5] = tim[11]; ss[6] = tim[12]; ss[7] = tim[13];
	v(3) = strtod(ss, end);
	STACKh(3, e)
	*(e->h) = v;
	return *(e->h);
}

//G87. Перевод секунд в строку вида 00:00:00.0000
//(часы:минуты:секунды)
Str& HMS(double t)
{
	int d = t / 86400, h, m, s, p, v1, v2;
	double q;
	Str s1(2), s2(2), s3(2), s4(3), S(13);
	v1 = 86400 * d; h = (t - v1) / 3600; s1 = Liter(h, 2);
	v2 = v1 + 3600 * h; m = (t - v2) / 60; s2 = Liter(m, 2);
	s = t - v2 - 60 * m; s3 = Liter(s, 2);
	q = t - v2 - 60 * m - s + 1e-12; p = q * 1000; s4 = Liter(p, 3);
	if (q < 0.1)s4(1) = '0'; if (q < 0.01)s4(2) = '0'; if (q < 0.001)s4(3) = '0';
	S = s1 + (char*)":" + s2 + (char*)":" + s3 + (char*)"." + s4;
	STACKSh(13, eS)
	*(eS->h) = S;
	return *(eS->h);
}

//G88. Перевод юлианской даты и времени в секундах
//в строку вида 00/00/0000 00:00:00.0000
//(день/месяц/год часы:минуты:секунды)
Str& DHMS(double ud, double t)
{
	Str S1(10), S2(13), S(24);
	S1 = DMY(ud); S2 = HMS(t); 
	S = S1 + (char*)" " + S2;
	STACKSh(24, eS)
	*(eS->h) = S;
	return *(eS->h);
}

//G89. Переход от геодезических координат (g=||B,L,H||) к гринвичским
Vect& GeodGrinv(Vect& g)
{
	double R, sB, cB, sL, cL;
	Vect x(3);
	sB = sin(g(1)); cB = cos(g(1)); sL = sin(g(2)); cL = cos(g(2));
	R = aE_ / sqrt(1 - (2 * acm_ - acm_ * acm_) * sB * sB);
	x(1) = (R + g(3)) * cB * cL; x(2) = (R + g(3)) * cB * sL;
	x(3) = (R * (1 - acm_) * (1 - acm_) + g(3)) * sB;
	STACKh(3, g)
	*(g->h) = x;
	return *(g->h);
}

//G90. Переход от гринвичских координат к геодезическим (g=||B,L,H||)
Vect& GrinvGeod(Vect& x)
{
	double zr, fi, cfi, r, R, e, ecfi, eps, dlt;
	Vect g(3);
	r = !x; zr = x(3) / r; fi = asin(zr); cfi = cos(fi); e = 2 * acm_ - acm_ * acm_;
	ecfi = 1 - e * cfi * cfi; dlt = e * zr * cfi / ecfi;
	R = aE_ * sqrt(1 - e) / sqrt(ecfi); eps = R * dlt / r;
	g(1) = fi + eps; g(2) = atan(x(2) / x(1)); g(3) = (r - R) * (1 - 0.5 * eps * dlt);
	STACKh(3, x)
	*(x->h) = g;
	return *(x->h);
}

//G91. Поправка к минимальному углу возвышения над горизонтом
//(условию радиовидимости), учитывающей эллипсоидальность
//Земли и высоту H над поверхностью земного эллипсоида.
//fip - географическая (геоцентрическая) широта пункта
double dGamma(double fip, double H)
{
	double fi, G1, G2, sf, cf, div, dG;
	fi = atan(tan(fip) / (1 - acm_) / (1 - acm_));
	sf = sin(fi); sf *= sf; cf = cos(fi); cf *= cf;
	div = sqrt(1 - (2 * acm_ - acm_ * acm_) * sf);
	G1 = aE_ / div + H; G2 = aE_ * (1 - acm_) * (1 - acm_) / div + H;
	dG = acos((G1 * cf + G2 * sf) / rE_);
	return dG;
}

//G92. Поправка к минимальному углу возвышения над горизонтом
//(условию радиовидимости), учитывающей эллипсоидальность
//Земли и высоту H над поверхностью земного эллипсоида.
//B - геодезическая широта пункта
void dGamma(double B, double H, double& dG)
{
	double G1, G2, sB, cB, div;
	sB = sin(B); sB *= sB; cB = cos(B); cB *= cB;
	div = sqrt(1 - (2 * acm_ - acm_ * acm_) * sB);
	G1 = aE_ / div + H; G2 = aE_ * (1 - acm_) * (1 - acm_) / div + H;
	dG = acos((G1 * cB + G2 * sB) / rE_);
}

//G93. Определение времени входа и выхода из зоны радиовидимости наземного пункта 
//Входные параметры:	                                                      
//GamMin- ограничение на угол возвышения назенной антенны;                 
//Lrn - географическая долгота восходящего узла;                           
//xos	- массив элементов орбиты, в последовательности:	               
//xos(1) - момент прохождения перигея;	                                  
//xos(2) - эксцентриситет;	                                              
//xos(3) - большая полуось;	                                              
//xos(4) - аргумент перигея;	                                              
//xos(5) - наклонение;	                                                  
//xos(6) - долгота восходящего узла(не используется);	                  
//q		- массив геодезических координат наземного пункта:        	      
//q(1)	- геодезическая широта;                                           
//q(2)	- геодезическая долгота;                                          
//q(3)	- высота над поверхностью земного эллипсоида.                     
//Выходные параметры:	                                                      
//T0 - время входа в зону радиовидимости от момента прохождения перигея;   
//Tf - время выхода из зоны радиовидимости от момента прохождения перигея; 
void TimeRadio(double F(Vect& G, double e), Vect& xos, Vect& q, double eps,
			   double GamMin, double Lrn, double& T0, double& Tf)
{
	int j; Vect G(16);
	double per, sq, som, com, somg, comg, si, ci, sf, rp, hy, hh, ff,
		   e0, e1, de1, g1, dLam, u1, u2, uu1, uu2, e, a, Lam,
		   ps, px, py, pz, qx, qy, qz, k1, g2, rp2, cf, rrc, rrs;
	Tf = 0; T0 = Tf; e = xos(2); a = xos(3);
	Lam = q(2); si = pi - xos(4) / 2; per = sqrt(a / mu_) * a;
	sq = sqrt((1 - e) / (1 + e)); si = 2 * atan(tan(si) * sq);
	ps = per * omE_; si = (si - e * sin(si)) * ps + Lrn;
	som = sin(xos(4)); com = cos(xos(4)); somg = sin(si); comg = cos(si);
	si = sin(xos(5)); ci = cos(xos(5)); k1 = (1 - acm_); k1 = k1 * k1;
	cf = atan(tan(q(1)) / k1); sf = sin(cf); cf = cos(cf);
	px = sf * sf; py = cf * cf;	g2 = aE_ / (py + k1 * px); g1 = g2 + q(3); g2 = g2 * k1 + q(3);
	rp = sqrt(g1 * g1 * py + g2 * g2 * px);	hy = (g1 * py + g2 * px) / rp;
	if (hy > 1)hy = 1; hy = acos(hy) + GamMin; g2 = g2 / g1 * sf; g1 = g1 * a;
	k1 = sqrt(1 - e * e); pz = som * si;
	px = com * comg - som * somg * ci; py = com * somg + som * comg * ci;
	qz = k1 * com * si; qx = -k1 * (som * comg + com * somg * ci);
	qy = k1 * (-som * somg + com * comg * ci);
	rrc = rp * cos(hy); rrc = rrc * rrc;
	rrs = rp * sin(hy) / g1; rp2 = rrc / g1; k1 = a * e;
	hh = 6 / per; u1 = sf / si; if (u1 > 1)u1 = 1;
	u1 = asin(u1); u2 = pi - u1;
	u1 = (u1 - xos(4)) / 2; u2 = (u2 - xos(4)) / 2;
	uu1 = 2 * atan(tan(u1) * sq); uu2 = 2 * atan(tan(u2) * sq);
	e0 = uu1 - pi; ff = uu1 + pi;
	if (Lam > Lrn) { dLam = Lam - Lrn; if (dLam > pi)dLam = pi2 + Lrn - Lam; }
	else { dLam = Lrn - Lam; if (dLam > pi)dLam = pi2 + Lam - Lrn; }
	G(1) = e; G(2) = a; G(3) = Lam; G(4) = ps; G(5) = px; G(6) = py; G(7) = pz; G(8) = qx;
	G(9) = qy; G(10) = qz; G(11) = k1; G(12) = g2; G(13) = rp2; G(14) = cf; G(15) = rrc; G(16) = rrs;
	if (((F(G, e0) * F(G, ff)) < 0) || ((Lrn - Lam) > pi_2))
	{
		e0 = uu2 - pi; ff = uu2 + pi;
		if (F(G, e0) * F(G, ff) < 0)
		{
			e0 = 0; ff = pi2;
			if (F(G, e0) * F(G, ff) < 0)return;
		}
	}
	if (BiSec(F, G, e0, uu1, eps, eps, T0))goto met1;
	if (BiSec(F, G, uu1, ff, eps, eps, Tf))goto met1;
	goto mek;
met1:
	if (BiSec(F, G, e0, uu2, eps, eps, T0))goto met2;
	if (BiSec(F, G, uu2, ff, eps, eps, Tf))goto met2;
	goto mek;
met2:
	de1 = pi / 12;
	for (j = 1; j <= 18; ++j)
	{
		e1 = e0 + j * de1;
		if (BiSec(F, G, e0, e1, eps, eps, T0))continue;
		if (BiSec(F, G, e1, ff, eps, eps, Tf))continue;
		goto mek;
	}
	u1 = F(G, uu1); u2 = F(G, uu2);
	if (fabs(u1) < fabs(u2))u1 = uu1; else u1 = uu2;
	if (fabs(F(G, u1 + hh)) > fabs(F(G, u1 - hh)))hh = -hh;
	u2 = u1; e0 = 0; e1 = pi2;
	for (j = 1; j <= 255; ++j)
	{
		u1 = u2; u2 = u1 + hh; uu1 = F(G, u1); uu2 = F(G, u2);
		if ((uu1 * uu2) < 0)
		{
			if (BiSec(F, G, e0, u2, eps, eps, T0))return;
			if (BiSec(F, G, u2, e1, eps, eps, Tf))return;
			goto mek;
		}
		if ((fabs(uu2) > fabs(uu1)) && ((uu1 * uu2) > 0))return;
	}
	return;
mek:
	T0 = xos(1) + (T0 - e * sin(T0)) * per;
	Tf = xos(1) + (Tf - e * sin(Tf)) * per;
}

//G94. Определение времени входа и выхода из зоны радиовидимости наземного пункта
void TimeRadio(double F(Vect& G, double e), double eps, double GamMin,
	double ud, double t0, double tf, Vect& xos, Vect& q,
	int& count, Vect& tauP, Vect& T0, Vect& Tf)
{
	int i = 0, N = T0.m; static bool start = 1;
	static double P; double t, tg0, tg, Lrn, T0_ = 0, Tf_ = 0;
	Vect UD0(N), UDf(N); if (start) { P = Pers(xos(3)); start = 0; }
	count = 0; tg0 = t0 - int(t0 / 86400) * 86400; UD0 = ud; UDf = ud;
	Do(t, t0, tf, P)
	{
		i++; if (i > N)return;
		tg = t - int(t / 86400) * 86400; if (tg < tg0)ud++; tg0 = tg;
		Lrn = xos(6) - StimeVer(ud, xos(1)); //Геогр. долгота ВУ
		TimeRadio(F, xos, q, eps, GamMin, Lrn, T0_, Tf_); xos(1) += P;
		if (Tf_ != T0_)
		{
			count++; tauP(count) = xos(1); T0(count) = T0_; Tf(count) = Tf_;
			tg = T0_ - int(T0_ / 86400) * 86400;
			if (tg < tg0)UD0(count) = ud + 1; else UD0(count) = ud;
			tg0 = Tf_ - int(Tf_ / 86400) * 86400;
			if (tg0 < tg)UDf(count) = UD0(count) + 1; else UDf(count) = UD0(count);
			if (T0(i) <= t0)T0(i) = t0; if (Tf(i) > tf)Tf(i) = tf;
		}
	}
}

//G95. Вычисление истинного звёздного времени.			
//Входные параметры: ud - юлианская дата; t - время по Гринвичу в с.
//Выходные параметры: S - истинное звездное время в радианах
//						  для текущей эпохи udt	
double StimeVer(double ud, double t)
{
	double S, eps, eps0, hpsi, heps, udt = ud + t / 86400; 
	Nutn(udt, hpsi, heps);
	eps0 = Eclipte(udt); //Средний наклон эклиптики к экватору
	eps = eps0 + heps;   //Истинный наклон эклиптики к экватору
	S = Stime(udt) + hpsi * cos(eps);
	S -= floor(S / pi2) * pi2; //нормализация
	return S;
}

//G96. Определение сидерического периода обращения по ВС в J2000
double Pers(Vect& x)
{
	double r = !x(1, 3), v = x(4, 6) * x(4, 6), a = r / (2 * mu_ - r * v);
	return 2 * pi * mu_ * sqrt(a * a * a);
}

//G97. Определение сидерического периода обращения по большой полуоси
double Pers(double a) 
{ 
	return 2 * pi * sqrt(a * a * a / mu_); 
}

//G98. Вычисление долго и короткопериодических членов нутации истинного полюса
//в долготе и наклоне.                          
//Входные параметры: ud - юлианская дата.
//Выходные параметры:	                                                 
//			hpsi - нутация в долготе, рад;
//			heps - нутация в наклоне, рад
void Nutn(const double ud, double& hpsi, double& heps)
{
	int kl[106] = { 0,0,-2,2,-2,1,0,2,0,0,0,0,0,2,0,0,0,0,0,-2,0,2,
				 0,1,2,0,0,0,-1,0,0,1,0,1,1,-1,0,1,-1,-1,1,0,2,1,
				 2,0,-1,-1,1,-1,1,0,0,1,1,2,0,0,1,0,1,2,0,1,0,1,
				 1,1,-1,-2,3,0,1,-1,2,1,3,0,-1,1,-2,-1,2,1,1,-2,-1,1,
				 2,2,1,0,3,1,0,-1,0,0,0,1,0,1,1,2,0,0 };
	int kl1[106] = { 0,0,0,0,0,-1,-2,0,0,1,1,-1,0,0,0,2,1,2,-1,0,-1,0,
				  1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,1,1,-1,0,0,0,0,0,0,0,-1,0,1,0,0,1,
				  0,-1,-1,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,
				  0,0,0,0,0,0,1,-1,0,0,1,0,-1,1,0,0,0,1 };
	int kf[106] = { 0,0,2,-2,2,0,2,-2,2,0,2,2,2,0,2,0,0,2,0,0,2,0,
				  2,0,0,-2,-2,0,0,2,2,0,2,2,0,2,0,0,0,2,2,2,0,2,
				  2,2,2,0,0,2,0,2,2,2,0,2,0,2,2,0,0,2,0,-2,0,0,
				  2,2,2,0,2,2,2,2,0,0,0,2,0,0,2,2,0,2,2,2,4,0,
				  2,2,0,4,2,2,2,0,-2,2,0,-2,2,0,-2,0,2,0 };
	int kd[106] = { 0,0,0,0,0,-1,-2,0,-2,0,-2,-2,-2,-2,-2,0,0,-2,0,2,-2,-2,
				  -2,-1,-2,2,2,0,1,-2,0,0,0,0,-2,0,2,0,0,2,0,2,0,-2,
				  0,0,0,2,-2,2,-2,0,0,2,2,-2,2,2,-2,-2,0,0,-2,0,1,0,
				  0,0,2,0,0,2,0,-2,0,0,0,1,0,-4,2,4,-4,-2,2,4,0,-2,
				  -2,2,2,-2,-2,-2,0,2,0,-1,2,-2,0,-2,2,2,4,1 };
	int kom[106] = { 1,2,1,0,2,0,1,1,2,0,2,2,1,0,0,0,1,2,1,1,1,1,
				  1,0,0,1,0,2,1,0,2,0,1,2,0,2,0,1,1,2,1,2,0,2,
				  2,0,1,1,1,1,0,2,2,2,0,2,1,1,1,1,0,1,0,0,0,0,
				  0,2,2,1,2,2,2,1,1,2,0,2,2,0,2,2,0,2,1,2,2,0,
				  1,2,1,2,2,0,1,1,1,2,0,0,1,1,0,0,2,0 };
	double ahci[106] =
	{ -171996,2062,46,11,-3,-3,-2,1,-13187,1426,-517,217,129,48,-22,17,
	 -15,-16,-12,-6,-5,4,4,-4,1,1,-1,1,1,-1,-2274,712,-386,-301,-158,
	 123,63,63,-58,-59,-51,-38,29,29,-31,26,21,16,-13,-10,-7,7,-7,
	 -8,6,6,-6,-7,6,-5,5,-5,-4,4,-4,-3,3,-3,-3,-2,-3,-3,2,-2,
	 2,-2,2,2,1,-1,1,-2,-1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,1,
	 -1,-1,-1,-1,-1,-1,-1,1,-1,1 };
	double bhci[106] =
	{ -174.2,0.2,0,0,0,0,0,0,-1.6,-3.4,1.2,-0.5,0.1,0,0,-0.1,0,0.1,0,0,0,
	 0,0,0,0,0,0,0,0,0,-0.2,0.1,-0.4,0,0,0,0,0.1,-0.1,0.9,0,0,0,0,
	 0,0,0,0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	double aepc[106] =
	{ 92025,-895,-24,0,1,0,1,0,5736,54,224,-95,-70,1,0,0,9,7,6,3,3,-2,
	 -2,0,0,0,0,0,0,0,977,-7,200,129,-1,-53,-2,-33,32,26,27,16,-1,-12,
	 13,-1,-10,-8,7,5,0,-3,3,3,0,-3,3,3,-3,3,0,3,0,0,0,0,
	 0,1,1,1,1,1,-1,1,-1,1,0,-1,-1,0,-1,1,0,-1,1,1,0,0,
	 -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	double bepc[106] =
	{ 8.9,0.5,0,0,0,0,0,0,-3.1,-0.1,-0.6,0.3,0,0,0,0,0,0,0,0,0,0,
	 0,0,0,0,0,0,0,0,-0.5,0,0,-0.1,0,0,0,0,0,0,0,0,0,0,
	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	 0,0,0,0,0,0,0,0,0,0,0,0,0.0,0,0,0,0,0. };
	double t = (ud - 2451545) / 36525;
	//t- время прошедшее от эпохи J 2000.0 в юлианских столетиях
	double t2 = t * t;
	double t3 = t2 * t;
	//fls-средняя аномалия Солнца
	double fls = 6.24003594 + 628.30195602 * t - 2.7974e-6 * t2 - 5.82e-8 * t3;
	//f- средний аргумент широты Луны
	double f = 1.62790193 + 8433.46615831 * t - 6.42717e-5 * t2 - 5.33e-8 * t3;
	//ds- разность средних долгот Луны и Солнца
	double ds = 5.19846951 + 7771.37714617 * t - 3.34085e-5 * t2 + 9.21e-8 * t3;
	//om- средняя долгота восх. узла орбиты Луны на эклиптике
	double om = 2.182438624 - 33.757045936 * t + 3.61429e-5 * t2 - 3.88e-8 * t3;
	//fm-средняя аномалия Луны
	double fm = 2.355548393 + 8328.69142288 * t + 1.517952e-4 * t2 + 3.103e-7 * t3;
	double apg; hpsi = 0; heps = 0;
	for (int i = 0; i <= 105; ++i)
	{
		apg = kl[i] * fm + kl1[i] * fls + kf[i] * f + kd[i] * ds + kom[i] * om;
		hpsi = hpsi + (ahci[i] + bhci[i] * t) * sin(apg);
		heps = heps + (aepc[i] + bepc[i] * t) * cos(apg);
	}
	hpsi = 4.848136811e-10 * hpsi; heps = 4.848136811e-10 * heps;
}

//G99. Вычисление среднего наклона экватора к эклиптике с учётом прецессии     
//Входные параметры: ud - юлианская дата.
//Выходные параметры: средний наклон эклиптики, рад
double Eclipte(const double ud)
{
	double t = (ud - 2451545) / 36525, t2 = t * t;
	return 0.4090928042 - 0.2269655e-3 * t - 0.29e-8 * t2 + 0.88e-8 * t2 * t;
}

//G100. След матрицы
double Trace(Matr& A)
{
	if (A.m != A.n) { Ecrir(cout, "\n\n Trace: argument - not square matrix"); END }
	int i; double t = 0; DO(i, 1, A.m)t += A(i, i); return t;
}

//G101. Модификация BiSec для использования в функции TimeRadio.
//Определение нуля функции Func методом дихотомии
//(деления отрезка пополам) на отрезке [t0,tf], 
//если на [t0,tf] - нечетное количество нулей.
//epsf, epst - критерии завершения вычислений
//по значению функции и значению аргумента соответственно;
//kod - код завершения вычислений:
//	    если kod=1, то на [t0,tf] нуля функции не существует
//	    или они в четном количестве
bool BiSec(double Func(Vect& G, double t), Vect& G, double t0,
		   double tf, double epsf, double epst, double& t)
{
	bool kod; double y, z;
	kod = 0; y = Func(G, t0);
	if (fabs(y) <= epsf) { t = t0; return kod; };
	z = Func(G, tf); if (fabs(z) <= epsf) { t = tf; return kod; };
	if (y * z > 0) { kod = 1; t = 0; return kod; }
iter:
	t = 0.5 * (t0 + tf); y = Func(G, t);
	if (y * z > 0)tf = t; else t0 = t;
	if (fabs(t0 - tf) >= epst)goto iter;
	return kod;
}

//G102. Определение оскулирующих элементов орбиты
//по прямоугольным инерциальным координатам.
//xos(1)=tau, xos(2)=e, xos(3)=a, xos(4)=omega, xos(5)=i, xos(6)=Omega
Vect& ElemOrb(double mu, double t, Vect& x)
{
	double mc, mc2, r, q, si, u, p, teta, T;
	Vect c(3), rx(3), rv(3); VSTACK(6, xos)
		 rx = x(1, 3); rv = x(4, 6); r = !rx; c = VecProd(rx, rv);
	mc2 = c * c; mc = sqrt(mc2); q = c(3) / mc; xos(5) = acos(q); si = sin(xos(5));
	if (si == 0) { xos(6) = 0; u = atan2(x(2) * q, x(1)); }
	else
	{
		xos(6) = atan2(c(1), -c(2));
		u = atan2(x(3) / si, x(1) * cos(xos(6)) + x(2) * sin(xos(6)));
	}
	xos(3) = mu * r / (2 * mu - r * rv * rv); p = mc2 / mu;
	if (p / xos(3) >= 1)xos(2) = 0; else xos(2) = sqrt(1 - p / xos(3));
	if (xos(2) == 0) { xos(4) = 0; teta = u; }
	else { teta = atan2(rx * rv * mc / mu, p - r); xos(4) = u - teta; }
	xos(1) = t - TimePrg(mu, 0, xos(2), xos(3), teta);
	//Нормализация времени прохождения перицентра:
	if (xos(1) < 0) { T = pi2 * sqrt(xos(3) * xos(3) * xos(3) / mu); xos(1) += T; }
	if (xos(4) < 0)xos(4) = xos(4) + pi2; if (xos(6) < 0)xos(6) = xos(6) + pi2;
	STACKh(6, x)
	*(x->h) = xos;
	return *(x->h);
}

//G103. Определение времени прохождения углового расстояния,
//равного истинной аномалии.
//Внимание: поскольку 0<=teta<=2*pi, то значению teta<0
//ставится в соответствие угол teta+2*pi		  
double TimePrg(double mu, double tau, double e, double a, double teta)
{
	double T, q, TP; T = sqrt(a * a * a / mu);
	if (e > 1)
	{
		q = sqrt(e * e - 1) * sin(teta) / (1 + e * cos(teta));
		TP = tau + (e * q - log(q + sqrt(q * q + 1))) * T;
		return TP;
	}
	q = sqrt(1 - e * e) * sin(teta) / (1 + e * cos(teta));
	TP = tau + (Asin(q, e + cos(teta)) - e * q) * T;
	return TP;
}

//G104. Переход от оскулирующих элементов орбиты
//к прямоугольным инерциальным координатам.
//x(1)=tau, x(2)=e, x(3)=a, x(4)=omega, x(5)=i, x(6)=OMEGA
Vect& Kepler(double mu, double eps, double t, Vect& x)
{
	double som, com, si, ci, teta, ct, r, u, su, cu;
	Vect y(6);
	som = sin(x(6)); com = cos(x(6)); si = sin(x(5)); ci = cos(x(5));
	teta = Anom(mu, eps, x(1), x(2), x(3), t); ct = cos(teta);
	r = fabs(x(3) * (1 - x(2) * x(2))) / (1 + x(2) * ct);
	u = x(4) + teta; su = sin(u); cu = cos(u);
	y(1) = r * (cu * com - som * su * ci); y(2) = r * (cu * som + su * com * ci);
	y(3) = r * su * si; r = sqrt(mu / r * (2 - fabs(r / x(3)) * sign(1 - x(2))));
	u = u + 1.570796326794896 - atan(x(2) * sin(teta) / (1 + x(2) * ct));
	su = sin(u); cu = cos(u); y(4) = r * (cu * com - som * su * ci);
	y(5) = r * (cu * som + su * com * ci); y(6) = r * su * si;
	STACKh(6, x)
	*(x->h) = y;
	return *(x->h);
}

//G105. Переход от оскулирующих элементов орбиты
//к прямоугольным инерциальным координатам.
//x(1)=tau, x(2)=e, x(3)=a, x(4)=omega, x(5)=i, x(6)=OMEGA
Vect& Kepler(double mu, double eps, double t, Vect& x, double& u)
{
	double som, com, si, ci, teta, ct, r, uv, su, cu;
	Vect y(6);
	som = sin(x(6)); com = cos(x(6)); si = sin(x(5)); ci = cos(x(5));
	teta = Anom(mu, eps, x(1), x(2), x(3), t); ct = cos(teta);
	r = fabs(x(3) * (1 - x(2) * x(2))) / (1 + x(2) * ct);
	u = x(4) + teta; su = sin(u); cu = cos(u);
	y(1) = r * (cu * com - som * su * ci); y(2) = r * (cu * som + su * com * ci);
	y(3) = r * su * si; r = sqrt(mu / r * (2 - fabs(r / x(3)) * sign(1 - x(2))));
	uv = u + 1.570796326794896 - atan(x(2) * sin(teta) / (1 + x(2) * ct));
	su = sin(uv); cu = cos(uv); y(4) = r * (cu * com - som * su * ci);
	y(5) = r * (cu * som + su * com * ci); y(6) = r * su * si;
	STACKh(6, x)
	*(x->h) = y;
	return *(x->h);
}

//G106. Вычисление матрицы перехода от J2000
//к J-ТЭ с учётом прецессии и нутации
//Входные параметры:
//ud - юлианская дата; t - время по Гринвичу в с.
//Выходные параметры:
//A(3,3) - матрица перехода от J2000 к J-ТЭ
Matr& Reduct(double ud, double t)
{
	double eps, eps0, hpsi, heps, udt = ud + t / 86400;
	Matr A(3, 3);
	Nutn(udt, hpsi, heps);
	eps0 = Eclipte(udt);			//Cредний наклон эклиптики к экватору
	eps = eps0 + heps;				//Истинный наклон эклиптики к экватору
	A = MatNut(udt) * MatPrec(udt); //C нутацией
	STACKH(3, 3, E)
	*(E->H) = A;
	return *(E->H);
}

//G107. Вычисление матрицы перехода от J2000
//к J-ТЭ с учётом прецессии и нутации.
//Входные параметры:	                                  
//ud - юлианская дата; t - время по Гринвичу в с.         
//Выходные параметры:	                                  
//A(6,6) - матрица перехода от J2000 к J-ТЭ			
Matr& ReductG(double ud, double t)
{
	Matr A(3, 3), B(6, 6);
	A = Reduct(ud, t);
	B = A; B(A, 4, 4);
	STACKH(6, 6, E)
	*(E->H) = B;
	return *(E->H);
}

//G108. Вычисление истинного звёздного времени и матрицы перехода от J2000
//к J-ТЭ с учётом прецессии и нутации.
//Входные параметры:	                                             
//ud - юлианская дата; t - время по Гринвичу в с.                 
//Выходные параметры:	                                             
//S - истинное звездное время в радианах для текущей эпохи udt;	
//A - матрица перехода от J2000 к J-ТЭ							
Matr& Reduct(double ud, double t, double& S)
{
	double eps, eps0, hpsi, heps, udt = ud + t / 86400;
	Matr A(3, 3);
	Nutn(udt, hpsi, heps);
	eps0 = Eclipte(udt); //Средний наклон эклиптики к экватору
	eps = eps0 + heps;     //Истинный наклон эклиптики к экватору
	S = Stime(udt) + hpsi * cos(eps); //Истинное звездное время
	S -= floor(S / pi2) * pi2; //Нормализация
	A = MatNut(udt) * MatPrec(udt); //СMSTACK(3, 3, A)  нутацией
	STACKH(3, 3, E)
	*(E->H) = A;
	return *(E->H);
}

/*//G109. Вычисление матрицы поворота относительно одной из заданных осей.
//Входные параметры                                          	  
//k  - номер оси, вокруг которой производится вращение:           
//     1 - ось X,                                           
//     2 - ось Y,                                           
//     3 - ось Z                                            
//u  - угол поворота в рад                                    
//Выходные параметры:								              
//A(3,3) - матрица направляющих косинусов,описывающих вращение
Matr& MatRot(int k, double u)
{
	double s = sin(u), c = cos(u);
	Ecrir(cout,"\n			MatRot MatRot MatRot MatRot MatRot");
	Matr A(3, 3);
	switch (k)
	{
		case 1:
		{ A(1, 1) = 1; A(2, 3) = s; A(2, 2) = c; A(3, 2) = -s; A(3, 3) = c; return A; }
		case 2:
		{ A(1, 1) = c; A(1, 3) = -s; A(2, 2) = 1; A(3, 1) = s; A(3, 3) = c; return A; }
		case 3:
		{ A(1, 1) = c; A(1, 2) = s; A(2, 2) = c; A(2, 1) = -s; A(3, 3) = 1; return A; }
	}
	STACKH(3, 3, E)
	*(E->H) = A;		
	return *(E->H);
}*/

//G109. Вычисление матрицы поворота относительно одной из заданных осей.
//Входные параметры                                          	  
//k  - номер оси, вокруг которой производится вращение:           
//     1 - ось X,                                           
//     2 - ось Y,                                           
//     3 - ось Z                                            
//u  - угол поворота в рад                                    
//Выходные параметры:								              
//A(3,3) - матрица направляющих косинусов,описывающих вращение
Matr& MatRot(int k, double u)
{
	double s = sin(u), c = cos(u);
	Matr A(3, 3);
	switch (k)
	{
		case 1:
		A(1, 1) = 1; A(2, 3) = s; A(2, 2) = c; A(3, 2) = -s; A(3, 3) = c; break;
		case 2:
		A(1, 1) = c; A(1, 3) = -s; A(2, 2) = 1; A(3, 1) = s; A(3, 3) = c; break;
		case 3:
		A(1, 1) = c; A(1, 2) = s; A(2, 2) = c; A(2, 1) = -s; A(3, 3) = 1;
	}
	STACKH(3, 3, E)
	*(E->H) = A;
	return *(E->H);
}

//G110. Вычисление матрицы перехода от J2000 к ГСК
//Входные параметры:	                                         
//ud - юлианская дата; t - гринвичское время в с				
//Выходные параметры:	                                         
//B(3,3) - матрица перехода от J2000 к ГСК					
Matr& Mat2000Gr(double ud, double t)
{
	double S; 
	Matr A(3, 3), B(3, 3);
	A = Reduct(ud, t, S);
	B = MatRot(3, S) * A;
	STACKH(3, 3, E)
	*(E->H) = B;
	return *(E->H);
}

//G111. Вычисление матрицы перехода от проекций скорости в J2000
//к проекциям скорости в ГСК.
//Входные параметры:	                                                   
//ud - юлианская дата; t - гринвичское время в с						
//Выходные параметры:	                                                   
//C(3,6) - матрица перехода от проекций скорости в J2000				
//         к проекциям скорости в ГСК								
Matr& Mat2000GrV(double ud, double t)
{
	double s, ss, cs; 
	Matr A(3, 3), B(3, 3), Q(3, 3), C(3, 6);
	B = Reduct(ud, t, s); ss = sin(s); cs = cos(s);
	Q(1, 1) = -ss; Q(1, 2) = cs; Q(2, 1) = -cs; Q(2, 2) = -ss; Q = omE_ * Q;
	A = MatRot(3, s) * B; B = Q * B; //B - производная от A (B=dA/dt)
	C = B; C(A, 1, 4); 
	STACKH(3, 6, E)
	*(E->H) = C;
	return *(E->H);
}

//G112. Вычисление матрицы перехода от проекций скорости в ГСК 
//к проекциям скорости в J2000. 						 
//Входные параметры:	                                 
//ud - юлианская дата; t - гринвичское время в с							
//Выходные параметры:	                                                    
//C(3,6) - матрица перехода от проекций скорости в J2000					
//         к проекциям скорости в ГСК									
Matr& MatGr2000V(double ud, double t)
{
	double s, ss, cs; 
	Matr A(3, 3), B(3, 3), Q(3, 3), C(3, 6);
	B = Reduct(ud, t, s); ss = sin(s); cs = cos(s);
	Q(1, 1) = -ss; Q(1, 2) = cs; Q(2, 1) = -cs; Q(2, 2) = -ss; Q = omE_ * Q;
	A = MatRot(3, s) * B; B = Q * B; //B - производная от A
	C = !B; C(!A, 1, 4); 
	STACKH(3, 6, E)
	*(E->H) = C;
	return *(E->H);
}

//G113. Вычисление матрицы перехода от вектора состояния в J2000
//к вектору состояния в ГСК.
//Входные параметры:	                                                
//ud - юлианская дата; t - гринвичское время в с						
//Выходные параметры:	                                                
//C(6,6) - матрица перехода от вектора состояния в J2000				
//         к вектору состояния в гринвичской СК							
Matr& Mat2000GrG(double ud, double t)
{
	double s, ss, cs;
	Matr A(3, 3), B(3, 3), Q(3, 3), C(6, 6);
	B = Reduct(ud, t, s); ss = sin(s); cs = cos(s);
	Q(1, 1) = -ss; Q(1, 2) = cs; Q(2, 1) = -cs; Q(2, 2) = -ss; Q = omE_ * Q;
	A = MatRot(3, s) * B; B = Q * B; //B - производная от A
	C = A; C(B, 4, 1); C(A, 4, 4); 
	STACKH(6, 6, E)
	*(E->H) = C;
	return *(E->H);
}

//G114. Вычисление матрицы перехода от вектора состояния в ГСК
//к вектору состояния в J2000.
//Входные параметры:	                                                     
//ud - юлианская дата; t - гринвичское время в с							 
//Выходные параметры:	                                                     
//C(6,6) - матрица перехода от вектора состояния в ГСК 					 
//         к вектору состояния в J2000									 
Matr MatGr2000G(double ud, double t)
{
	double s, ss, cs; 
	Matr A(3, 3), AT(3, 3), B(3, 3), Q(3, 3), C(6, 6);
	B = Reduct(ud, t, s); ss = sin(s); cs = cos(s);
	Q(1, 1) = -ss; Q(1, 2) = cs; Q(2, 1) = -cs; Q(2, 2) = -ss; Q = omE_ * Q;
	A = MatRot(3, s) * B; AT = !A; B = Q * B; //B - производная от A
	C = AT; C(!B, 4, 1); C(AT, 4, 4); 
	STACKH(6, 6, E)
	*(E->H) = C;
	return *(E->H);
}

//G115. Вычисление матрицы прецессии.                    
//Входные параметры: ud - юлианская дата.
//Выходные параметры: P - матрица прецессии
Matr& MatPrec(double ud)
{
	double dz, z, tet, sd, cd, sz, cz, st, ct, t, t2, t3;
	Matr P(3, 3);
	t = (ud - 2451545.0) / 36525;
	// t- время прошедшее от эпохи J 2000 в юлианских столетиях
	t2 = t * t; t3 = t2 * t;
	dz = 0.0111808609 * t + 0.146356e-5 * t2 + 0.872e-7 * t3;
	z = 0.0111808609 * t + 0.53072e-5 * t2 + 0.883e-7 * t3;
	tet = 0.0097171735 * t - 0.20685e-5 * t2 + 0.2028e-6 * t3;
	sd = sin(dz); cd = cos(dz); sz = sin(z); cz = cos(z);
	st = sin(tet); ct = cos(tet);
	P(1, 1) = cd * cz * ct - sd * sz; P(1, 2) = -sd * cz * ct - cd * sz; P(1, 3) = -cz * st;
	P(2, 1) = cd * sz * ct + sd * cz; P(2, 2) = -sd * sz * ct + cd * cz; P(2, 3) = -sz * st;
	P(3, 1) = cd * st;          P(3, 2) = -sd * st;          P(3, 3) = ct;
	STACKH(3, 3, E)
	*(E->H) = P;
	return *(E->H);
}

//G116. Вычисление матрицы нутации.                         
//Входные параметры:
//ud - юлианская дата.	                                                 
//np - признак выбора формул для расчета нутации в долготе и наклоне:     
//     np=0 - значения нутации вычисляются по полным формулам;              
//     np=1 - значения нутации определяются с точностью до членов           
//     порядка 0".1 (или 1e-6).                                           
//nm - признак выбора формул для расчета матрицы нутации:                 
//     nm=0 - матрица нутации вычисляется по точным формулам;               
//     nm=1 - матрица нутации вычисляется по упрощенным формулам            
//     с точностью до малых членов порядка 1e-8                           
//Выходные параметры: N - матрица нутации
Matr MatNut(double ud)
{
	double epc, hpci, hepc, epc0, spci, cpci, s, c, sepc, cepc;
	Matr N(3,3);
	epc0 = Eclipte(ud); // средний наклон эклиптики к экватору
	Nutn(ud, hpci, hepc);
	epc = epc0 + hepc;    // истинный наклон эклиптики к экватору
	spci = sin(hpci); cpci = cos(hpci); s = sin(epc0);
	c = cos(epc0); sepc = sin(epc); cepc = cos(epc);
	N(1, 1) = cpci; N(1, 2) = -spci * c; N(1, 3) = -spci * s;
	N(2, 1) = spci * cepc; N(2, 2) = cpci * cepc * c + sepc * s;
	N(2, 3) = cpci * cepc * s - sepc * c; N(3, 1) = spci * sepc;
	N(3, 2) = cpci * sepc * c - cepc * s; N(3, 3) = cpci * sepc * s + cepc * c;
	STACKH(3, 3, E)
	*(E->H) = N;
	return *(E->H);
}

//G117. Вычисление гринвичского среднего звёздного времени              
//на заданное время суток.                  
//Входные параметры: udt - юлианская дата плюс
//						   текущая часть суток по гринвичу.
//Выходные параметры: звёздное время на заданное время суток, рад
double Stime(double udt)
{
	double jd = floor(udt + 0.5) - 0.5, //Юлианская дата на начало суток
	d = udt - 2451545, tau = d / 36525, S;
	S = 1.7533685592 + 0.0172027918051 * d + 6.2831853072 * (udt - jd) +
		6.7707139e-6 * tau * tau - 4.50876e-10 * tau * tau * tau;
	return S - floor(S / pi2) * pi2;
}

//G118. Вычисление среднего звёздного времени гринвичского меридиана на начало суток
//Входные параметры: ud - юлианская дата.
//Выходные параметры: звёздное временя, рад
double Snul(double ud)
{
	double jd = floor(ud + 0.5) - 0.5, d = ud - 2451545, tau = d / 36525, S;
	S = 1.7533685592 + 0.0172027918051 * d + 6.7707139e-6 * tau * tau -
		4.50876e-10 * tau * tau * tau;
	return S - floor(S / pi2) * pi2;
}

//G119. Перевод календарной даты в юлианскую 
//Входные параметры: d - день; m - месяц; y - год.
//Выходные параметры: юлианская дата на начало суток по Гринвичу
double Ulian(int d, int m, int y)
{
	double ud = floor(30.57 * m) + d;
	ud += floor((5 - floor(m / 3)) / 5) * (2 - floor((4 - fmod(y, 4)) / 4) +
		  floor((1e2 - fmod(y, 1e2)) / 1e2) - floor((400 - fmod(y, 4e2)) / 4e2));
	return ud += 1721027.5 + 365 * y + floor(y / 4) - floor(y / 1e2) + floor(y / 4e2);
}

//G120. Перевод календарной даты в юлианскую 
//Входные параметры: d - день; m - месяц; y - год.
//Выходные параметры: юлианская дата на полдень по Гринвичу
double Uday(int d, int m, int y)
{
	int c, ya; double ud;
	if (m > 2)m -= 3; else { m += 9; y--; }
	c = y / 100; ya = y - 100 * c;
	ud = 146097 * c / 4 + 1461 * ya / 4 + (153 * m + 2) / 5 + d + 1721119;
	return ud;
}

//G121. Полином степени n-1 (порядка n) от x 				              
double Polynom(Vect& c, double x)
{
	int k, n = c.m; double p = 1, P = 0;
	DO(k, 0, n - 1) { if (k > 0)p *= x; P += c(k + 1) * p; }
	return P;
}

//G122.
void PolinomM(double x, Matr& a, Vect& f)
{
	int i, j; 
	DO(i, 1, 3) 
	{ 
		f(i) = a(10, i); 
		dO(j, 9, 1, -1)f(i) = f(i) * x + a(j, i); 
	}
}

//G123.
void PolinomS(double x, Matr& a, Vect& f)
{
	int i, j; 
	DO(i, 1, 3) 
	{ 
		f(i) = a(11, i); 
		dO(j, 10, 1, -1)f(i) = f(i) * x + a(j, i); 
	}
}

//G124. Переход от декартовых координат к сферическим координатам
//(размерности векторов d и возвращаемого s могут равняться 3 или 6,
//т.е. d=[x,y,z] или d=[x,y,z,Vx,Vy,Vz];  
//s=[r,fi,lam] или s=[r,fi,lam,Vr,Vfi,Vlam] 
Vect& DecSpher(Vect& d)
{
	int m = d->m;
	double R, sl, cl;
	Vect s(m);
	STACKh(m, e)
	R = !d(1, 2); s(1) = !d(1, 3);
	s(2) = asin(d(3) / s(1)); sl = d(2) / R; cl = d(1) / R;
	s(3) = Asin(sl, cl);
	*(e->h) = s;
	if (m == 3)return *(e->h);
	s(4) = d(1, 3) * d(4, 6) / s(1);
	s(5) = (d(6) * s(1) - d(3) * s(4)) / s(1) / R;
	s(6) = (d(5) - d(2) * d(1, 2) * d(4, 5) / R / R) / d(1);
	if (sl >= 0 && cl < 0 || sl < 0 && cl <= 0)s(6) = -s(6);
	*(e->h) = s;
	return *(e->h);
}

//G125. Переход от сферических координат к декартовым координатам
//(размерности векторов s и возвращаемого d могут равняться 3 или 6,
//т.е. d=[x,y,z] или d=[x,y,z,Vx,Vy,Vz];  
//s=[r,fi,lam] или s=[r,fi,lam,Vr,Vfi,Vlam] 
Vect& SpherDec(Vect& s)
{
	int m = s->m;
	double sf = sin(s(2)), cf = cos(s(2)), sl = sin(s(3)), cl = cos(s(3));
	Vect d(m);
	STACKh(m, e)
	d(1) = s(1) * cf * cl; d(2) = s(1) * cf * sl; d(3) = s(1) * sf;
	*(e->h) = d;
	if (m == 3)return *(e->h);
	d(4) = s(4) * cf * cl - s(1) * (sf * cl * s(5) + cf * sl * s(6));
	d(5) = s(4) * cf * sl + s(1) * (cf * cl * s(6) - sf * sl * s(5));
	d(6) = s(4) * sf + s(1) * cf * s(5);
	*(e->h) = d;
	return *(e->h);
}

//G126. Вычисление значений полинома Лежандра порядка n
//для любого действительного аргумента x.
//Рекурсивная функция
double Legendre(int n, double x)
{
	if (n == 0)return 1; if (n == 1)return x;
	return double(2 * n - 1) / n * x * Legendre(n - 1, x) - double(n - 1) / n * Legendre(n - 2, x);
}

//G127. Функция Legendr вычисляет значения присоединенных функций
//Лежандра первого рода Pij(x) (i=2,3,...,n; j=0,1,...,i)
//для вещественных аргументов
void Legendr(int n, double x, Matr& P)
{
	int i, j, k; double y = sqrt(1 - x * x);
	P(1, 1) = 1; if (n == 0)return; P(2, 1) = x; P(2, 2) = y;
	DO(i, 2, n)
	{
		j = 2 * (i - 1) + 1;
		P(i + 1, 1) = (j * x * P(i, 1) - (i - 1) * P(i - 1, 1)) / i;
		P(i + 1, i + 1) = j * y * P(i, i);
	}
	DO(i, 2, n)
	{
		k = i + 1;
		DO(j, 2, k - 1)P(k, j) = x * P(i, j) + (i + j - 2) * y * P(i, j - 1);
	}
}

//G128. Матрица перехода от декартовой к сферической системе координат
Matr& MatDecSpher(double fi, double lam)
{
	double sf = sin(fi), sl = sin(lam), cf = cos(fi), cl = cos(lam);
	Matr M(3, 3);
	M(1, 1) = cf * cl;  M(1, 2) = cf * sl;  M(1, 3) = sf;
	M(2, 1) = -sf * cl; M(2, 2) = -sf * sl; M(2, 3) = cf;
	M(3, 1) = -sl;      M(3, 2) = cl;
	STACKH(3, 3, E)
	*(E->H) = M;
	return *(E->H);
}

//G129. Матрица перехода от декартовой к сферической системе координат
Matr MatDecSpher(Vect& x)
{
	double r = !x, R = !x(1, 2), sf = x(3) / r, cf = R / r,
		   sl = x(2) / R, cl = x(1) / R;
	Matr M(3, 3);
	M(1, 1) = cf * cl;  M(1, 2) = cf * sl;  M(1, 3) = sf;
	M(2, 1) = -sf * cl; M(2, 2) = -sf * sl; M(2, 3) = cf;
	M(3, 1) = -sl;      M(3, 2) = cl;
	STACKH(3, 3, E)
	*(E->H) = M;
	return *(E->H);
}

//G130. Статическая атмосфера (модель ABCD-59 [Нариманов, Тихонравов])
double StAtmos(double H)
{
	double ro, a0 = -17.295, h0 = 104000, h = H, kappa = 0.01594;
	//Приведенное в модели значение a0=-16.72 заменено на a0=-17.295
	ro = exp(a0 - kappa * sqrt(h - h0)); 
	return ro;
}

//G131. Вектор кинетического момента нормированный
Vect& VectKinMom(double i, double Om)
{
	double ci = cos(i), si = sin(i), cOm = cos(Om), sOm = sin(Om);
	Vect r1(3), r2(3), c(3);
	r1(1) = rE_ * cOm; r1(2) = rE_ * sOm; r1(3) = 0;
	r2(1) = -rE_ * sOm * ci; r2(2) = rE_ * cOm * ci; r2(3) = rE_ * si;
	c = VecProd(r1, r2); c = c / !c;
	STACKh(3, e)
	*(e->h) = c;
	return *(e->h);
}

//G132. Определение координат Луны в J2000
void Moon(double udt, Vect& r)
{
	bool startMoon = 1;
	int i, j, NInt, NStr;
	double DTJ0, D_DTJ, DTJ, NI, e;
	Matr KFT(3, 10), KF(10, 3);
	DTJ0 = 2449600; D_DTJ = udt - DTJ0;
	NInt = floor(D_DTJ / 20); NStr = NInt * 3;
	NI = NInt; DTJ = D_DTJ - NI * 20;
	if (startMoon)
	{
		DO(i, 1, NStr)Read(InMoon, e, e, e, e, e, e, e, e, e, e);
		ReadMatr(InMoon, KFT); InMoon.seekg(ios::beg);
	}
	KF = !KFT;	if (muM_ > 1e9)KF = KF * 1000;
	PolinomM(DTJ, KF, r);
}

//G133. Определение координат Солнца в J2000
void Sun(double udt, Vect& r)
{
	bool startSun = 1;
	int i, j, NInt, NStr;
	double DTJ0, D_DTJ, DTJ, NI, e;
	Matr KFT(3, 11), KF(11, 3);
	DTJ0 = 2449600; D_DTJ = udt - DTJ0;
	NInt = floor(D_DTJ / 100); NStr = NInt * 3; NI = NInt;
	DTJ = D_DTJ - NI * 100;
	if (startSun)
	{
		DO(i, 1, NStr)Read(InSun, e, e, e, e, e, e, e, e, e, e, e);
		ReadMatr(InSun, KFT); InSun.seekg(ios::beg);
	}
	KF = !KFT;
	if (muS_ > 1e15)KF = KF * 1000; PolinomS(DTJ, KF, r);
}

//G134. Матрица перехода из с.к. J2000 в орбитальную систему координат,
//оси которой направлены по трансверсали (1), радиусу-вектору (2),
//и бинормали к плоск. орбиты (3, дополняет до правой).
//Входные параметры:
//xos - вектор оскулирующих координат КА xos=[e,om,a,i,Om,u].
//Возвращаемые параметры:
//A(3,3) - матрица перехода из J2000 в орбитальную с.к.
Matr& Mat2000Orb(Vect& xos)
{
	double su = sin(xos(6)), cu = cos(xos(6)), si = sin(xos(4)), ci = cos(xos(4)),
		   sOm = sin(xos(5)), cOm = cos(xos(5));
	Matr A(3, 3);
	A(1, 1) = -cOm * su - sOm * cu * ci; A(1, 2) = -sOm * su + cOm * cu * ci; A(1, 3) = cu * si;
	A(2, 1) = cOm * cu - sOm * su * ci;  A(2, 2) = sOm * cu + cOm * su * ci;  A(2, 3) = su * si;
	A(3, 1) = -sOm * si;                 A(3, 2) = cOm * si;                  A(3, 3) = -ci;
	STACKH(3, 3, E)
	*(E->H) = A;
	return *(E->H);
}

//G135. Матрица перехода из с.к. J2000 в орбитальную с.к.,
//оси которой направлены по трансверсали (1), радиусу-вектору (2),
//и бинормали к плоск. орбиты (3, дополняет до правой).
//Входные параметры:
//xos - вектор оскулирующих координат КА xos=[e,om,a,i,Om,u].
//Возвращаемые параметры:
//B(6,6) - матрица перехода из с.к. J2000 в орбитальную с.к.
Matr& Mat2000OrbG(Vect& xos)
{
	Matr A(3, 3), B(6, 6);
	A = Mat2000Orb(xos);
	B = A; B(A, 4, 4); 
	STACKH(6, 6, E)
	*(E->H) = B;
	return *(E->H);
}

//G136. Матрица перехода из с.к. J2000 в орбитальную с.к.,
//оси которой направлены по трансверсали (1), радиусу-вектору (2),
//и бинормали к плоск. орбиты (3, дополняет до правой).
//Входные параметры:
//x(6) - вектор сосотояния КА в J2000.
//Выходные параметры:
//A(3,3) - матрица перехода из с.к. J2000 в орбитальную с.к.
Matr& MatAbsOrb(Vect& x)
{
	Vect xos(6); 
	Matr A(3, 3); 
	xos = AbsOsc(x); //xos=[e,om,a,i,Om,u]
	double su = sin(xos(6)), cu = cos(xos(6)), si = sin(xos(4)), ci = cos(xos(4)),
		   sOm = sin(xos(5)), cOm = cos(xos(5));
	A(1, 1) = -cOm * su - sOm * cu * ci; A(1, 2) = -sOm * su + cOm * cu * ci; A(1, 3) = cu * si;
	A(2, 1) = cOm * cu - sOm * su * ci;  A(2, 2) = sOm * cu + cOm * su * ci;  A(2, 3) = su * si;
	A(3, 1) = -sOm * si;           A(3, 2) = cOm * si;            A(3, 3) = -ci;
	STACKH(3, 3, E)
	*(E->H) = A;
	return *(E->H);
}

//G137. Матрица перехода из с.к. J2000 в орбитальную с.к.,
//оси которой направлены по трансверсали (1), радиусу-вектору (2),
//и бинормали к плоск. орбиты (3, дополняет до правой).
//Входные параметры:
//x(6) - вектор состояния КА (фазовый вектор) в J2000.
//Возвращаемые параметры:
//B(6,6) - матрица перехода из с.к. J2000 в орбитальную с.к.
Matr& MatAbsOrbG(Vect& x)
{
	double dt = 1, r = !x(1, 3); Vect x1(6), x2(6), w(3);
	Matr A(3, 3), V(3, 3), B(6, 6);
	A = MatAbsOrb(x); r = r * r * r; w = -mu_ * x(1, 3) / r;
	x1 = x(1, 3) - x(4, 6) * dt; x2 = x(1, 3) + x(4, 6) * dt;
	x1(x(4, 6) - w * dt, 4); x2(x(4, 6) + w * dt, 4);
	V = 0.5 * !(MatAbsOrb(x2) - MatAbsOrb(x1)) / dt;
	B = A; B(A, 4, 4); B(V, 4, 1);
	STACKH(6, 6, E)
	*(E->H) = B;
	return *(E->H);
}

//G138. Матрица перехода из орбитальной с.к. в с.к. J2000.
//Оси орбитальной с.к. направлены по трансверсали (1),
//радиусу-вектору (2) и бинормали к плоск. орбиты (3).
//Входные параметры:
//x(6) - вектор состояния (фазовый вектор) КА в J2000.
//Возвращаемые параметры:
//B(6,6) - матрица перехода из орбитальной с.к. в с.к. J2000 
Matr& MatOrbAbsG(Vect& x)
{
	double dt = 1, r = !x(1, 3); Vect x1(6), x2(6), w(3);
	Matr A(3, 3), V(3, 3), B(6,6);
	A = MatAbsOrb(x); r = r * r * r; w = -mu_ * x(1, 3) / r;
	x1 = x(1, 3) - x(4, 6) * dt; x2 = x(1, 3) + x(4, 6) * dt;
	x1(x(4, 6) - w * dt, 4); x2(x(4, 6) + w * dt, 4);
	V = 0.5 * !(MatAbsOrb(x2) - MatAbsOrb(x1)) / dt;
	A = !A; B = A; B(A, 4, 4); B(-A * V * A, 4, 1);
	STACKH(6, 6, E)
	*(E->H) = B;
	return *(E->H);
}

//G139. Интегрирование системы дифференциальных уравнений
//в прямом и обратном направлениях.
//Входные параметры:                                                      
//F		- вектор-функция правых частей системы уравнений f(t,x)
//eps	- относительная точность интегрирования;	                      
//dt0	- начальный шаг интегрирования;	                              
//t0	- начальное значение независимой переменной;	                              
//tf	- конечное значение независимой переменной;	                              
//x		- начальный вектор интегрируемых параметров.	                      
//Выходные параметры:
//dt0	- конечный шаг интегрирования;	                                  
//x		- конечный вектор интегрируемых параметров             	          
void BulSto(Vect& F(double t, Vect& x), double eps,
			double& dt0, double t0, double tf, Vect& x)
{
	if (t0 == tf)return;
	bool bb, bh, fin = 0, ex = 0, rvr = 0, konv;
	int i, j, k, r, sr, l, m, jj, kk, n = x.m;
	double a, b, b1, c, g, hh, fc, t = t0, ta, u, v;
	Vect d(7), dy(n), dz(n), ya(n), yl(n), ym(n);
	//s	- массив максимальных шагов интегрирования	                      
	static Vect s(n); Matr dt(n, 7), yg(8, n), yh(8, n);
	if ((tf - t0) * dt0 < 0)dt0 = -dt0;
	if (tf < t0)rvr = 1; if (tf <= t + dt0)ex = 1;
	if (rvr) { if (ex)ex = 0; else ex = 1; }
	if (ex) { dt0 = tf - t; hh = dt0; fin = 1; }
Restep:
	dz = F(t, x); bh = 0; ya = x;
Anf:
	a = dt0 + t; fc = 1.5; bb = 0; m = 1; r = 2; sr = 3; jj = 0; dt = 0;
	DO(j, 0, 9)
	{
		if (bb)
		{
			d(2) = 1.77777777777778;
			d(4) = 7.111111111111111;
			d(6) = 28.44444444444444;
		}
		else { d(2) = 2.25; d(4) = 9; d(6) = 36; }
		konv = j > 2; if (j > 6) { l = 6; d(7) = 64; fc = 0.6 * fc; }
		else { l = j; d(l + 1) = m * m; }
		m *= 2; g = dt0 / m; b = g * 2;
		if (bh && j < 8) { ym = yh(j + 1); yl = yg(j + 1); }
		else
		{
			kk = (m - 2) / 2; m--; yl = ya; ym = ya + g * dz;
			DO(k, 1, m)
			{
				dy = F(t + k * g, ym);
				DO(i, 1, n)
				{
					u = yl(i) + b * dy(i); yl(i) = ym(i);
					ym(i) = u; u = fabs(u); if (u > s(i))s(i) = u;
				}
				if (k == kk && k != 2) { jj++; yh(jj, ym); yg(jj, yl); }
			}
		}
		dy = F(a, ym);
		DO(i, 1, n)
		{
			v = dt(i, 1); ta = 0.5 * (ym(i) + yl(i) + g * dy(i)); c = ta; dt(i, 1) = c;
			DO(k, 1, l)
			{
				b1 = d(k + 1) * v; b = b1 - c; u = v;
				if (b != 0) { b = (c - v) / b; u = c * b; c = b1 * b; }
				v = dt(i, k + 1);
				dt(i, k + 1) = u;
				ta += u;
			}
			if (fabs(x(i) - ta) > eps * s(i))konv = 0;
			x(i) = ta;
		}
		if (konv)goto End; d(3) = 4; d(5) = 16; bb = bb ? 0 : 1;
		m = r; r = sr; sr = m * 2;
	}
	bh = bh ? 0 : 1; fin = 0; dt0 *= 0.5; goto Anf;
End:
	if (fin) { dt0 = hh; return; }
	dt0 *= fc; t = a; if (tf < (t + 1.1 * dt0))ex = 1; else ex = 0;
	if (rvr) { if (ex)ex = 0; else ex = 1; }
	if (ex) { hh = dt0; fin = 1; dt0 = tf - t; }
	goto Restep;
}

//G140. Получение равномерно распределённого на отрезке [0, 1] псевдослучайного числа 
double Rrandom(double& x)
{
	//100000000001<=x<=137438953471
	double r; //if (x < 0)x = StartRand();
	r = 137438953472e0; x = x * 3125; x = x - (int)(x / r) * r;
	x = x * 625 + 29044268343e0; x = x - (int)(x / r) * r;
	return x / r;
}

//G141. Получение нормально распределённого псевдослучайного числа 
//с нулевым мат. ожиданием и единичным с.к.о.
//100000000001<=x<=137438953471
double Nrandom(double& x)
{
	double s; s = 0; for (int i = 1; i <= 5; i++)s += Rrandom(x);
	s = 7.745966692414834e-1 * (2 * s - 5); return s + 0.01 * (pow(s, 3) - 3 * s);
}

//G142. Получение реализации случайного, нормально распределенного
//n-мерного вектора с математическим ожиданием x
//и матрицей ковариаций P
//100000000001<=rand<=137438953471
Vect& RandVec(double& rand, Vect& x0, Matr& P)
{
	int i, n = P.n; Matr A(n, n); VSTACK(n, x)
	A = P->RootMat(); for (i = 1; i <= n; i++)x(i) = Nrandom(rand);
	x = x0 + A * x; return x;
}

//G143. Инициализация датчиков ПСВ с использованием текущего времени
double StartRand()
{
	fstream out("StartRand.txt", ios::out | ios::in | ios::trunc);
	out.precision(12); time_t t; time(&t); Str S(24);
	char hour[2], min[2], sec[2]; int h, m, s;
	//S = ctime(&t); не работает в VS2019
	hour[0] = S(12); hour[1] = S(13); h = atoi(hour);
	min[0] = S(15);  min[1] = S(16);  m = atoi(min);
	sec[0] = S(18);  sec[1] = S(19);  s = atoi(sec); s = h * 3600 + m * 60 + s;
	double Start = floor(s / 86400e0 * 37438953470e0 + 100000000001e0);
	if (fmod(Start, 2) < 0.999999)Start++;
	Ecrir(cout, "\n StartRand=", Start);
	Ecrir(out, "\n StartRand=", Start);
	out.close();
	return Start;
}

//G144. Определение матрицы частных производных A(n,m)
//от функции Func размерности m по аргументу x размерности n
Matr& DifFun(Vect& Func(Vect& x), int m, Vect& x, Vect& dx)
{
	int i, n = x.m; double r;
	Vect f1(m), f2(m); MSTACK(n, m, A)
		DO(i, 1, n)
	{
		r = x(i); x(i) = r - dx(i); f1 = Func(x); x(i) = r + dx(i);
		f2 = Func(x); x(i) = r; A(i, 0.5 * (f2 - f1) / dx(i));
	}
	return A;
}

//G145. Определение транспонированной матрицы частных производных A(m,n)
//функции размерности m по аргументу x размерности n
Matr& DifFunT(Vect& Func(Vect& x), int m, Vect& x, Vect& dx)
{
	int i, n = x.m; double r; Vect f1(m), f2(m); MSTACK(m, n, A)
		DO(i, 1, n)
	{
		r = x(i); x(i) = r - dx(i); f1 = Func(x); x(i) = r + dx(i);
		f2 = Func(x); x(i) = r; A(0.5 * (f2 - f1) / dx(i), i);
	}
	return A;
}

//G146. Определение градиента фунции f(x) конечно-разностным методом
//в точке x0 с приращением dx
Vect& GradF(double f(Vect& x), Vect& x0, Vect& dx)
{
	int i, m = x0->m; double xi, f1, f2; VSTACK(m, grad)
		DO(i, 1, m)
	{
		xi = x0(i); x0(i) = xi + dx(i); f2 = f(x0);
		x0(i) = xi - dx(i); f1 = f(x0); x0(i) = xi;
		grad(i) = 0.5 * (f2 - f1) / dx(i);
	}
	return grad;
}

//G147. Определение производной от скалярной функции
//по скалярному аргументу в точке x
double DifFun(double Func(double x), double x, double dx)
{
	return 0.5 * (Func(x + dx) - Func(x - dx)) / dx;
}

//G148. Определение производной n-го порядка
//от скалярной функции скалярного аргумента
//(рекурсивная функция)
double dnf_dxn(double f(double x), int n, double x, double dx)
{
	double dfdx;
	if (n == 1)
	{
		dfdx = (f(x + dx) - f(x - dx)) / 2 / dx;
		return dfdx;
	}
	dfdx = (dnf_dxn(f, n - 1, x + dx, dx) - dnf_dxn(f, n - 1, x - dx, dx)) / 2 / dx;
	return dfdx;
}

//G149. Белый шум на [0,dt].
//rand - инициатор датчика ПСЧ,
//100000000001<=rand<=137438953471;
//D - интенсивность белого шума
double BlancBruit(double& rand, double D, double dt)
{
	double B = sqrt(D / dt), x = Nrandom(rand);
	return B * x;
}

//G150. Векторный белый шум на [t1,t2].
//Отрезок [t1,t2] разбивается на N частей.
//Для каждой из частей генерерируется векторный БШ.
//Из N векторов выбирается вектор, относящийся к отрезку,
//в пределах которого располагается t.
//Перед обращением к BlancBruit необходимо,
//чтобы marker=1. 
//rand - инициатор датчика ПСЧ,
//100000000001<=rand<=137438953471;
//D - матрица интенсивностей белого шума
//	  (диагональная)
Vect& BlancBruit(bool& marker, double& rand, int N,
				 Matr& D, double t, double t1, double t2)
{
	int i, j, n = D->m;
	double dt = (t2 - t1) / N;
	VSTACK(n, x)
	static Matr X(N, n);
	if (marker)
	{
		DO(i, 1, N)DO(j, 1, n)
			X(i, j) = BlancBruit(rand, D(j, j), dt);
		marker = 0;
	}
	DO(i, 1, N)if (t1 + (i - 1) * dt <= t && t <= t1 + i * dt)x = X(i);
	return x;
}

//G151. Нулевая матрица 
Matr& Zero(int m, int n) { MSTACK(m, n, A) return A; }

//G152. Нулевой вектор 
Vect& Zero(int m) { VSTACK(m, a) return a; }

//G153. Критерий Сильвестра
//Возвращает 1, если матрица положительно определенная;
//возвращает -1, если матрица отрицательно определенная;
//возвращает 0, если квадратичная форма знакопеременная
//или полупределённая
int Silvestr(Matr& P)
{
	int k, m = P->m; double d0 = P(1, 1), d; Matr P0(m, m);
	P0 = P;
	if (d0 > 0)
	{
		DO(k, 2, m)
		{
			P0->m = k; P0->n = k; //Это вполне корректно (Денис Звягинцев)
			d = det(P0);
			if (d0 * d > 0)d0 = d;
			else return 0;
		}
		return 1;
	}
	if (d0 < 0)
	{
		DO(k, 2, m)
		{
			P0->m = k; P0->n = k;
			d = det(P0);
			if (d0 * d < 0)d0 = d;
			else return 0;
		}
		return -1;
	}
	return 0;
}

//G154. Восстановление положительной определённости
//матрицы P с точностью gamma>0.
//Возвращает 1, если на входе матрица P>0;
//возвращает 2, если на входе матрица P неположительно
//				определённая, а на выходе P>0;
//возвращает 0, если на входе и выходе матрица P 
//неположительно определённая,
//т.е. восстановление положительной
//определённости не произошло.
int Repositive(double eps, double gamma, Matr& P)
{
	if (Silvestr(P) == 1) return 1;
	int i, n = P->m; Vect Lam(n); Matr S(n, n);
	Jacobi(P, S, Lam, eps);
	DO(i, 1, n)if (Lam(i) <= 0)Lam(i) = gamma;
	P = S * Dmat(Lam) * !S;
	if (Silvestr(P) == 1)return 2;
	else return 0;
}

//G155. Собственные значения и собственные векторы симметричной матрицы
void Jacobi(Matr& B, Matr& S, Vect& lambda, double eps)
{
	double norm1 = 0, norm2 = 0, mu = 0, in = 0, thr = 0, omega = 0,
		v1 = 0, v2 = 0, v3 = 0, sint = 0, cost = 0;
	int i = 0, j = 0, k = 0, ind = 0, p = 0, q = 0, n = lambda.m; Matr A(n, n);
	A = B; for (i = 1; i <= n; i++) S(i, i) = 1;
	for (i = 2; i <= n; i++) for (j = 1; j <= i - 1; j++) in = in + 2 * pow(A(i, j), 2);
	if (in == 0) goto lab1;
	norm1 = sqrt(in); thr = norm1; norm2 = eps / n * norm1; ind = 0;
lab2:
	thr = thr / n;
lab3:
	for (q = 2; q <= n; q++)
		for (p = 1; p <= q - 1; p++)
		{
			if (fabs(A(p, q)) >= thr)
			{
				ind = 1; v1 = A(p, p); v2 = A(p, q); v3 = A(q, q); mu = 0.5 * (v1 - v3);
				if (mu == 0) omega = -1;
				else omega = -sign(mu) * v2 / sqrt(pow(v2, 2) + pow(mu, 2));
				sint = omega / sqrt(2 * (1 + sqrt(1 - pow(omega, 2))));
				cost = sqrt(1 - pow(sint, 2));
				for (i = 1; i <= n; i++)
				{
					if (i != p && i != q)
					{
						in = A(i, p); mu = A(i, q); A(i, q) = in * sint + mu * cost;
						A(q, i) = A(i, q); A(i, p) = in * cost - mu * sint; A(p, i) = A(i, p);
					}
					in = S(i, p); mu = S(i, q);
					S(i, q) = in * sint + mu * cost;
					S(i, p) = in * cost - mu * sint;
				}
				mu = pow(sint, 2); omega = pow(cost, 2); in = sint * cost;
				A(p, p) = v1 * omega + v3 * mu - 2 * v2 * in;
				A(q, q) = v1 * mu + v3 * omega + 2 * v2 * in;
				A(q, p) = (v1 - v3) * in + v2 * (omega - mu);
				A(p, q) = A(q, p);
			}
		}
	if (ind == 1) { ind = 0; goto lab3; }
	if (thr > norm2) goto lab2;
lab1:
	for (i = 1; i <= n; i++)
	{
		lambda(i) = 0;
		for (j = 1; j <= n; j++) for (k = 1; k <= n; k++)
			lambda(i) = lambda(i) + S(j, i) * B(j, k) * S(k, i);
	}
}

//G156. Разложение матрицы P на матрицу собсвенных векторов S и вектор собственных значений
void Decomp(fstream& out, double eps, Matr& P, Matr& S, Vect& Lam)
{
	Jacobi(P, S, Lam, fabs(eps));
	if (eps < 0)
	{
		Ecrir(out,"\n\n            Decomposition of matrix: Silvestr(P)=", Silvestr(P));
		MatrOut(out,"P",P,"S",S,10);
		VecOut(out,"Lam",Lam,10);
	}
}

//G157. Вектор с компонентами от 1 до n
Vect& VN(int n) { int i; VSTACK(n, v) DO(i, 1, n)v(i) = i; return v; }

//G158.
Vect& VS(double s) { VSTACK(1, v) v(1) = s; return v; }

//G159.
Vect& VS(int m, double s)
{
	int i; VSTACK(m, v) DO(i, 1, m)v(i) = s;
	return v;
}

//G160. Перезапись строк матрицы в обратном порядке
void Invert(Matr& A)
{
	int k, m = A->m, n = A->n; Matr B(m, n);
	DO(k, 1, m) B(k, A(m - k + 1)); A = B;
}

Vect& Adding(Vect& x, double y) //G161. Прямое сложение вектора и скаляра
{
	int mx = x->m, n = mx + 1;
	VSTACK(n, xy)
	xy = x; xy(n) = y;
	return xy;
}

Vect& Adding(double y, Vect& x) //G162. Прямое сложение скаляра и вектора
{
	int mx = x->m, n = mx + 1;
	VSTACK(n, xy)
	xy(1) = y; xy(2, x);
	return xy;
}

Vect Adding(Vect& x, Vect& y, Vect& z) //G163. Прямое сложение трёх векторов
{
	int mx = x.m, my = y.m, mz = z.m, n = mx + my + mz; VSTACK(n, xyz)
		xyz(x, 1); xyz(y, mx + 1); xyz(z, mx + my + 1);
	return xyz;
}

//G164. Разложение вектора x на фрагменты x1, x2
void Scatter(Vect& x, Vect& x1, Vect& x2)
{
	int i, m1 = x1->m, m2 = x2->m;
	if (m1 + m2 > x->m)
	{
		Ecrir(cout, "\n\n Scatter(Vect& x,Vect& x1,Vect& x2):\n");
		Ecrir(cout, " Amount of dimensionality x1,x2 exceeds dimensionality x");
		END
	}
	DO(i, 1, m1)x1(i) = x(i);
	DO(i, m1 + 1, m1 + m2)x2(i - m1) = x(i);
}

//G165. Разложение вектора x на фрагменты x1, x2, x3
void Scatter(Vect& x, Vect& x1, Vect& x2, Vect& x3)
{
	int i, m1 = x1->m, m2 = x2->m, m3 = x3->m;
	if (m1 + m2 + m3 > x->m)
	{
		Ecrir(cout, "\n\n Scatter(Vect& x,Vect& x1,Vect& x2,Vect& x3):\n");
		Ecrir(cout, " Amount of dimensionality x1,x2,x3 exceeds dimensionality x");
		END
	}
	DO(i, 1, m1)x1(i) = x(i);
	DO(i, m1 + 1, m1 + m2)x2(i - m1) = x(i);
	DO(i, m1 + m2 + 1, m1 + m2 + m3)x3(i - m1 - m2) = x(i);
}

//G166. Факториал (рекурсивная)
long int Fact(int x)
{
	if (x == 0 || x == 1) return 1;
	return x * Fact(x - 1);
}

//G167. Поэлементное извлечение корня вектора 
Vect& Sqrt(Vect& x)
{
	int k, n = x->m;
	VSTACK(n,s);
	DO(k, 1, n)
	{
		if (x(k) < 0)
		{
			Ecrir(cout,"\n Vect& Sqrt(Vect& x): x(k)<0! x(k)=",x(k)," k=",k,"\n");
			END
		}
		s(k) = sqrt(x(k));
	}
	return s;
}

//G168. Поэлементное извлечение корня из матрицы
Matr& Sqrt(Matr& A)
{
	int i, j, m=A->m, n=A->n;
	MSTACK(m,n,S);
	DO(i, 1, m)DO(j, 1, n)
	{
		if (A(i,j) < 0)
		{
			Ecrir(cout, "\n Matr& Sqrt(Matr& A): A(i,j)<0! A(i,j)=", A(i, j),
				  " i=", i, " j=", j, "\n");
			END
		}
		S(i, j) = A(i, j);
	}
	return S;
}

//G169. Удаляет из матрицы I-ю строку
Matr Cut(int I, Matr& P)
{
	int i, m1 = P->m - 1, n = P->n; MSTACK(m1, n, Q)
		DO(i, 1, m1)if (i < I)Q(i, P(i)); else Q(i, P(i + 1));
	return Q;
}

//G170. Удаляет из матрицы J-й столбец
Matr Cut(Matr& P, int J)
{
	int i, m = P->m, n1 = P->n - 1; MSTACK(m, n1, Q)
		DO(i, 1, n1)if (i < J)Q(P[i], i); else Q(P[i + 1], i);
	return Q;
}

//G171. Удаляет из матрицы I-ю строку и J-й столбец
Matr Cut(Matr& P, int I, int J)
{
	int i, j, k, l, m1 = P->m - 1, n1 = P->n - 1; MSTACK(m1, n1, Q)
		DO(i, 1, m1)
	{
		if (i < I)k = i; else k = i + 1;
		DO(j, 1, n1)
		{
			if (j < J)l = j; else l = j + 1;
			Q(i, j) = P(k, l);
		}
	}
	return Q;
}

//G172. Удаляет из вектора p k-й элемент
Vect Cut(Vect& p, int k)
{
	int i, m1 = p->m - 1; VSTACK(m1, q)
		DO(i, 1, m1)q(i) = i < k ? p(i) : p(i + 1);
	return q;
}

//Переводит целое десятичное x в строку символов,
//соответствующую значению x в двоичной системе
void DecBit(int x, char bit[])
{
	int i;
	for (i = 0; i < 16; i++)
	{
		if (x & 0x8000)bit[i] = '1';
		else bit[i] = '0';
		x <<= 1;
	}
	bit[i] = 0;
}

//Переводит целое десятичное x в строку символов,
//соответствующую значению x в двоичной системе
char* DecBit(int x)
{
	int i = 0, size = 0, xx = abs(x);
	size = NumberBase(xx, 2);
	char* bin = new char[size + 2];
	int x_2 = xx;
	for (i = size; i >= 1; i--)
	{
		if (x_2 % 2 == 0) bin[i] = '0'; else bin[i] = '1';
		x_2 = x_2 / 2;
	}
	bin[size + 1] = 0;
	if (x < 0) bin[0] = '-'; else bin[0] = ' ';
	return bin;
}

//Переводит целое десятичное x в строку символов,
//соответствующую значению x в системе с
//заданным основанием base<=32
char* DecBase(int x, int base)
{
	int i = 0, size = 0, xx = abs(x);
	size = NumberBase(xx, base);
	char numeral[32] =
	{
		'0','1','2','3','4','5','6','7','8','9',
		'A','B','C','D','E','F','G','H','I','J',
		'K','L','M','N','O','P','Q','R','S','T',
		'U','V'
	};
	char* b = new char[size + 2];
	int x_b = xx;
	for (i = size; i >= 1; i--)
	{
		for (int j = 0; j < 32; j++)
			if (x_b % base == j) { b[i] = numeral[j]; x_b = x_b / base; break; }
	}
	b[size + 1] = 0;
	if (x < 0)b[0] = '-'; else b[0] = '+';
	return b;
}

//Переводит целое десятичное x в строку символов,
//соответствующую значению x в системе с
//заданным основанием base<=32
char* DecBase(double x, int base)
{
	int i = 0, size = 0, xx = abs(x);
	size = NumberBase(xx, base);
	char numeral[32] =
	{
		'0','1','2','3','4','5','6','7','8','9',
		'A','B','C','D','E','F','G','H','I','J',
		'K','L','M','N','O','P','Q','R','S','T',
		'U','V'
	};
	char* b = new char[size + 2];
	int x_b = xx;
	for (i = size; i >= 1; i--)
	{
		for (int j = 0; j < 32; j++)
			if (x_b % base == j) { b[i] = numeral[j]; x_b = x_b / base; break; }
	}
	b[size + 1] = 0;
	if (x < 0)b[0] = '-'; else b[0] = '+';
	return b;
}

//Возвращает количество разрядов целого
//десятичного x в системе с заданным основанием base
int NumberBase(int x, int base)
{
	if (x == 0)return 1;
	else
	{
		int xx = 1, size = 0;
		while (xx <= x) { xx = base * xx; size++; }
		return size;
	}
}

