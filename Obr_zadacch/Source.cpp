#define _USE_MATH_DEFINES
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

/// <summary>
/// Вывод матрицы в файл
/// </summary>
/// <param name="A"> - Матрица для вывода</param>
/// <param name="fileName"> - Имя файла, куда выводить (Относительное или абсолютное)</param>
void writeMatrix(vector<vector<double>> A, string fileName)
{
	ofstream stream;
	stream.open(fileName);
	int n = A.size();
	for (int q = 0; q < n; q++)
	{
		for (int p = 0; p < n; p++)
		{
			stream << A[q][p] << " ";
		}
		stream << endl;
	}
	stream.close();
}

/// <summary>
/// Вывод матрицы в консоль
/// </summary>
/// <param name="A"> - Матрица для вывода</param>
void writeMatrix(vector<vector<double>> A)
{
	int n = A.size();
	for (int q = 0; q < n; q++)
	{
		for (int p = 0; p < n; p++)
		{
			cout << A[q][p] << " ";
		}
		cout << endl;
	}
}

/// <summary>
/// Вывод вектора в файл
/// </summary>
/// <param name="v"> - Вектор для вывода</param>
/// <param name="fileName"> - Имя файла, куда выводить (Относительное или абсолютное)</param>
void writeVector(vector<double> v, string fileName)
{
	ofstream stream;
	stream.open(fileName);
	int n = v.size();
	for (int i = 0; i < n; i++)
		stream << v[i] << endl;
	stream.close();
}

/// <summary>
/// Вывод вектора в консоль
/// </summary>
/// <param name="v"> - Вектор для вывода</param>
void writeVector(vector<double> v)
{
	int n = v.size();
	for (int i = 0; i < n; i++)
		cout << v[i] << endl;
}

class solvn
{
public:
	solvn() {};
	~solvn() {};

	/// <summary>
	/// Инициализация модели (Чтение фходных данных)
	/// </summary>
	void model_init();

	/// <summary>
	/// Решение обратной задачи
	/// </summary>
	void model_solve();

private:
	int n_params = 3;
	int max_iter = 100;
	double eps_func = 1E-8;
	double eps_dif = 1E-14;
	//double alpha = 1E-20;
	double alpha = 1E-8;
	double sigma = 0;

	int flag_F; // Иницилизированна ли модель или нет

	double Fq;
	//double Acoeff_prev, Bcoeff_prev, Ccoeff_prev;

	struct CoeffStruct
	{
		double Acoeff, Bcoeff, Ccoeff;

		friend CoeffStruct operator+(CoeffStruct C1, CoeffStruct C2)
		{
			CoeffStruct tmpCoeff;
			tmpCoeff.Acoeff = C1.Acoeff + C2.Acoeff;
			tmpCoeff.Bcoeff = C1.Bcoeff + C2.Bcoeff;
			tmpCoeff.Ccoeff = C1.Ccoeff + C2.Ccoeff;

			return tmpCoeff;
		}

		friend CoeffStruct operator*(CoeffStruct C1, double aC)
		{
			CoeffStruct tmpCoeff;
			tmpCoeff.Acoeff = C1.Acoeff * aC;
			tmpCoeff.Bcoeff = C1.Bcoeff * aC;
			tmpCoeff.Ccoeff = C1.Ccoeff * aC;

			return tmpCoeff;
		}

		CoeffStruct operator=(CoeffStruct C1)
		{
			this->Acoeff = C1.Acoeff;
			this->Bcoeff = C1.Bcoeff;
			this->Ccoeff = C1.Ccoeff;

			return *this;
		}
	};

	vector<double>  Fval,
		Fval_true, // Истинные значения функции на сетке
		Tval; // Узлы сетки

	CoeffStruct TrueCoeffs, // Истинные значения искомых коэффициентов
		ResCoeffs,
		DiffCoeffs,
		PrevCoeffs;

	/// <summary>
	/// Зануление матрицы СЛАУ и правой части
	/// </summary>
	/// <param name="A"> - Матрица СЛАУ</param>
	/// <param name="b"> - Вектор правой части</param>
	void Clear_SLAU(vector<vector<double>>& A, vector<double>& b);

	/// <summary>
	/// Решение СЛАУ методом Гаусса (Недореализован)
	/// </summary>
	/// <param name="A"> - Матрица СЛАУ</param>
	/// <param name="b"> - Вектор правой части</param>
	/// <param name="x"> - Вектор для записи решения</param>
	/// <param name="N"> - Размерность СЛАУ</param>
	/// <returns></returns>
	int Gauss(vector<vector<double>>& A, vector<double>& b, vector<double>& x, int N);

	/// <summary>
	/// Расчёт нормы вектора
	/// </summary>
	/// <param name="point1"> - Вектор, норму которого необходимо расчитать</param>
	/// <returns> - Норма вектора</returns>
	double Get_norm(vector<double>& point1);
	//------------------------------------

	/// <summary>
	/// Расчёт приращения коэффициента для численного расчёта производных целевой функции от коэффициентов
	/// </summary>
	void Calc_dCoeff();
	double Get_func(vector<double> MFval, vector<double> MFval_true);

	/// <summary>
	/// Аналитический расчёт производных целевой функции по искомым коэффициентам
	/// </summary>
	/// <param name="t"> - Точка, значение производной функции по параметру в которой нужно рассчитать</param>
	/// <param name="MyCoeffs"> - Коэффициенты, по которым дифференциируется функция</param>
	/// <param name="numCoeff"> - Номер параметра по которому ищеться производная</param>
	/// <returns></returns>
	double dF_dt(double t, CoeffStruct MyCoeffs, int numCoeff);

	/// <summary>
	/// Численный расчёт производных целевой функции по искомым коэффициентам
	/// </summary>
	/// <param name="t"> - Точка, значение производной функции по параметру в которой нужно рассчитать</param>
	/// <param name="MyCoeffs"> - Коэффициенты, по которым дифференциируется функция</param>
	/// <param name="numCoeff"> - Номер параметра по которому ищеться производная</param>
	/// <returns></returns>
	double F_dCoeff(double t, CoeffStruct MyCoeffs, int numCoeff);

	/// <summary>
	/// Целевая функция
	/// </summary>
	/// <param name="t"> - Значение аргумента, по которому необходимо вычислить значение функции</param>
	/// <param name="MyCoeffs"> - Коэффициенты для вычисления целевой функции</param>
	/// <returns></returns>
	double Formul(double t, CoeffStruct MyCoeffs);
	/// <summary>
	/// Вроде как тоже целевая функция
	/// </summary>
	/// <param name="t"> - Значение аргумента в котором необходимо узнать значение функции</param>
	/// <param name="MyCoeffs"> - Значение коэффициентов целевой функции</param>
	/// <returns> - Значение целевой функции</returns>
	double Receive_F(double t, CoeffStruct MyCoeffs);

	/// <summary>
	/// Вроде как тоже целевая функция с зависимостью от решения
	/// </summary>
	/// <param name="t"> - Значение аргумента в котором необходимо узнать значение функции</param>
	/// <param name="MyCoeffs"> - Значение коэффициентов целевой функции</param>
	/// <returns> - Значение целевой функции</returns>
	double Receive_F(double t, CoeffStruct MyCoeffs, vector<double>& b);
	int Input_data(vector<double>& MFval, vector<double>& MFval_true, vector<double>& MTval, string path);
	int Direct_task(vector<double>& MFval, vector<double>& MTval, CoeffStruct MyCoeffs, string path);
	//int Direct_task(vector<double>& MFval, vector<double>& MTval, CoeffStruct MyCoeff, string path, vector<double> b);
	int Inverse_task(string path);
};

//----------------------------------------
double solvn::Get_func(vector<double> MFval, vector<double> MFval_true)
{
	double func = 0;
	//for (int i = 0; i < MFvec.size(); i++)
	//{
		//func += (MFvec[i].Tval - MFvec[i].Tval_prev);
	//}
	//func *= alpha;

	/*for (int i = 0; i < MFvec.size(); i++)
	{
		MFvec[i].Fval = MFvec[i].Tval* MFvec[i].Tval*Acoeff+ MFvec[i].Tval*Bcoeff+Ccoeff;
	}*/


	for (int i = 0; i < MFval.size(); i++)
	{
		func += (MFval[i] - MFval_true[i]) * (MFval[i] - MFval_true[i]) / (MFval_true[i] * MFval_true[i]);
	}

	func += alpha * (PrevCoeffs.Acoeff - ResCoeffs.Acoeff) * (PrevCoeffs.Acoeff - ResCoeffs.Acoeff);
	func += alpha * (PrevCoeffs.Bcoeff - ResCoeffs.Bcoeff) * (PrevCoeffs.Acoeff - ResCoeffs.Acoeff);
	func += alpha * (PrevCoeffs.Ccoeff - ResCoeffs.Ccoeff) * (PrevCoeffs.Acoeff - ResCoeffs.Acoeff);

	return func;
}
double solvn::dF_dt(double t, CoeffStruct MyCoeffs, int numCoeff)
{
	double dF = 0;

	switch (numCoeff)
	{
	case 0:
		//dF = 2 * MyCoeffs.Acoeff * t;
		dF = t * t;
		//dF = t * t * MyCoeffs.Acoeff;
		//dF = cos(MyCoeffs.Acoeff) * t;
		break;
	case 1:
		//dF = 2 * MyCoeffs.Bcoeff;
		//dF = t;
		//dF = -sin(MyCoeffs.Bcoeff);
		break;
	case 2:
		//dF = 1;
		break;
	default:
		break;
	}

	return dF;
}

void solvn::Calc_dCoeff()
{
	DiffCoeffs = (ResCoeffs * 0.005);

	if (DiffCoeffs.Acoeff < eps_func)
		DiffCoeffs.Acoeff = eps_func;
	if (DiffCoeffs.Bcoeff < eps_func)
		DiffCoeffs.Bcoeff = eps_func;
	if (DiffCoeffs.Ccoeff < eps_func)
		DiffCoeffs.Ccoeff = eps_func;
}

double solvn::Formul(double t, CoeffStruct MyCoeffs)
{
	//return (MyCoeffs.Acoeff * MyCoeffs.Acoeff * t + MyCoeffs.Bcoeff * MyCo-effs.Bcoeff);
	//return (t * t * MyCoeffs.Acoeff + t * MyCoeffs.Bcoeff + MyCoeffs.Ccoeff);
	//return (MyCoeffs.Acoeff * MyCoeffs.Acoeff * t * t);
	//return (log(MyCoeffs.Acoeff) * t * t);
	//return (sin(MyCoeffs.Acoeff)*t+cos(MyCoeffs.Bcoeff));
	//return sin(MyCoeffs.Acoeff)* t + cos(MyCoeffs.Bcoeff);
	return MyCoeffs.Acoeff * MyCoeffs.Acoeff * t + MyCoeffs.Bcoeff * MyCoeffs.Bcoeff + MyCoeffs.Ccoeff;
}

double solvn::F_dCoeff(double t, CoeffStruct MyCoeffs, int numCoeff)
{
	double dF = 0;
	Calc_dCoeff();

	CoeffStruct SummCoeff = (ResCoeffs + DiffCoeffs);
	CoeffStruct t_Coeff = ResCoeffs;

	switch (numCoeff)
	{
	case 0:
		// A
		t_Coeff.Acoeff += DiffCoeffs.Acoeff;
		dF = (Formul(t, t_Coeff) - Formul(t, ResCoeffs)) / DiffCoeffs.Acoeff;
		break;
	case 1:
		// B
		t_Coeff.Bcoeff += DiffCoeffs.Bcoeff;
		dF = (Formul(t, t_Coeff) - Formul(t, ResCoeffs)) / DiffCoeffs.Bcoeff;
		break;
	case 2:
		// С
		t_Coeff.Ccoeff += DiffCoeffs.Ccoeff;
		dF = (Formul(t, t_Coeff) - Formul(t, ResCoeffs)) / DiffCoeffs.Ccoeff;
		break;
	default:
		break;
	}

	return dF;
}


double solvn::Receive_F(double t, CoeffStruct MyCoeffs)
{
	//MFval = MyCoeffs.Acoeff * t * t + MyCoeffs.Bcoeff * t + MyCoeffs.Ccoeff;
	//MFval = (t * t * MyCoeffs.Acoeff + t * MyCoeffs.Bcoeff + MyCoeffs.Ccoeff);
	//MFval = (MyCoeffs.Acoeff * MyCoeffs.Acoeff *t*t);
	//return (log(MyCoeffs.Acoeff) * t * t);
	//MFval = (sin(MyCoeffs.Acoeff)*t+cos(MyCoeffs.Bcoeff));
	//return 0;
	//return sin(MyCoeffs.Acoeff) * t + cos(MyCoeffs.Bcoeff);
	return MyCoeffs.Acoeff * MyCoeffs.Acoeff * t + MyCoeffs.Bcoeff * MyCoeffs.Bcoeff + MyCoeffs.Ccoeff;
}

double solvn::Receive_F(double t, CoeffStruct MyCoeffs, vector<double>& b)
{

	//MFval = (MyCoeffs.Acoeff + b[0]) * t * t + (MyCoeffs.Bcoeff+ b[1]) * t + My-Coeffs.Ccoeff + b[2];
	//MFval = (t * t * (MyCoeffs.Acoeff + b[0]) + t * (MyCoeffs.Bcoeff + b[1]) + My-Coeffs.Ccoeff + b[2]);
	//return (t * t * log(MyCoeffs.Acoeff + b[0]));
	//MFval = ((sin(MyCoeffs.Acoeff + b[0]))*t+ (cos(MyCoeffs.Bcoeff+ b[1])));
	//return 0;
	//return sin(MyCoeffs.Acoeff + b[0]) * t + cos(MyCoeffs.Bcoeff + b[1]);
	return MyCoeffs.Acoeff * MyCoeffs.Acoeff * t + MyCoeffs.Bcoeff * MyCoeffs.Bcoeff + MyCoeffs.Ccoeff;
}

int solvn::Input_data(vector<double>& MFval, vector<double>& MFval_true, vector<double>& MTval, string path)
{
	ifstream inf;
	int tmp_int, coeff, n_c;
	double tmp_dbl;

	// Чтение сетки
	inf.open(path + "tnet.txt", ios_base::in);
	if (inf)
	{
		inf >> tmp_int; // Количество узлов
		//=========================
		MFval.resize(tmp_int);
		MFval_true.resize(tmp_int);
		MTval.resize(tmp_int);
		//=========================

		//Чтение узлов
		for (int i = 0; i < tmp_int; i++)
			inf >> Tval[i];
	}
	inf.close();
	inf.clear();

	// Чтение параметров решателя
	inf.open(path + "params.txt", ios_base::in);
	if (inf)
	{

		inf >> n_params; // Количество параметров
		inf >> max_iter; // Максимальное количество итераций при решении обратной задачи
		inf >> eps_func; // Минимальное значение функционала
		inf >> eps_dif; // Минимальоне значение нормы вектора b (Может быть невязки?)

	}
	inf.close();
	inf.clear();

	// Считываем способ получения истинных значений функции (true - из файла, false - программно, просто подставляем истинные значения параметров в уравнение)
	inf.open(path + "fvalues.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "fvalues.txt";
		return 1;
	}
	inf >> flag_F; // Способ получения истинных значений функции (true - из файла, false - программно, просто подставляем истинные значения параметров в уравнение)
	//============================
	//MFvec.resize(tmp_int);
	//============================

	if (flag_F) {
		for (int i = 0; i < tmp_int; i++)
		{
			inf >> Fval_true[i]; // Истинные значения функции на сетке
		}
	}
	else
	{
		inf.close();
		inf.clear();
		// Считываем истинные коэффициенты из файла
		inf.open(path + "coeffs.txt", ios_base::in);

		inf >> TrueCoeffs.Acoeff; // Истинное значение коэффициента при t^2 
		inf >> TrueCoeffs.Bcoeff; // Истинное значение коэффициента при t
		inf >> TrueCoeffs.Ccoeff; // Истинное значение свободного коэффициента 
	}
	inf.close();
	inf.clear();

	// Считываем начальное приближение коэффициентов
	inf.open(path + "startv.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "start_value.txt";
	}
	else
	{
		inf >> ResCoeffs.Acoeff; // Начальное приближение коэффициента при t^2 
		inf >> ResCoeffs.Bcoeff; // Начальное приближение коэффициента при t
		inf >> ResCoeffs.Ccoeff; // Начальное приближение свободного коэффициента 
	}

	inf.close();
	inf.clear();

	cout << "-----True Parameters-----" << endl;
	cout << "True A = " << TrueCoeffs.Acoeff << endl;
	cout << "True B = " << TrueCoeffs.Bcoeff << endl;
	cout << "True C = " << TrueCoeffs.Ccoeff << endl;
	cout << endl;

	cout << "-----Init mean of Parameters-----" << endl;
	cout << "Init A = " << ResCoeffs.Acoeff << endl;
	cout << "Init B = " << ResCoeffs.Bcoeff << endl;
	cout << "Init C = " << ResCoeffs.Ccoeff << endl;
	cout << endl;

	return 0;
}

int solvn::Direct_task(vector<double>& MFval, vector<double>& MTval, CoeffStruct MyCoeffs, string path)
{
	ofstream opf;
	for (int i = 0; i < MFval.size(); i++)
	{
		MFval[i] = Receive_F(MTval[i], MyCoeffs);
	}
	opf.open(path + "funcsres.txt", ios_base::app);

	for (int i = 0; i < MFval.size(); i++)
	{
		opf << MFval[i] << "\t";
	}

	opf << "\n";
	opf.close();
	opf.clear();

	return 0;
}

/*
int solvn::Direct_task(vector<double>& MFval, vector<double>& MTval, CoeffStruct MyCoeffs, string path, vector<double> dP)
{
	ofstream opf;
	for (int i = 0; i < MFval.size(); i++)
	{
		MFval[i] = Receive_F(MTval[i], MyCoeffs, dP);
	}
	opf.open(path + "funcsres.txt", ios_base::app);

	for (int i = 0; i < MFval.size(); i++)
	{
		opf << MFval[i] << "\t";
	}
	opf << "\n";
	opf.close();
	opf.clear();

	return 0;
}
*/

int solvn::Inverse_task(string path)
{
	ofstream opf(path + "results.txt", ios_base::out | ios_base::trunc);
	vector<vector<double>> A(n_params, vector<double>(n_params));
	vector<double> b(n_params);
	vector<double> dP(n_params); // - diff of Parameters - Синоним для повышения читабельности кода и экономии памяти

	PrevCoeffs.Acoeff = ResCoeffs.Acoeff;
	PrevCoeffs.Bcoeff = ResCoeffs.Bcoeff;
	PrevCoeffs.Ccoeff = ResCoeffs.Ccoeff;

	cout << "Launch Inverse_task(...)" << endl;

	double func = 1E+30, prev_func, dif = 1E+30, alpha_loc;
	int iter_a;
	//if (n_params == Fval.size())
	if (true)
	{
		cout << "Pass check (n_params == Fval.size())" << endl;
		Direct_task(Fval, Tval, ResCoeffs, path);
		func = Get_func(Fval, Fval_true);
		cout << "func = " << func << endl;

		for (int i = 0; i < Fval.size(); i++)
			cout << "FVal[" << i << "] = " << Fval[i] << endl;
		for (int i = 0; i < Fval_true.size(); i++)
			cout << "FVal_true[" << i << "] = " << Fval_true[i] << endl;

		prev_func = func;

		opf << "iter 0: " << '\n';
		opf << "func = " << func << '\n' << "norma = " << 0 << '\n';
		opf << setprecision(17) << "A" << " = " << ResCoeffs.Acoeff << '\n';
		opf << setprecision(17) << "B" << " = " << ResCoeffs.Bcoeff << '\n';
		opf << setprecision(17) << "C" << " = " << ResCoeffs.Ccoeff << '\n';
		opf << "----------------------------------------------------------------" << endl;

		for (int iter = 1; iter < max_iter && func > eps_func && dif > eps_dif; iter++)
		{
			cout << "Pass for1" << endl;

			iter_a = 0;
			alpha_loc = alpha;
			//do
			//{
				// очистить СЛАУ
			Clear_SLAU(A, b);
			// собрать СЛАУ
			for (int q = 0; q < n_params; q++) // 3 == A B C == n_params
			{
				for (int p = 0; p < n_params; p++) // 3 == A B C
				{
					for (int i = 0; i < Fval.size(); i++)
					{
						//A[q][p] += dF_dt(Tval[i], ResCoeffs, q) * dF_dt(Tval[i], ResCoeffs, p);
						//F_dCoeff(double t, CoeffStruct MyCo-effs, int numCoeff):
						A[q][p] += F_dCoeff(Tval[i], ResCoeffs, q) * F_dCoeff(Tval[i], ResCoeffs, p);
						//A[q][p] += dF_dt(Tval[i], ResCoeffs, q) * dF_dt(Tval[i], ResCoeffs, p);
					}
				}
				for (int j = 0; j < Fval.size(); j++)
				{
					//b[q] -=  dF_dt(Tval[j], ResCoeffs, q) * (Fval[j] - Fval_true[j]);
					//F_dCoeff(double t, CoeffStruct MyCoeffs, int numCoeff):
					b[q] -= F_dCoeff(Tval[j], ResCoeffs, q) * (Fval[j] - Fval_true[j]);
					//b[q] -= dF_dt(Tval[j], ResCoeffs, q) * (Fval[j] - Fval_true[j]);
				}
				b[0] -= alpha * (PrevCoeffs.Acoeff - ResCoeffs.Acoeff);
				//b[1] -= alpha * (PrevCoeffs.Bcoeff - ResCoeffs.Bcoeff);
				//b[2] -= alpha * (PrevCoeffs.Ccoeff - ResCoeffs.Ccoeff);
			}

			// Добавляем регуляризацию
			for (int i = 0; i < n_params; i++)
				A[i][i] += alpha_loc;

			cout << "----------A[q][p]-----------" << endl;
			writeMatrix(A);
			cout << "----------------------------" << endl;
			cout << endl;

			cout << "------------b[q]------------" << endl;
			writeVector(b);
			cout << "----------------------------" << endl;
			cout << endl;

			// решить СЛАУ
			Gauss(A, b, dP, n_params);

			cout << "------------x[q]------------" << endl;
			writeVector(dP);
			cout << "----------------------------" << endl;
			cout << endl;
			cout << endl;

			dif = Get_norm(dP);

			PrevCoeffs.Acoeff = ResCoeffs.Acoeff;
			PrevCoeffs.Bcoeff = ResCoeffs.Bcoeff;
			PrevCoeffs.Ccoeff = ResCoeffs.Ccoeff;

			ResCoeffs.Acoeff += dP[0];
			ResCoeffs.Bcoeff += dP[1];
			//ResCoeffs.Ccoeff += dP[2];

			Direct_task(Fval, Tval, ResCoeffs, path);

			prev_func = func;
			func = Get_func(Fval, Fval_true);
			iter_a++;
			alpha_loc *= 10;

			//} while (prev_func < func || iter_a < 10);

			opf << "iter " << iter << ": " << '\n';
			opf << setprecision(17) << "A" << " = " << ResCoeffs.Acoeff << '\n';
			opf << setprecision(17) << "B" << " = " << ResCoeffs.Bcoeff << '\n';
			opf << setprecision(17) <<"C"  << " = " << ResCoeffs.Ccoeff << '\n';

			opf << "func = " << func << '\n' << "norma = " << dif << '\n';
			opf << "----------------------------------------------------------------" << endl;
		}
	}
	else
		cout << "n_params != MFvec.size()\n";
	opf.close();
	opf.clear();
	return 0;
}
//------------------------------------------


void solvn::Clear_SLAU(vector<vector<double>>& A, vector<double>& b)
{
	for (int q = 0; q < n_params; q++)
	{
		for (int p = 0; p < n_params; p++)
		{
			A[q][p] = 0;
		}
		b[q] = 0;
	}

}

int solvn::Gauss(vector<vector<double>>& A, vector<double>& b, vector<double>& x, int N)
{
	// Гаусс
   // приведение к треугольному виду
	double t;
	for (int k = 0; k < N - 1; k++)
	{
		for (int i = k + 1; i < N; i++)
		{
			t = A[i][k] / A[k][k];

			b[i] -= t * b[k];
			for (int j = k + 1; j < N; j++)
			{
				A[i][j] -= t * A[k][j];
			}
		}
	}
	b[N - 1] /= A[N - 1][N - 1];
	x[N - 1] = b[N - 1];
	// решение СЛАУ с треугольной матрицей
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
		{
			sum += A[k][j] * x[j];
		}
		x[k] = (b[k] - sum) / A[k][k];
	}

	return 0;
}

double solvn::Get_norm(vector<double>& point1)
{
	double r = 0;

	for (int i = 0; i < point1.size(); i++)
	{
		r += point1[i] * point1[i];
	}
	r = sqrt(r);
	return r;
}

void solvn::model_init()
{
	string path = "";
	Input_data(Fval, Fval_true, Tval, path);
}

void solvn::model_solve()
{
	string path = "";
	if (!flag_F)
		Direct_task(Fval_true, Tval, TrueCoeffs, path);
	Inverse_task(path);
}

int main()
{
	solvn Solver;

	Solver.model_init();
	Solver.model_solve();

	return 0;
}