#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

class solvn
{
public:
	solvn() {};
	~solvn() {};

	void model_init();
	void model_solve();

private:
	const double PI = 3.1415926535897932;
	int n_params = 3;
	int max_iter = 100;
	double eps_func = 1E-8;
	double eps_dif = 1E-14;
	double alpha = 1E-20;
	double sigma = 0;

	double Acoeff, Bcoeff, Ccoeff, Fq;
	double Acoeff_prev, Bcoeff_prev, Ccoeff_prev;
	struct FuncStruct
	{
		double Fval;
		double Fval_true;
		double Tval;
	};
	vector<FuncStruct> Fvec;

	struct Receiver
	{
		vector<double> M, N;
		double V, V_true;


		Receiver()
		{
			M.resize(3);
			N.resize(3);
			V = 0;
			V_true = 0;
		}
	};

	struct Source
	{
		vector<double> A, B;
		double I, I_prev;

		Source()
		{
			A.resize(3);
			B.resize(3);
			I = 0;
			I_prev = 0;
		}
	};

	vector<Receiver> MN;
	vector<Source> AB;

	double r(vector<double>& point1, vector<double>& point2);
	double Get_func(vector<Source>& AB, vector<Receiver>& MN);
	double dV_dI(vector<Source>& AB, Receiver& MN, double sigma, int id_I);
	int Get_V_on_receiver(vector<Source>& AB, Receiver& MN, double sigma);
	int Get_V_on_receiver(vector<Source>& AB, Receiver& MN, double sigma, vector<double> b);
	int Input_data(vector<Source>& AB, vector<Receiver>& MN, double& sigma, string path);
	int Direct_task(vector<Source>& AB, vector<Receiver>& MN, double sigma, string path);
	int Direct_task(vector<Source>& AB, vector<Receiver>& MN, double sigma, string path, vector<double> b);
	void Clear_SLAU(vector<vector<double>>& A, vector<double>& b);
	int Gauss(vector<vector<double>>& A, vector<double>& b, int N);
	double Get_norm(vector<double>& point1);
	int Inverse_task(vector<Source>& AB, vector<Receiver>& MN, double sigma, string path);
	//---------------------------------------------------------------------------------------------------------------
	double Get_func(vector<FuncStruct> MFvec);
	double dF_dt(double t, int numCoeff);
	int Get_F_on_receiver(FuncStruct &MF, double t);
	int Get_F_on_receiver(FuncStruct &MF, double t, vector<double> &b);
	int Input_data(vector<FuncStruct> &MFvec, string path);
	int Direct_task(vector<FuncStruct> &MFvec, string path);
	int Direct_task(vector<FuncStruct> &MFvec, string path, vector<double> b);
	int Inverse_task(vector<FuncStruct> &MFvec, string path);
};

double solvn::r(vector<double>& point1, vector<double>& point2)
{
	double r = 0;
	
	for (int i = 0; i < point1.size(); i++)
	{
		r += (point1[i] - point2[i]) * (point1[i] - point2[i]);
	}
	r = sqrt(r);
	return r;
}

double solvn::Get_func(vector<Source>& AB, vector<Receiver>& MN)
{
	double func = 0;
	for (int i = 0; i < AB.size(); i++)
	{
		func += (AB[i].I - AB[i].I_prev);
	}
	func *= alpha;
	for (int i = 0; i < MN.size(); i++)
	{
		func += (MN[i].V - MN[i].V_true) * (MN[i].V - MN[i].V_true) / (MN[i].V_true * MN[i].V_true);
	}

	return func;
}
//---------------------------------------------------------------------------------------------------------------
double solvn::Get_func(vector<FuncStruct> MFvec)
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


	for (int i = 0; i < MFvec.size(); i++)
	{
		func += (MFvec[i].Fval - MFvec[i].Fval_true) * (MFvec[i].Fval - MFvec[i].Fval_true) / (MFvec[i].Fval_true * MFvec[i].Fval_true);
	}

	return func;
}

double solvn::dF_dt(double t, int numCoeff)
{
	double dF = 0;

	switch (numCoeff)
	{
	case 0:
		dF = 2 * Acoeff;
		//dF = t * t;
		//dF = cos(Acoeff);
		break;
	case 1:
		dF = 2 * Bcoeff;
		//dF = t;
		//dF = -sin(Bcoeff);
		break;
	case 2:
		//dF = 1;
		break;
	default:
		break;
	}

	return dF;
}

int solvn::Get_F_on_receiver(FuncStruct &MF, double t)
{
	MF.Fval = (Acoeff* Acoeff +Bcoeff*Bcoeff);
	//MF.Fval = (t * t * Acoeff + t * Bcoeff + Ccoeff);
	//MF.Fval = (Acoeff*t*t);
	//MF.Fval = (sin(Acoeff)+cos(Bcoeff));

	return 0;
}

int solvn::Get_F_on_receiver(FuncStruct &MF, double t, vector<double> &b)
{
	MF.Fval = (Acoeff + b[0])* (Acoeff + b[0]) + (Bcoeff+ b[1])* (Bcoeff + b[1]);
	//MF.Fval = (t * t * (Acoeff+b[0]) + t * (Bcoeff+b[1]) + Ccoeff + b[2]);
	//MF.Fval = (t * t * (Acoeff + b[0]) );
	//MF.Fval = ((sin(Acoeff + b[0]))+ (cos(Bcoeff+ b[1])));
	return 0;
}

int solvn::Input_data(vector<FuncStruct> &MFvec, string path)
{
	ifstream inf;
	int tmp_int, flag_F, coeff, n_c;
	double tmp_dbl;


	inf.open(path + "tnet.txt", ios_base::in);
	if (inf)
	{
		inf >> tmp_int;
		//===================================================================================================
		MFvec.resize(tmp_int);
		//===================================================================================================
		for (int i = 0; i < tmp_int; i++)
			inf >> MFvec[i].Tval;
	}
	inf.close();
	inf.clear();

	inf.open(path + "params.txt", ios_base::in);
	if (inf)
	{
		inf >> n_params;
		inf >> max_iter;
		inf >> eps_func;
		inf >> eps_dif;

	}
	inf.close();
	inf.clear();

	inf.open(path + "fvalues.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "fvalues.txt";
		return 1;
	}
	inf >> flag_F;
	inf >> tmp_int;
	//===================================================================================================
	//MFvec.resize(tmp_int);
	//===================================================================================================
	if (flag_F) {
		for (int i = 0; i < tmp_int; i++)
		{
			inf >> MFvec[i].Fval_true;
		}
	}
	inf.close();
	inf.clear();

	inf.open(path + "startv.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "start_value.txt";
	}
	else
	{
		inf >> Acoeff;
		inf >> Bcoeff;
		inf >> Ccoeff;
	}
	inf.close();
	inf.clear();

	return 0;
}

int solvn::Direct_task(vector<FuncStruct> &MFvec, string path)
{
	ofstream opf;
	for (int i = 0; i < MFvec.size(); i++)
	{
		Get_F_on_receiver(MFvec[i], MFvec[i].Tval);
	}
	opf.open(path + "funcsres.txt", ios_base::app);
	for (int i = 0; i < MFvec.size(); i++)
	{
		opf << MFvec[i].Fval << "\t";
	}
	opf << "\n";
	opf.close();
	opf.clear();

	return 0;
}

int solvn::Direct_task(vector<FuncStruct>& MFvec, string path, vector<double> b)
{
	ofstream opf;
	for (int i = 0; i < MFvec.size(); i++)
	{
		Get_F_on_receiver(MFvec[i], MFvec[i].Tval, b);
	}
	opf.open(path + "funcsres.txt", ios_base::app);
	for (int i = 0; i < MFvec.size(); i++)
	{
		opf << MFvec[i].Fval << "\t";
	}
	opf << "\n";
	opf.close();
	opf.clear();

	return 0;
}

int solvn::Inverse_task(vector<FuncStruct>& MFvec, string path)
{
	ofstream opf(path + "results.txt", ios_base::out | ios_base::trunc);
	vector<vector<double>> A(n_params, vector<double>(n_params));
	vector<double> b(n_params);

	double func = 1E+30, prev_func, dif = 1E+30, alpha_loc;
	int iter_a;
	if (n_params == MFvec.size())
	{
		Direct_task(MFvec, path);
		func = Get_func(MFvec);
		prev_func = func;
		opf << "iter 0: " << '\n';
		opf << "func = " << func << '\n' << "norma = " << 0 << '\n';
		opf << "A" << " = " << Acoeff << '\n';
		opf << "B" << " = " << Bcoeff << '\n';
		opf << "C" << " = " << Ccoeff << '\n';
		opf << "----------------------------------------------------------------" << endl;

		for (int iter = 1; iter < max_iter && func > eps_func && dif > eps_dif; iter++)
		{
			iter_a = 0;
			alpha_loc = alpha;
			do
			{
				// очистить СЛАУ
				Clear_SLAU(A, b);
				// собрать СЛАУ
				for (int q = 0; q < n_params; q++) // 3 == A B C == n_params
				{
					for (int p = 0; p < n_params; p++) // 3 == A B C
					{
						for (int i = 0; i < MFvec.size(); i++)
						{
							A[q][p] += 1 / MFvec[i].Fval_true * dF_dt(MFvec[i].Tval, q) * dF_dt(MFvec[i].Tval, p); // уточнить по времени
						}
						if (p == q)
						{
							A[q][p] += alpha_loc;
						}
					}
					for (int j = 0; j < MFvec.size(); j++)
					{
						b[q] -= 1 / MFvec[j].Fval_true * dF_dt(MFvec[j].Tval, q) * (MFvec[j].Fval - MFvec[j].Fval_true);
					}
				}

				cout << "----------A[q][p]-----------" << endl;
				for (int q = 0; q < n_params; q++)
				{
					for (int p = 0; p < n_params; p++)
					{
						cout << A[q][p] << " ";
					}
					cout << endl;
				}
				cout << "------------b[q]------------" << endl;
				for (int q = 0; q < n_params; q++)
				{
					cout << b[q] << endl;
				}
				// решить СЛАУ
				Gauss(A, b, n_params);

				dif = Get_norm(b);
				Direct_task(MFvec, path, b);
				prev_func = func;
				func = Get_func(MFvec);
				iter_a++;
				alpha_loc *= 10;

				cout << "------------x[q]------------" << endl;
				for (int q = 0; q < n_params; q++)
				{
					cout << b[q] << endl;
				}

				Acoeff += b[0];;
				Bcoeff += b[1];;
				Ccoeff += b[2];;
			} while (prev_func < func || iter_a < 10);

			
			/*
			for (int i = 0; i < AB.size(); i++)
			{
				AB[i].I += b[i];
			}
			*/
			opf << "iter " << iter << ": " << '\n';
			/*
			for (int i = 0; i < AB.size(); i++)
			{
				opf << "I" << i + 1 << " = " << AB[i].I << '\n';
				AB[i].I_prev = AB[i].I;
			}
			*/
			opf << "A" << " = " << Acoeff << '\n';
			opf << "B" << " = " << Bcoeff << '\n';
			opf << "C" << " = " << Ccoeff << '\n';

			Acoeff_prev = Acoeff;
			Bcoeff_prev = Bcoeff;
			Ccoeff_prev = Ccoeff;

			opf << "func = " << func << '\n' << "norma = " << dif << '\n';
			opf << "----------------------------------------------------------------" << endl;
		}
	}
	else
		cout << "n_params != AB.size()\n";
	opf.close();
	opf.clear();
	return 0;
}
//---------------------------------------------------------------------------------------------------------------
double solvn::dV_dI(vector<Source>& AB, Receiver& MN, double sigma, int id_I)
{
	return (((1. / r(MN.M, AB[id_I].B) - 1. / r(MN.M, AB[id_I].A)) -
		(1. / r(MN.N, AB[id_I].B) - 1. / r(MN.N, AB[id_I].A)))) / (2 * PI * sigma);
}

int solvn::Get_V_on_receiver(vector<Source>& AB, Receiver& MN, double sigma)
{
	double V = 0;

	for (int i = 0; i < AB.size(); i++)
	{
		V += AB[i].I * ((1. / r(MN.M, AB[i].B) - 1. / r(MN.M, AB[i].A)) -
			(1. / r(MN.N, AB[i].B) - 1. / r(MN.N, AB[i].A)));
	}
	V /= 2 * PI * sigma;
	MN.V = V;
	return 0;
}
int solvn::Get_V_on_receiver(vector<Source>& AB, Receiver& MN, double sigma, vector<double> b)
{
	double V = 0;

	for (int i = 0; i < AB.size(); i++)
	{
		V += (AB[i].I + b[i]) * ((1. / r(MN.M, AB[i].B) - 1. / r(MN.M, AB[i].A)) -
			(1. / r(MN.N, AB[i].B) - 1. / r(MN.N, AB[i].A)));
	}
	V /= 2 * PI * sigma;
	MN.V = V;
	return 0;
}

int solvn::Input_data(vector<Source>& AB, vector<Receiver>& MN, double& sigma, string path)
{
	ifstream inf;
	int tmp_int, flag_I, flag_V, n, net, coeff, n_c;
	double tmp_dbl;


	inf.open(path + "t_net.txt", ios_base::in);
	if (inf)
	{
		inf >> n;
		for (int i = 0; i < n; i++)
			inf >> net;
	}
	inf.close();
	inf.clear();

	inf.open(path + "coeff.txt", ios_base::in);
	if (inf)
	{
		inf >> n_c;
		for (int i = 0; i < n; i++)
			inf >> coeff;
	}
	inf.close();
	inf.clear();


	inf.open(path + "params.txt", ios_base::in);
	if (inf)
	{
		inf >> n_params;
		inf >> max_iter;
		inf >> eps_func;
		inf >> eps_dif;

	}
	inf.close();
	inf.clear();

	inf.open(path + "recivers.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "recivers.txt";
		return 1;
	}
	inf >> flag_V;
	inf >> tmp_int;
	MN.resize(tmp_int);
	for (int i = 0; i < tmp_int; i++)
	{
		for (int j = 0; j < MN[i].M.size(); j++)
			inf >> MN[i].M[j];
		for (int j = 0; j < MN[i].N.size(); j++)
			inf >> MN[i].N[j];
		if (flag_V)
			inf >> MN[i].V_true;
	}
	inf.close();
	inf.clear();

	inf.open(path + "sources.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "sources.txt";
		return 1;
	}
	inf >> flag_I;
	inf >> tmp_int;
	AB.resize(tmp_int);
	for (int i = 0; i < tmp_int; i++)
	{
		for (int j = 0; j < AB[i].A.size(); j++)
			inf >> AB[i].A[j];
		for (int j = 0; j < AB[i].B.size(); j++)
			inf >> AB[i].B[j];
		if (flag_I)
			inf >> AB[i].I;
	}
	inf.close();
	inf.clear();

	inf.open(path + "startv.txt", ios_base::in);
	if (!inf)
	{
		cout << "Can't open " << path + "start_value.txt";
	}
	else
	{
		if (flag_I)
			inf >> sigma;
		else
		{
			for (int i = 0; i < AB.size(); i++)
			{
				inf >> AB[i].I;
				AB[i].I_prev = AB[i].I;
			}
		}
	}
	inf.close();
	inf.clear();

	inf.open(path + "sigma.txt", ios_base::in);
	if (inf)
	{
		inf >> tmp_dbl;
		if (tmp_dbl < 0)
			sigma = 0;
		else sigma = tmp_dbl;
	}
	inf.close();
	inf.clear();

	inf.open(path + "volts.txt", ios_base::out | ios_base::trunc);
	inf.close();
	inf.clear();

	return 0;
}

int solvn::Direct_task(vector<Source>& AB, vector<Receiver>& MN, double sigma, string path)
{
	ofstream opf;
	for (int i = 0; i < MN.size(); i++)
	{
		Get_V_on_receiver(AB, MN[i], sigma);
	}
	opf.open(path + "volts.txt", ios_base::app);
	for (int i = 0; i < MN.size(); i++)
	{
		opf << MN[i].V << "\t";
	}
	opf << "\n";
	opf.close();
	opf.clear();

	return 0;
}

int solvn::Direct_task(vector<Source>& AB, vector<Receiver>& MN, double sigma, string path, vector<double> b)
{
	ofstream opf;
	for (int i = 0; i < MN.size(); i++)
	{
		Get_V_on_receiver(AB, MN[i], sigma, b);
	}
	opf.open(path + "volts.txt", ios_base::app);
	for (int i = 0; i < MN.size(); i++)
	{
		opf << MN[i].V << "\t";
	}
	opf << "\n";
	opf.close();
	opf.clear();

	return 0;
}


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

int solvn::Gauss(vector<vector<double>>& A, vector<double>& b, int N)
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
	// решение СЛАУ с треугольной матрицей
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
		{
			sum += A[k][j] * b[j];
		}
		b[k] = (b[k] - sum) / A[k][k];
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

int solvn::Inverse_task(vector<Source>& AB, vector<Receiver>& MN, double sigma, string path)
{
	ofstream opf(path + "results.txt", ios_base::out | ios_base::trunc);
	vector<vector<double>> A(n_params, vector<double>(n_params));
	vector<double> b(n_params);

	double func = 1E+30, prev_func, dif = 1E+30, alpha_loc;
	int iter_a;
	if (n_params == AB.size())
	{
		Direct_task(AB, MN, sigma, path);
		func = Get_func(AB, MN);
		prev_func = func;
		opf << "iter 0: " << '\n';
		for (int i = 0; i < AB.size(); i++)
		{
			opf << "I" << i + 1 << " = " << AB[i].I << '\n';
		}
		opf << "func = " << func << '\n' << "norma = " << 0 << '\n';

		for (int iter = 1; iter < max_iter && func > eps_func && dif > eps_dif; iter++)
		{
			iter_a = 0;
			alpha_loc = alpha;
			do
			{
				// очистить СЛАУ
				Clear_SLAU(A, b);
				// собрать СЛАУ
				for (int q = 0; q < n_params; q++)
				{
					for (int p = 0; p < n_params; p++)
					{
						for (int i = 0; i < MN.size(); i++)
						{
							A[q][p] += 1 / MN[i].V_true * dV_dI(AB, MN[i], sigma, q) * dV_dI(AB, MN[i], sigma, p);
						}
						if (p == q)
						{
							A[q][p] += alpha_loc;
						}
					}
					for (int j = 0; j < MN.size(); j++)
					{
						b[q] -= 1 / MN[j].V_true * dV_dI(AB, MN[j], sigma, q) * (MN[j].V - MN[j].V_true);
					}
				}
				// решить СЛАУ
				Gauss(A, b, n_params);

				dif = Get_norm(b);
				Direct_task(AB, MN, sigma, path, b);
				prev_func = func;
				func = Get_func(AB, MN);
				iter_a++;
				alpha_loc *= 10;
			} while (prev_func < func || iter_a < 10);

			for (int i = 0; i < AB.size(); i++)
			{
				AB[i].I += b[i];
			}

			opf << "iter " << iter << ": " << '\n';
			for (int i = 0; i < AB.size(); i++)
			{
				opf << "I" << i + 1 << " = " << AB[i].I << '\n';
				AB[i].I_prev = AB[i].I;
			}
			opf << "func = " << func << '\n' << "norma = " << dif << '\n';
		}
	}
	else
		cout << "n_params != AB.size()\n";
	opf.close();
	opf.clear();
	return 0;
}

void solvn::model_init()
{
	string path = "";
	Input_data(Fvec, path);
}

void solvn::model_solve()
{
	string path = "";
	//Direct_task(AB, MN, sigma, path);
	//Inverse_task(AB, MN, sigma, path);
	//Direct_task(Fvec, path);
	Inverse_task(Fvec, path);
}

int main()
{
	solvn Solver;

	Solver.model_init();
	Solver.model_solve();

	return 0;
}