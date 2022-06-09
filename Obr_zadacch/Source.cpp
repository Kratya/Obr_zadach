#define _USE_MATH_DEFINES
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <filesystem>
#include <windows.h>
#include <tchar.h>
#include <thread>

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

vector<double> operator+(vector<double> c1, vector<double> c2)
{
	vector<double> vTmpCoeff;
	const int n = c1.size();

	vTmpCoeff.resize(n);

	for (int i = 0; i < n; i++)
	{
		vTmpCoeff[i] = c1[i] + c2[i];
	}

	return vTmpCoeff;
}

vector<double> operator+(vector<double> c1, double aC)
{
	vector<double> vTmpCoeff;
	const int n = c1.size();

	vTmpCoeff.resize(n);

	for (int i = 0; i < n; i++)
	{
		vTmpCoeff[i] = c1[i] + aC;
	}

	return vTmpCoeff;
}

vector<double> operator*(vector<double> c1, double aC)
{
	vector<double> vTmpCoeff;
	const int n = c1.size();

	vTmpCoeff.resize(n);

	for (int i = 0; i < n; i++)
	{
		vTmpCoeff[i] = c1[i] * aC;
	}

	return vTmpCoeff;
}

int makeDir(string dir)
{
	error_code err;
	if (std::filesystem::exists(dir))
		return std::filesystem::is_directory(std::filesystem::status(dir));
	else
		return std::filesystem::create_directories(dir, err);
}

class paramWell
{
public:

	int num_well; // Номер скважины
	string name;
	int n_pump; //
	string type_name;
	int n_times;
	vector<int> last_day;
	vector<double> prev_theta; // Для расчетов предыдущая тета
	vector<double> theta;
	vector<int> p1;
	vector<int> p2;
	vector<int> p3;
	int stat; // Тип скважины

	paramWell() {};
	paramWell(int _num_well, string _name, int _n_pump, string _type_name, int _n_times,
		vector<int> _last_day, vector<double> _theta, vector<int> _p1, vector<int> _p2, vector<int> _p3)
	{
		this->num_well = _num_well;
		this->name = _name;
		this->n_pump = _n_pump;
		this->type_name = _type_name;
		this->n_times = _n_times;
		this->last_day = _last_day;
		this->theta = _theta;
		this->p1 = _p1;
		this->p2 = _p2;
		this->p3 = _p3;

		prev_theta.resize(_n_times);
	};
	~paramWell() {};
};

int read_well_param(vector<paramWell>& ourWells, int& num_all_well, string path)
{
	string r_name;
	int r_n_pump;
	string r_type_name;
	int r_n_times;
	vector<int> r_last_day;
	vector<double> r_theta;
	vector<int> r_p1;
	vector<int> r_p2;
	vector<int> r_p3;

	ifstream in(path + "/wellsCONS.txt");
	if (in.is_open())
	{
		in >> num_all_well;
		if (num_all_well > 0)
		{
			for (int i = 0; i < num_all_well; i++)
			{
				in >> r_name;
				in >> r_n_pump;
				in >> r_type_name;
				in >> r_n_times;
				if (r_n_times > 0)
				{
					r_last_day.resize(r_n_times);
					r_theta.resize(r_n_times);
					r_p1.resize(r_n_times);
					r_p2.resize(r_n_times);
					r_p3.resize(r_n_times);

					for (int j = 0; j < r_n_times; j++)
					{
						in >> r_last_day[j];
						in >> r_theta[j];
						in >> r_p1[j];
						in >> r_p2[j];
						in >> r_p3[j];
					}
				}

				ourWells.push_back(paramWell(i + 1, r_name, r_n_pump, r_type_name, r_n_times, r_last_day, r_theta, r_p1, r_p2, r_p3));
			}
		}
	}
	else
	{
		cout << "err open wellsCONS.txt" << endl;
		return 1;
	}

	return 0;
}

int write_well_param(vector<paramWell> ourWells, int num_all_well, string path)
{
	ofstream out(path + "/wellsCONS.txt");
	if (out.is_open())
	{
		if (num_all_well > 0)
		{
			out << num_all_well << "\n";
			for (int i = 0; i < num_all_well; i++)
			{
				out << ourWells[i].name << "\n";
				out << ourWells[i].n_pump << "\n";
				out << ourWells[i].type_name << "\n";
				out << ourWells[i].n_times << "\n";

				for (int j = 0; j < ourWells[i].n_times; j++)
				{
					out << ourWells[i].last_day[j] << "\n";
					out << ourWells[i].theta[j] << "\n";
					out << ourWells[i].p1[j] << "\n";
					out << ourWells[i].p2[j] << "\n";
					out << ourWells[i].p3[j] << "\n";
				}
			}
		}
	}
	else
	{
		cout << "err open wellsCONS.txt" << endl;
		return 1;
	}

	return 0;
}

class solvn
{
public:
	solvn() {};
	~solvn() {};

	/// <summary>
	/// Инициализация модели (Чтение входных данных)
	/// </summary>
	void model_init();

	/// <summary>
	/// Решение обратной задачи
	/// </summary>
	void model_solve();
	void checkGauss();

private:
	int max_iter = 100;
	double eps_func = 1E-8;
	double eps_dif = 1E-14;
	double alpha = 1E-8;
	double delta = 0.005;
	double dThetaProc = 0.1;
	double absThetaMax = 100;
	double vStar = -100;
	//string name_model;
	int is_polymer_flag = 0;

	// Переменные для отладки
	int dthetcount = 0;
	int dfunccount = 0;
	int dvalscount = 0;
	int d_iter = 0;

	int alphas_flag;
	vector<double>	vAlphas;
	vector<int> vFxThetas;

	int nWells; // количество скважин
	int a_dim; // размерность, зависящая от кол-ва мощностей
	vector<int> numWell;
	vector<int> numTheta;

	vector<paramWell> ourWells;
	vector<vector<double>> fWellsAll, TWellsAll;
	string pathComlex = "Filtr3D.exe";
	string pathDatas = "./output";
	string pathWellsConf = "./properties";

	/// <summary>
	/// Зануление матрицы СЛАУ и правой части
	/// </summary>
	/// <param name="A"> - Матрица СЛАУ</param>
	/// <param name="b"> - Вектор правой части</param>
	void clear_SLAU(vector<vector<double>>& A, vector<double>& b);

	/// <summary>
	/// Решение СЛАУ методом Гаусса
	/// </summary>
	/// <param name="A"> - Матрица СЛАУ</param>
	/// <param name="b"> - Вектор правой части</param>
	/// <param name="x"> - Вектор для записи решения</param>
	/// <param name="N"> - Размерность СЛАУ</param>
	/// <returns></returns>
	int gauss_plus(vector<vector<double>>& A, vector<double>& b, vector<double>& x, int N);
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

	/// <summary>
	/// Аналитический расчёт производных целевой функции по искомым коэффициентам
	/// </summary>
	/// <param name="t"> - Точка, значение производной функции по параметру в которой нужно рассчитать</param>
	/// <param name="MyCoeffs"> - Коэффициенты, по которым дифференциируется функция</param>
	/// <param name="numCoeff"> - Номер параметра по которому ищеться производная</param>
	/// <returns></returns>

	/// <summary>
	/// Численный расчёт производных целевой функции по искомым коэффициентам
	/// </summary>
	/// <param name="t"> - Точка, значение производной функции по параметру в которой нужно рассчитать</param>
	/// <param name="MyCoeffs"> - Коэффициенты, по которым дифференциируется функция</param>
	/// <param name="numCoeff"> - Номер параметра по которому ищеться производная</param>
	/// <returns></returns>

	/// <summary>
	/// Целевая функция
	/// </summary>
	/// <param name="t"> - Значение аргумента, по которому необходимо вычислить значение функции</param>
	/// <param name="MyCoeffs"> - Коэффициенты для вычисления целевой функции</param>
	/// <returns></returns>
	/// <summary>
	/// Вроде как тоже целевая функция
	/// </summary>
	/// <param name="t"> - Значение аргумента в котором необходимо узнать значение функции</param>
	/// <param name="MyCoeffs"> - Значение коэффициентов целевой функции</param>
	/// <returns> - Значение целевой функции</returns>

	/// <summary>
	/// Вроде как тоже целевая функция с зависимостью от решения
	/// </summary>
	/// <param name="t"> - Значение аргумента в котором необходимо узнать значение функции</param>
	/// <param name="MyCoeffs"> - Значение коэффициентов целевой функции</param>
	/// <returns> - Значение целевой функции</returns>
	//Функции на чтение
	int load_data(vector<vector<double>>& MFVal, vector<vector<double>>& MTval, string pathProg, string pathData);
	int fLineParams(string path);
	int fLineVals(vector<double>& MFval_true, vector<double>& MTval, string path, string tFile);
	int fChkFxThetas(string path);
	int genFList(vector<string>& fWOList, vector<string>& fPolyList);
	int read_all_vals(vector<vector<double>>& AllFVals, vector<vector<double>>& AllTVals, string path, vector<string> woFiles, vector<string> polyFiles);
	//------------------------------------------------------------------------------------------------------------------
	int opti_task(string path);
	// Расчет производных, интеграл и функционал
	double trapez(vector<double> tVal, vector<double> fVal);
	double deriv(double fVal, double fValNew, double MFThetas);
	double functional(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff);
	double functional(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff, 
					  vector<string> fWOList, vector<string> fPolyList, int iter);
	int wdebugfunc(double dfunc, int d_iter);
	int wdebugtheta(vector<paramWell> fourWells, int d_iter);
	int wdebugvals(int d_iter);
	int checkMaxMinTheta(vector<paramWell>& fourWells, vector<double> deltaP);
	void prep_deriv_folder();
	void prep_data_deriv(vector<double>& vecFunc, int i, vector<string> fWOFiles, vector<string> fPolyFiles);

};

//----------------------------------------
double solvn::deriv(double fVal, double fValNew, double MFThetas)
{
	double result = 0;

	result = (fValNew - fVal) / abs(MFThetas * delta);

	return result;
}

double solvn::trapez(vector<double> tVal, vector<double> fVal)
{
	const int n = fVal.size();
	double result = 0, x1, x2;

	//cout << "calc trapez:" << endl;

	for (int h = 0; h < n - 1; h++)
	{
		x1 = tVal[h];
		x2 = tVal[h + 1];
		result += 0.5 * (x2 - x1) * (fVal[h] + fVal[h + 1]);
		cout << "x1:" << x1 << " x2:" << x2 << "fval[" << h << "]:" << fVal[h] << "fval[" << h+1 << "]:" << fVal[h + 1] << endl;

	}

	cout << "result = " << result << endl;

	return result;

}

double solvn::functional(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff)
{
	int log = 0;
	double result = 0;
	vector<double> part_res;
	//cout << "start calc functional" << endl;
	//for (int i = 0; i < 3 * nWells; i++)
	//{
	//	cout << "Coeff[" << i << "] = " << MFCoeff[i] << endl;
	//}

	for (int i = 0; i < nWells; i++)
	{
		for (int j = 0; j < ourWells[i].n_times; j++) {
			if (ourWells[i].stat)
			{
				for (int j = 0; j < ValWellsAll[2 * i + 1].size(); j++)
				{
					ValWellsAll[2 * i + 1][j] = (ValWellsAll[2 * i + 1][j] + vStar);
					ValWellsAll[2 * i + 1][j] = abs(ValWellsAll[2 * i + 1][j]);
				}
			}
		}
	}

	int num_values;
	if (is_polymer_flag)
	{
		num_values = 3 * nWells;
	}
	else
	{
		num_values = 2 * nWells;
	}


	part_res.resize(num_values);
	for (int i = 0; i < num_values; i++)
	{
		part_res[i] = trapez(TWellsAll[i], ValWellsAll[i]);
	}

	// Подбор альф
	if (alphas_flag == 2)
	{
		for (int i = 0; i < num_values; i++)
			MFCoeff[i] = 1;
		alphas_flag = 0;
	}

	for (int i = 0; i < num_values; i++)
	{
		result += pow(MFCoeff[i] * part_res[i], 2);
	}
	cout << "functional result = " << result << endl;

	// debug thetas ------------------------------>
	//makeDir("./NewDebug/debug_i" + to_string(d_iter));
	//ofstream dfunc("./NewDebug/debug_i" + to_string(d_iter) + "/" + "dfunc_i" + to_string(d_iter) + "_n" + to_string(dfunccount) + ".txt");
	//dfunc << result << endl;
	//dfunccount++;	
	//dfunc.close();
	// <-------------------------------------------

	return result;
}

double solvn::functional(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff,
						 vector<string> fWOList, vector<string> fPolyList, int iter)
{
	int log = 0;
	double result = 0;
	vector<double> part_res;
	//cout << "start calc functional" << endl;
	//for (int i = 0; i < 3 * nWells; i++)
	//{
	//	cout << "Coeff[" << i << "] = " << MFCoeff[i] << endl;
	//}
	
	for (int i = 0; i < nWells; i++) 
	{
		for (int j = 0; j < ourWells[i].n_times; j++) {
			if (ourWells[i].stat)
			{
				for (int j = 0; j < ValWellsAll[2 * i + 1].size(); j++)
				{
					ValWellsAll[2 * i + 1][j] = (ValWellsAll[2 * i + 1][j] + vStar);
					ValWellsAll[2 * i + 1][j] = abs(ValWellsAll[2 * i + 1][j]);
				}
			}
		}
	}

	int num_values;
	if (is_polymer_flag)
	{
		num_values = 3 * nWells;
	}
	else
	{
		num_values = 2 * nWells;
	}

	part_res.resize(num_values);
	for (int i = 0; i < num_values; i++)
	{
		part_res[i] = trapez(TWellsAll[i], ValWellsAll[i]);
	}

	// Подбор альф
	if (alphas_flag == 2)
	{
		for (int i = 0; i < num_values; i++)
			MFCoeff[i] = 1;
		alphas_flag = 0;
	}

	// Значение интеграла каждого слагаемого записываем в файл:
	makeDir("./debug_data/");
	ofstream out;
	double integr_res = 0;
	for (int i = 0; i < nWells; i++)
	{
		integr_res += pow(MFCoeff[2 * i] * part_res[2 * i], 2);
	}
	out.open("./debug_data/integr_s1.txt", std::ios::app);
	out << iter << ". " << integr_res << endl;
	out.close();
	integr_res = 0;

	for (int i = 0; i < nWells; i++)
	{
		integr_res += pow(MFCoeff[2 * i + 1] * part_res[2 * i + 1], 2);
	}
	out.open("./debug_data/integr_s2.txt", std::ios::app);
	out << iter << ". " << integr_res << endl;
	out.close();
	integr_res = 0;

	if (is_polymer_flag) {
		int j = 2 * nWells;
		for (int i = 0; i < nWells; i++)
		{
			integr_res += pow(MFCoeff[i + j] * part_res[i + j], 2);
		}
		out.open("./debug_data/integr_c3.txt", std::ios::app);
		out << iter << ". " << integr_res << endl;
		out.close();
		integr_res = 0;
	}
	
	for (int i = 0; i < 2 * nWells; i++)
	{
		out.open("./debug_data/integr_" + fWOList[i], std::ios::app);
		out << iter << ". " << pow(MFCoeff[i] * part_res[i], 2) << endl;
		out.close();
	}

	if (is_polymer_flag) {
		int j = 2 * nWells;
		for (int i = j; i < j + nWells; i++)
		{
			out.open("./debug_data/integr_" + fPolyList[i - j], std::ios::app);
			out << iter << ". " << pow(MFCoeff[i] * part_res[i], 2) << endl;
			out.close();
		}
	}
	

	for (int i = 0; i < num_values; i++)
	{
		result += pow(MFCoeff[i] * part_res[i], 2);
	}
	cout << "functional result = " << result << endl;

	// debug thetas ------------------------------>
	//makeDir("./NewDebug/debug_i" + to_string(d_iter));
	//ofstream dfunc("./NewDebug/debug_i" + to_string(d_iter) + "/" + "dfunc_i" + to_string(d_iter) + "_n" + to_string(dfunccount) + ".txt");
	//dfunc << result << endl;
	//dfunccount++;	
	//dfunc.close();
	// <-------------------------------------------

	return result;
}

void solvn::prep_deriv_folder()
{
	makeDir("./debug_data/");
	for (int i = 0; i < a_dim; i++)
	{
		makeDir("./debug_data/deriv" + to_string(i) + "/");
		//std::filesystem::copy("./", "./debug_data/deriv" + to_string(i), std::filesystem::copy_options::recursive); // "./fcopyf/" ???
		std::filesystem::copy("./graph", "./debug_data/deriv" + to_string(i) + "/graph", std::filesystem::copy_options::recursive);
		std::filesystem::copy("./mesh", "./debug_data/deriv" + to_string(i) + "/mesh", std::filesystem::copy_options::recursive);
		std::filesystem::copy("./output_balance_Q", "./debug_data/deriv" + to_string(i) + "/output_balance_Q", std::filesystem::copy_options::recursive);
		std::filesystem::copy("./output3D", "./debug_data/deriv" + to_string(i) + "/output3D", std::filesystem::copy_options::recursive);
		std::filesystem::copy("./properties", "./debug_data/deriv" + to_string(i) + "/properties", std::filesystem::copy_options::recursive);
		std::filesystem::copy("./tables", "./debug_data/deriv" + to_string(i) + "/tables", std::filesystem::copy_options::recursive);
		std::filesystem::copy("./alignment.txt", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./concrt140.dll", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./DllInitAndDrawText.dll", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./FilesConfig.cfg", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./filtr.cfg", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./filtr.log", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./Filtr3D.exe", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./libiomp5md.dll", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./lines.txt", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./mkl_avx.dll", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./mkl_core.dll", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./mkl_def.dll", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./mkl_intel_thread.dll", "./debug_data/deriv" + to_string(i));
		//std::filesystem::copy("./Pressure Change.txt", "./debug_data/deriv" + to_string(i));
		std::filesystem::copy("./time_points.txt", "./debug_data/deriv" + to_string(i));
		//std::filesystem::copy("./" + name_model, "./debug_data/deriv" + to_string(i)); // ????
		//std::filesystem::copy("./output", "./debug_data/deriv" + to_string(i) + "/output", std::filesystem::copy_options::recursive);// ????
	}
}

void solvn::prep_data_deriv(vector<double>& vecFunc, int i, vector<string> fWOFiles, vector<string> fPolyFiles)
{
	vector<vector<double>> fWellsAllNew, TWellsAllNew;
	vector<paramWell> plDeltaWells = ourWells;
	double tTheta;

	std::filesystem::path cwd = std::filesystem::current_path();
	cout << cwd.string() << endl;

	int num_values;
	if (is_polymer_flag)
	{
		num_values = 3 * nWells;
	}
	else
	{
		num_values = 2 * nWells;
	}

	fWellsAllNew.resize(num_values);
	TWellsAllNew.resize(num_values);
	//for (int i = 0; i < a_dim; i++)
	//{ 

	tTheta = plDeltaWells[numWell[i]].theta[numTheta[i]];
	plDeltaWells[numWell[i]].theta[numTheta[i]] = plDeltaWells[numWell[i]].theta[numTheta[i]] +
		abs(plDeltaWells[numWell[i]].theta[numTheta[i]] * delta);
	cout << "write1 deriv well params" << endl;
	write_well_param(plDeltaWells, nWells, "./debug_data/deriv" + to_string(i) + "/properties/");
	cout << "launch deriv filtr3d" << endl;
	//system(("\\debug_data\\deriv" + to_string(i) + "\\" + pathComlex).c_str());
	//-------------------------------------------------------------------------------------------------------------------------
	// ("\\debug_data\\deriv" + to_string(i) + "\\" + pathComlex).c_str()
	//string tpath = "\\debug_data\\deriv" + to_string(i) + "\\" + pathComlex;
	//std::wstring stemp = std::wstring(tpath.begin(), tpath.end());
	//system(("c: && cd C:/Users/krave/OneDrive/Рабочий стол/МОЙ ДИПЛОМ/ИССЛЕДОВАНИЯ/тест с полимером и другой нефтенасыщ/debug_data/deriv" + to_string(i) + "/ &&" + pathComlex).c_str());
	system(("c: && cd " + cwd.string() + "\\debug_data\\deriv" + to_string(i) + " \\ && " + pathComlex).c_str());
	//ShellExecute(NULL, L"open", stemp.c_str(), NULL, NULL, SW_SHOWDEFAULT);
	//-------------------------------------------------------------------------------------------------------------------------

	cout << "read_all_vals deriv" << endl;
	read_all_vals(fWellsAllNew, TWellsAllNew, "./debug_data/deriv" + to_string(i) + "/output", fWOFiles, fPolyFiles);
	vecFunc[i] = functional(fWellsAllNew, TWellsAllNew, vAlphas);

	plDeltaWells[numWell[i]].theta[numTheta[i]] = tTheta;
	cout << "write2 deriv well params" << endl;
	write_well_param(plDeltaWells, nWells, "./debug_data/deriv" + to_string(i) + "/properties/");
	//}
}

int solvn::wdebugfunc(double dfunc, int d_iter)
{
	makeDir("./debug_data/");
	ofstream f_func("./debug_data/debug_func.txt", std::ios::app);
	f_func << d_iter << " " << dfunc << '\n';

	return 0;
}

int solvn::wdebugtheta(vector<paramWell> fourWells, int d_iter)
{
	int sz = fourWells.size();
	makeDir("./debug_data/");
	ofstream f_theta;
	for (int i = 0; i < sz; i++)
	{
		f_theta.open("./debug_data/debug_theta_w" + to_string(i + 1) + ".txt", std::ios::app);

		for (int j = 0; j < fourWells[i].n_times; j++)
		{
			f_theta << d_iter << " " << fourWells[i].theta[j] << " ";
		}
		f_theta << '\n';
		f_theta.close();
	}

	return 0;
}

int solvn::wdebugvals(int d_iter)
{
	makeDir("./debug_data/");
	ofstream deb_vals; //f_vals;
	ifstream f_vals;
	deb_vals.open("./debug_data/debug_vals.txt", std::ios::app);
	/*if (d_iter == 0)
	{
		deb_vals << "iter: ";
		for (int i = 0; i < nWells; i++)
		{
			deb_vals << ourWells[i].name + "_s1 ";
			deb_vals << ourWells[i].name + "_s2 ";
		}

		for (int i = 0; i < nWells; i++)
		{
			deb_vals << ourWells[i].name + "_c3 ";
		}
		deb_vals << "\n";
	}*/
	f_vals.open("./output/Well1_s2.txt");
	deb_vals << d_iter << ": ";
	double TVal, FVal;

	while (f_vals >> TVal >> FVal);
	deb_vals << FVal << " | ";
	deb_vals << "\n";

	f_vals.close();
	deb_vals.close();
	return 0;
}

int solvn::fLineParams(string path)
{
	ifstream in;

	in.open(path + "params.txt", ios_base::in);
	if (in.is_open())
	{
		in >> max_iter; // Максимальное количество итераций при решении обратной задачи
		in >> eps_func; // Минимальное значение функционала
		in >> eps_dif; // Минимальоне значение нормы вектора b
		in >> delta;
		in >> vStar;
		//in >> name_model;
		in >> is_polymer_flag;

		in.close();

		in.close();
		in.clear();
	}
	else
		return 1;

	return 0;
}

int solvn::fChkFxThetas(string path)
{
	ifstream in;
	int flag = 0;
	vFxThetas.resize(a_dim);

	for (int i = 0; i < a_dim; i++)
		vFxThetas[i] = 0;

	in.open(path + "fxthetas.txt", ios_base::in);
	if (in.is_open())
	{
		in >> flag;
		if (flag)
		{
			for (int i = 0; i < a_dim; i++)
			{
				in >> vFxThetas[i];
			}
		}

		in.close();
		in.clear();

		for (int i = 0; i < a_dim; i++)
		{
			cout << "vFxThetas[" << i << "]" << vFxThetas[i] << endl;
		}
	}
	else
	{
		cout << "err read fxthetas.txt" << endl;
		system("pause");
		return -1;
	}

	return flag;
}

int solvn::fLineVals(vector<double>& MFval_true, vector<double>& MTval, string path, string tFile)
{
	double Fval, Tval;
	ifstream in(path+ "/" + tFile);

	if (in.is_open())
	{
		if (MFval_true.size() != 0)
		{
			MFval_true.clear();
			MFval_true.shrink_to_fit();
		}
		if (MTval.size() != 0)
		{
			MTval.clear();
			MTval.shrink_to_fit();
		}

		while (in >> Tval >> Fval)
		{
			MFval_true.push_back(Fval);
			MTval.push_back(Tval);
		}
		
		in.close();
		in.clear();
	}
	else
	{
		cout << "err read vals" << endl;
		return 1;
	}
	// debug info --->
	//makeDir("./NewDebug/debug_i" + to_string(d_iter));
	//ofstream dvals("./NewDebug/debug_i" + to_string(d_iter) + "/" + tFile + "_i" + to_string(d_iter) + "_n" + to_string(dvalscount) + ".txt");
	//int numvals = MFval_true.size();
	//for (int i = 0; i < numvals; i++)
	//{
	//	dvals << MTval[i] << " : " << MFval_true[i] << "\n";
	//}
	//dvalscount++;
	//dvals.close();
	// <--------------

	return 0;
}

int solvn::genFList(vector<string>& fWOList, vector<string>& fPolyList)
{
	string filename;

	if (fWOList.size() != 0)
	{
		fWOList.clear();
		fWOList.shrink_to_fit();
	}

	if (fPolyList.size() != 0)
	{
		fPolyList.clear();
		fPolyList.shrink_to_fit();
	}

	for (int i = 0; i < nWells; i++)
	{
		filename = "avg_m_" + ourWells[i].name + "_s1.txt";
		fWOList.push_back(filename);
		filename = "avg_m_" + ourWells[i].name + "_s2.txt";
		fWOList.push_back(filename);

		filename = "m_M_c3_s1_" + ourWells[i].name + ".txt";
		fPolyList.push_back(filename);
	}

	return 0;
}

int solvn::read_all_vals(vector<vector<double>>& AllFVals, vector<vector<double>>& AllTVals, string path, vector<string> woFiles, vector<string> polyFiles)
{
	int j = 2 * nWells;

	// Считываем файлы с водой и нефтью
	for (int i = 0; i < j; i++)
	{
		fLineVals(AllFVals[i], AllTVals[i], path, woFiles[i]);
	}

	// Считываем файлы с полиномом
	if (is_polymer_flag) {
		for (int i = j; i < j + nWells; i++)
		{
			fLineVals(AllFVals[i], AllTVals[i], path, polyFiles[i - j]);
		}
	}

	return 0;
}

int solvn::load_data(vector<vector<double>>& MFVal, vector<vector<double>>& MTval, string pathProg, string pathData)
{
	vector<string> fWOList, fPolyList;

	fLineParams("");
	genFList(fWOList, fPolyList);

	int num_values;
	if (is_polymer_flag)
	{
		num_values = 3 * nWells;
	}
	else
	{
		num_values = 2 * nWells;
	}

	MFVal.resize(num_values);
	MTval.resize(num_values);

	read_all_vals(MFVal, MTval, pathData, fWOList, fPolyList);
	
	for (int i = 0; i < nWells; i++)
	{
		for (int j = 0; j < ourWells[i].n_times; j++)
		{
			if (ourWells[i].theta[j] > 0)
			{
				ourWells[i].stat = 0;
			}
			else
			{
				ourWells[i].stat = 1;
			}
		}
	}

	for (int i = 0; i < nWells; i++)
	{
		ourWells[i].prev_theta = ourWells[i].theta;
		a_dim += ourWells[i].n_times;
	}

	numWell.resize(a_dim);
	numTheta.resize(a_dim);

	int dp = 0;
	for (int i = 0; i < nWells; i++)
	{
		for (int j = 0; j < ourWells[i].n_times; j++)
		{
			numWell[i + j + dp] = i;
			numTheta[i + j + dp] = j;
		}
		dp += ourWells[i].n_times - 1;
	}

	vAlphas.resize(num_values);
	ifstream in;
	in.open("alphas.txt", ios_base::in);
	if (in.is_open())
	{
		in >> alphas_flag;
		switch (alphas_flag)
		{
		case 0:
			for (int i = 0; i < num_values; i++)
				vAlphas[i] = 1;
			break;
		case 1:
			for (int i = 0; i < num_values; i++)
				in >> vAlphas[i];
			break;
		default:
			break;
		}
	}
	else
	{
		cout << "err open alphas.txt" << endl;
		system("pause");
	}
		

	
	//for (int i = 0; i < 3 * nWells; i++)
	//	vAlphas[i] = 0;
	//vAlphas[6] = 1;
	//vAlphas[7] = 1;
	//vAlphas[8] = 1;

	return 0;
}

int solvn::opti_task(string path)
{
	ofstream sum_data_f;
	vector<string> tWOFiles, tPolyFiles; // Список генерируемых файлов с данными
	vector<double> rising_func;
	int num_of_rising = 0;
	rising_func.resize(5);

	genFList(tWOFiles, tPolyFiles);
	prep_deriv_folder();

	vector<vector<double>> A(a_dim, vector<double>(a_dim));
	vector<double> b(a_dim);
	vector<double> dP(a_dim); // - diff of Parameters - Синоним для повышения читабельности кода и экономии памяти
	vector<thread> func_tread; // Массив потоков

	cout << "Launch optimization task(...)" << endl;

	double func = 1E+30, prev_func, funcWdelta, dif = 1E+30, alpha_loc;
	vector<double> vecFuncW;
	int iter_a;
	if(true)
	{
		func = functional(fWellsAll, TWellsAll, vAlphas, tWOFiles, tPolyFiles, d_iter);
		wdebugfunc(func, d_iter);
		wdebugvals(d_iter);

		cout << "func = " << func << endl;
		
		prev_func = func;

		vecFuncW.resize(a_dim);

		cout << "max_iter = " << max_iter << " func = " << func << " dif = " << dif << endl;
		for (int iter = 1; iter < max_iter && func > eps_func && dif > eps_dif; iter++)
		{
			iter_a = 0;
			alpha_loc = alpha;
			//do
			//{
				// очистить СЛАУ
			clear_SLAU(A, b);
			// собрать СЛАУ 
			// Расчет всех f(p[i] + dp[i])
			// 
			// debug func for threads
			//for (int i = 0; i < a_dim; i++)
			//	prep_data_deriv(vecFuncW, i, tWOFiles, tPolyFiles);
			
			// run threads
			//for (int i = 0; i < a_dim; i++)
			//	func_tread.push_back(thread(&solvn::prep_data_deriv, this, std::ref(vecFuncW), i, tWOFiles, tPolyFiles));
			//for (int i = 0; i < a_dim; i++)
			//	func_tread[i].join();

			int ia = a_dim % 4;
			int ib = 0;
			for (ib = 0; ib < a_dim - ia;)
			{
				for (int i = 0; i < 4; i++)
					func_tread.push_back(thread(&solvn::prep_data_deriv, this, std::ref(vecFuncW), ib + i, tWOFiles, tPolyFiles));
				for (int i = 0; i < 4; i++)
					func_tread[ib + i].join();
				ib += 4;
			}

			for(int i = 0; i < ia; i++)
				func_tread.push_back(thread(&solvn::prep_data_deriv, this, std::ref(vecFuncW), ib + i, tWOFiles, tPolyFiles));
			for (int i = 0; i < ia; i++)
				func_tread[ib + i].join();

			for (int q = 0; q < a_dim; q++)
			{
				for (int p = 0; p < a_dim; p++)
				{
					A[q][p] += deriv(func, vecFuncW[q], ourWells[numWell[q]].theta[numTheta[q]]) *
								deriv(func, vecFuncW[p], ourWells[numWell[p]].theta[numTheta[p]]);
				}

				b[q] -= deriv(func, vecFuncW[q], ourWells[numWell[q]].theta[numTheta[q]]) * (func - 0);

				for(int i = 0; i < a_dim; i++)
					b[i] += alpha_loc * (ourWells[numWell[q]].prev_theta[numTheta[q]] - 
										ourWells[numWell[q]].theta[numTheta[q]]); // Или "-"
			}

			// Добавляем регуляризацию
			for (int i = 0; i < a_dim; i++)
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
			gauss_plus(A, b, dP, a_dim);

			cout << "------------dP[q]------------" << endl;
			writeVector(dP);
			cout << "----------------------------" << endl;
			cout << endl;
			cout << endl;

			dif = Get_norm(dP);

			for(int i = 0; i < a_dim; i++)
				ourWells[numWell[i]].prev_theta[numTheta[i]] = ourWells[numWell[i]].theta[numTheta[i]];

			checkMaxMinTheta(ourWells, dP);

			write_well_param(ourWells, nWells, "./properties/");
			wdebugtheta(ourWells, d_iter + 1);

			system((pathComlex).c_str());

			if (func > prev_func)
			{
				rising_func[num_of_rising] = func;
				num_of_rising++;
			}
			else if (num_of_rising != 0)
				num_of_rising = 0;
			if (num_of_rising == 5)
				break;

			prev_func = func;

			read_all_vals(fWellsAll, TWellsAll, pathDatas, tWOFiles, tPolyFiles);
			func = functional(fWellsAll, TWellsAll, vAlphas, tWOFiles, tPolyFiles, d_iter);

			sum_data_f.open("./debug_data/sum_val.txt", std::ios::app);
			if(sum_data_f.is_open())
			{
				double data_v = 0;
				sum_data_f << "/-----------------------------------------------------------------------------------------------/\n";
				for (int nf = 0; nf < 2 * nWells; nf++)
				{
					for (int vf = 0; vf < fWellsAll[nf].size(); vf++)
						data_v += fWellsAll[nf][vf];
					sum_data_f << tWOFiles[nf] << " sum_val = " << data_v << "\n";
					data_v = 0;
				}
				if (is_polymer_flag)
				{
					int nj = 2 * nWells;
					for (int nf = 0; nf < nWells; nf++)
					{
						for (int vf = 0; vf < fWellsAll[nf + nj].size(); vf++)
							data_v += fWellsAll[nf + nj][vf];
						sum_data_f << tPolyFiles[nf] << " sum_val = " << data_v << "\n";
						data_v = 0;
					}
				}
			}
			sum_data_f.close();

			iter_a++;
			alpha_loc *= 10;

			//debug info ------------------------------>
			dthetcount = 0;
			dfunccount = 0;
			dvalscount = 0;
			d_iter++;

			func_tread.clear(); // clear vector threads

			wdebugfunc(func, d_iter);
			wdebugvals(d_iter);
		}
	}
	else
		cout << "nWells != MFvec.size()\n";
	return 0;
}
//------------------------------------------


void solvn::clear_SLAU(vector<vector<double>>& A, vector<double>& b)
{
	for (int q = 0; q < a_dim; q++)
	{
		for (int p = 0; p < a_dim; p++)
		{
			A[q][p] = 0;
		}
		b[q] = 0;
	}

}

int solvn::gauss_plus(vector<vector<double>>& A, vector<double>& b, vector<double>& x, int N)
{
	int c, err = 0, sum = 0;
	for (int i = 0; i < N; i++)
	{
		if (A[i][i] == 0)
		{
			c = 1;
			while ((i + c) < N && A[i + c][i] == 0)
				c++;
			if ((i + c) == N)
			{
				err = 1;
				break;
			}
			for (int j = i, k = 0; k < N; k++)
			{
				swap(A[j][k], A[j + c][k]);
			}
			swap(b[i], b[i + c]);
		}

		for (int j = 0; j < N; j++) {

			if (i != j) {
				float pro = A[j][i] / A[i][i];

				for (int k = 0; k < N; k++)
					A[j][k] = A[j][k] - (A[i][k]) * pro;
				b[j] = b[j] - (b[i]) * pro;
			}
		}
	}

	if (err == 1)
	{
		err = 3;
		int tj = 0;

		for (int i = 0; i < N; i++)
		{
			sum = 0;
			for (tj = 0; tj < N; tj++)
				sum = sum + A[i][tj];
			if(tj < N)
				if (sum == A[i][tj])
					err = 2;
		}
	}

	if (err == 2)
	{
		cout << "Infinite Solutions Exists" << endl;
		return 2;
	}
	else if (err == 3)
	{
		cout << "No Solution Exists" << endl;
		return 3;
	}
	else {
		for (int i = 0; i < N; i++)
			x[i] = b[i] / A[i][i];
	}

	return 0;
}

void solvn::checkGauss()
{
	vector<vector<double>> A =  {{0, 2, 1},
								 {1, 1, 2},
							     {2, 1, 1}};
	vector<double> b = {4, 6, 7};
	vector<double> x = {0, 0, 0};

	gauss_plus(A, b, x, 3);

	for (int i = 0; i < 3; i++)
		cout << "x[" << i << "] = " << x[i] << endl;
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

int solvn::checkMaxMinTheta(vector<paramWell>& fourWells, vector<double> deltaP)
{
	double max_theta;
	for (int i = 0; i < a_dim; i++)
	{
		if (!vFxThetas[i]) {
			max_theta = abs(ourWells[numWell[i]].theta[numTheta[i]] * dThetaProc);
			if (abs(deltaP[i]) < max_theta)
				ourWells[numWell[i]].theta[numTheta[i]] += deltaP[i];
			else
				ourWells[numWell[i]].theta[numTheta[i]] += max_theta * deltaP[i] / abs(deltaP[i]);
		}
	}

	for (int i = 0; i < a_dim; i++)
	{
		if (fourWells[numWell[i]].theta[numTheta[i]] > absThetaMax && fourWells[numWell[i]].stat == 0)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 1;
		}
		if (fourWells[numWell[i]].theta[numTheta[i]] < -absThetaMax && fourWells[numWell[i]].stat == 1)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 2;
		}
		if (fourWells[numWell[i]].theta[numTheta[i]] < 0 && fourWells[numWell[i]].stat == 0)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 3;
		}
		if (fourWells[numWell[i]].theta[numTheta[i]] > 0 && fourWells[numWell[i]].stat == 1)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 4;
		}
	}
	return 0;
}

void solvn::model_init()
{
	vector<string> tFileList;

	fLineParams("");
	read_well_param(ourWells, nWells, "./properties/");
	load_data(fWellsAll, TWellsAll, pathComlex, pathDatas);
	fChkFxThetas("");
}

void solvn::model_solve()
{
	string path = "";
	opti_task(path);
}

int main()
{
	solvn Solver;
	system("chcp 1251");
	//Solver.checkGauss();

	Solver.model_init();
	Solver.model_solve();
	
	cout << "---***Done!!!***---" << endl;
	return 0;
}

