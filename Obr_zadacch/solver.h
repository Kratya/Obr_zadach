#pragma once

#include "params.h"

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

	double alpha = 1E-7;
	double alpha_gamma = 1E-07;

	//double delta = 0.005;
	double delta = 0.05;
	double dThetaProc = 0.3;
	double absThetaMax = 100;
	double vStar = -100;

	int is_polymer_flag = 0;


	// Переменные для отладки
	int dthetcount = 0;
	int dfunccount = 0;
	int dvalscount = 0;
	int d_iter = 0;

	int alphas_flag;
	int omega_flag = 0;

	vector<double>	vAlphas;
	vector<int> vFxThetas;

	int nWells; // количество скважин
	int a_dim; // размерность, зависящая от кол-ва мощностей
	vector<int> numWell;
	vector<int> numTheta;

	vector<paramWell> ourWells;
	paramPhase ourOmega;
	paramBord ourPress;
	vector<string> ourPressFiles;

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
	int genFList(vector<string>& fWOList, vector<string>& fPolyList, vector<string>& fPressList);
	int read_all_vals(vector<vector<double>>& AllFVals, vector<vector<double>>& AllTVals, string path, vector<string> woFiles, vector<string> polyFiles);
	//------------------------------------------------------------------------------------------------------------------
	int opti_task(string path);
	// Расчет производных, интеграл и функционал
	double trapez(vector<double> tVal, vector<double> fVal);
	double deriv(double fVal, double fValNew, double Old_Theta, double New_Theta);

	int wdebugfunc(double dfunc, int d_iter);
	int wdebugtheta(vector<paramWell> fourWells, int d_iter);
	int wdebugvals(vector<double>gamma_reg, int d_iter);
	int wdebugalphas(vector<double> gamma_reg, vector<int> error_gamma, vector<int> fix_th_iter, int d_iter, int sub_iter);
	int wdebugsumvals(vector<string> poly_files, int d_iter);
	int wdedugoutprop(int d_iter);
	int wdebugsummacc(int d_iter);
	int checkMaxMinPress(vector<string> fPressList, paramBord tPress, vector<double>& gamma_reg, vector<int> &error_gamma, int num_wells, vector<int> &fix_th_iter, int d_iter, int sub_iter);
	int checkMaxMinTheta(vector<paramWell>& fourWells, vector<double> deltaP, vector<double>& gamma_reg, vector<int> &error_gamma, vector<int> &fix_th_iter, int d_iter, int sub_iter);
	void prep_deriv_folder();

	vector<vector<double>> parts_func(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff);
	vector<vector<double>> parts_func(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff,
		vector<string> fWOList, vector<string> fPolyList, int iter);
	void prep_data_deriv(vector<vector<vector<double>>>& vecFunc, vector<double>& new_thetas, int i, vector<string> fWOFiles, vector<string> fPolyFiles);
	void prep_data_deriv_omega(vector<vector<vector<double>>>& vecFunc, vector<double>& new_thetas, vector<string> fWOFiles, vector<string> fPolyFiles);
	double functional(vector<vector<double>> parts_result, vector<double> MFCoeff);
};

int test_phase_rw(paramPhase& stfPhase);