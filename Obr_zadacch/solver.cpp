#include "solver.h"

//----------------------------------------
double solvn::deriv(double fVal, double fValNew, double Old_Theta, double New_Theta)
{
	double result = 0;

	//result = (fValNew - fVal) / abs(MFThetas * delta);
	result = (fValNew - fVal) / abs(New_Theta - Old_Theta);

	return result;
}

double solvn::trapez(vector<double> tVal, vector<double> fVal)
{
	const int n = fVal.size();
	double result = 0, x1, x2;

	for (int h = 0; h < n - 1; h++)
	{
		x1 = tVal[h];
		x2 = tVal[h + 1];
		result += 0.5 * (x2 - x1) * (fVal[h] + fVal[h + 1]);
		cout << "x1:" << x1 << " x2:" << x2 << "fval[" << h << "]:" << fVal[h] << "fval[" << h + 1 << "]:" << fVal[h + 1] << endl;

	}

	cout << "result = " << result << endl;

	return result;

}

void solvn::clear_SLAU(vector<vector<double>>& A, vector<double>& b)
{
	int new_dim = a_dim;
	if (omega_flag)
	{
		new_dim += 1;
	}

	for (int q = 0; q < new_dim; q++)
	{
		for (int p = 0; p < new_dim; p++)
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
			if (tj < N)
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

vector<vector<double>> solvn::parts_func(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff)
{
	int log = 0;
	vector<double> parts_f;
	vector<vector<double>> parts_result;

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


	parts_f.resize(num_values);
	parts_result.resize(nWells);
	for (int i = 0; i < nWells; i++)
	{
		parts_result[i].resize(4);
	}

	for (int i = 0; i < num_values; i++)
	{
		parts_f[i] = trapez(TWellsAll[i], ValWellsAll[i]);
	}

	// Подбор альф
	if (alphas_flag == 2)
	{
		for (int i = 0; i < num_values; i++)
			MFCoeff[i] = 1;
		alphas_flag = 0;
	}

	for (int i = 0; i < nWells; i++)
	{
		if (parts_f[i * 2] > 0)
			parts_result[i][2] = parts_f[i * 2];
		else
			parts_result[i][0] = parts_f[i * 2];

		parts_result[i][3] = parts_f[i * 2 + 1];
	}
	if (is_polymer_flag)
	{
		for (int i = 0; i < nWells; i++)
		{
			parts_result[i][1] = parts_f[2 * nWells + i];
		}
	}
	else
	{
		for (int i = 0; i < nWells; i++)
		{
			parts_result[i][1] = 0;
		}
	}

	return parts_result;
}

vector<vector<double>> solvn::parts_func(vector<vector<double>> ValWellsAll, vector<vector<double>> TWellsAll, vector<double>& MFCoeff,
	vector<string> fWOList, vector<string> fPolyList, int iter)
{
	int log = 0;
	vector<double> parts_f;
	vector<vector<double>> parts_result;

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

	parts_f.resize(num_values);
	parts_result.resize(nWells);
	for (int i = 0; i < nWells; i++)
	{
		parts_result[i].resize(4);
	}

	for (int i = 0; i < num_values; i++)
	{
		parts_f[i] = trapez(TWellsAll[i], ValWellsAll[i]);
	}

	for (int i = 0; i < nWells; i++)
	{
		if (parts_f[i * 2] > 0)
			parts_result[i][2] = parts_f[i * 2];
		else
			parts_result[i][0] = parts_f[i * 2];

		parts_result[i][3] = parts_f[i * 2 + 1];
	}
	if (is_polymer_flag)
	{
		for (int i = 0; i < nWells; i++)
		{
			parts_result[i][1] = parts_f[2 * nWells + i];
		}
	}
	else
	{
		for (int i = 0; i < nWells; i++)
		{
			parts_result[i][1] = 0;
		}
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
		integr_res += pow(MFCoeff[0] * parts_result[i][0], 2);
	}

	out.open("./debug_data/integr_s1_out.txt", std::ios::app);
	out << iter << ". " << integr_res << endl;
	out.close();
	integr_res = 0;

	for (int i = 0; i < nWells; i++)
	{
		integr_res += pow(MFCoeff[2] * parts_result[i][2], 2);
	}
	out.open("./debug_data/integr_s1_ext.txt", std::ios::app);
	out << iter << ". " << integr_res << endl;
	out.close();
	integr_res = 0;

	for (int i = 0; i < nWells; i++)
	{
		integr_res += pow(MFCoeff[1] * parts_result[i][1], 2);
	}
	out.open("./debug_data/integr_s2.txt", std::ios::app);
	out << iter << ". " << integr_res << endl;
	out.close();
	integr_res = 0;

	if (is_polymer_flag) {
		int j = 2 * nWells;
		for (int i = 0; i < nWells; i++)
		{
			integr_res += pow(MFCoeff[3] * parts_result[i][3], 2);
		}
		out.open("./debug_data/integr_c3.txt", std::ios::app);
		out << iter << ". " << integr_res << endl;
		out.close();
		integr_res = 0;
	}

	if (omega_flag)
	{
		out.open("./debug_data/omega.txt", std::ios::app);
		out << iter << ". " << "P = " << std::fixed << std::setprecision(8) << ourOmega.paramPoly << ", H2O = " << std::fixed << std::setprecision(8) << ourOmega.paramH2O << endl;
		out.close();
	}

	return parts_result;
}

double solvn::functional(vector<vector<double>> parts_result, vector<double> MFCoeff)
{
	double func = 0;

	for (int i = 0; i < nWells; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			func += pow(MFCoeff[j] * parts_result[i][j], 2);
		}
	}

	return func;
}

void solvn::prep_data_deriv(vector<vector<vector<double>>>& vecFunc, vector<double>& new_thetas, int i, vector<string> fWOFiles, vector<string> fPolyFiles)
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

	tTheta = plDeltaWells[numWell[i]].theta[numTheta[i]];
	plDeltaWells[numWell[i]].theta[numTheta[i]] = plDeltaWells[numWell[i]].theta[numTheta[i]] +
		abs(plDeltaWells[numWell[i]].theta[numTheta[i]] * delta);
	cout << "write1 deriv well params" << endl;
	write_well_param(plDeltaWells, nWells, "./debug_data/deriv" + to_string(i) + "/properties/");
	cout << "launch deriv filtr3d" << endl;
	system(("c: && cd " + cwd.string() + "\\debug_data\\deriv" + to_string(i) + " \\ && " + pathComlex).c_str());

	cout << "read_all_vals deriv" << endl;
	read_all_vals(fWellsAllNew, TWellsAllNew, "./debug_data/deriv" + to_string(i) + "/output", fWOFiles, fPolyFiles);
	vecFunc[i] = parts_func(fWellsAllNew, TWellsAllNew, vAlphas);
	new_thetas[i] = plDeltaWells[numWell[i]].theta[numTheta[i]];

	plDeltaWells[numWell[i]].theta[numTheta[i]] = tTheta;
	cout << "write2 deriv well params" << endl;
	write_well_param(plDeltaWells, nWells, "./debug_data/deriv" + to_string(i) + "/properties/");
}

void solvn::prep_data_deriv_omega(vector<vector<vector<double>>>& vecFunc, vector<double>& new_thetas, vector<string> fWOFiles, vector<string> fPolyFiles)
{
	vector<vector<double>> fWellsAllNew, TWellsAllNew;
	paramPhase tOmega;

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

	tOmega = ourOmega;

	tOmega.paramPoly += tOmega.paramPoly * delta;
	tOmega.paramH2O = 1 - tOmega.paramPoly;

	cout << "write1 deriv omega params" << endl;
	write_phase_param(tOmega, "./debug_data/deriv_omega/properties/");
	cout << "launch deriv filtr3d" << endl;
	system(("c: && cd " + cwd.string() + "\\debug_data\\deriv_omega" + " \\ && " + pathComlex).c_str());

	cout << "read_all_vals deriv" << endl;
	read_all_vals(fWellsAllNew, TWellsAllNew, "./debug_data/deriv_omega/output", fWOFiles, fPolyFiles);
	vecFunc[a_dim] = parts_func(fWellsAllNew, TWellsAllNew, vAlphas);
	new_thetas[a_dim] = tOmega.paramPoly;

	cout << "write2 deriv well params" << endl;
	write_phase_param(ourOmega, "./debug_data/deriv_omega/properties/");

	// debug this
	//ofstream dout;
	//dout.open("./debug_data/deriv_omega/debug_deriv_omega.txt", std::ios::app);
	//dout << ourOmega.paramPoly << "; " << tOmega.paramPoly << endl;
	//dout.close();
	// end debug
}



int solvn::opti_task(string path)
{
	ofstream sum_data_f;
	vector<string> tWOFiles, tPolyFiles; // Список генерируемых файлов с данными
	vector<double> rising_func;
	int num_of_rising = 0;
	rising_func.resize(5);

	genFList(tWOFiles, tPolyFiles, ourPressFiles);
	prep_deriv_folder();

	vector<thread> func_tread; // Массив потоков

	vector<vector<vector<double>>> vecFuncW;
	vector<vector<double>> vec_func;
	vector<double> gamma_reg;
	vector<double> new_thetas;

	double func = 1E+30, prev_func, funcWdelta, dif = 1E+30;
	int res_check_th = 10, res_check_P = 10, iters_of_check = 0;

	int new_dim = 0;
	if (omega_flag)
	{
		new_dim = a_dim + 1;
	}
	else
	{
		new_dim = a_dim;
	}

	vector<vector<double>> A(new_dim, vector<double>(new_dim));
	vector<double> b(new_dim);
	vector<double> dP(new_dim); // - diff of Parameters

	gamma_reg.resize(new_dim);
	new_thetas.resize(new_dim);

	for (int i = 0; i < a_dim; i++)
	{
		gamma_reg[i] = alpha;
	}
	if (omega_flag)
		gamma_reg[a_dim] = alpha_gamma;

	int iter_a;
	if (true)
	{
		vec_func = parts_func(fWellsAll, TWellsAll, vAlphas, tWOFiles, tPolyFiles, d_iter);
		func = functional(vec_func, vAlphas);

		wdebugfunc(func, d_iter);
		wdebugvals(gamma_reg, d_iter);
		wdebugsumvals(tPolyFiles, d_iter);

		cout << "func = " << func << endl;

		prev_func = func;

		vecFuncW.resize(new_dim);
		vector<int> error_gamma;
		error_gamma.resize(new_dim);

		cout << "max_iter = " << max_iter << " func = " << func << " dif = " << dif << endl;
		for (int iter = 1; iter < max_iter && func > eps_func && dif > eps_dif; iter++)
		{
			iter_a = 0;
			//alpha_loc = alpha;


			// очистить СЛАУ
			clear_SLAU(A, b);
			// собрать СЛАУ 
			// Расчет всех f(p[i] + dp[i])

			// run threads
			int ia = a_dim % 4;
			int ib = 0;
			for (ib = 0; ib < a_dim - ia;)
			{
				for (int i = 0; i < 4; i++)
					func_tread.push_back(thread(&solvn::prep_data_deriv, this, std::ref(vecFuncW), std::ref(new_thetas), ib + i, tWOFiles, tPolyFiles));
				for (int i = 0; i < 4; i++)
					func_tread[ib + i].join();
				ib += 4;
			}

			for (int i = 0; i < ia; i++)
				func_tread.push_back(thread(&solvn::prep_data_deriv, this, std::ref(vecFuncW), std::ref(new_thetas), ib + i, tWOFiles, tPolyFiles));
			for (int i = 0; i < ia; i++)
				func_tread[ib + i].join();

			//___________________Можно инт. в потоки---------------------------:
			if (omega_flag)
			{
				prep_data_deriv_omega(vecFuncW, new_thetas, tWOFiles, tPolyFiles);
			}
			//-----------------------------------------------------------------
			vector<int> fix_th_iter = vFxThetas;

			for (int i = 0; i < new_dim; i++)
			{
				error_gamma[i] = 0;
			}

			while ((res_check_th != 0 || res_check_P != 0) && iters_of_check < 200)
			{
				iters_of_check++;

				for (int q = 0; q < a_dim; q++)
				{
					for (int p = 0; p < a_dim; p++)
					{
						for (int k = 0; k < 4; k++)
						{
							for (int w = 0; w < nWells; w++)
							{
								A[q][p] += vAlphas[k] * deriv(vec_func[w][k], vecFuncW[q][w][k], ourWells[numWell[q]].theta[numTheta[q]], new_thetas[q]) *
									deriv(vec_func[w][k], vecFuncW[q][w][k], ourWells[numWell[p]].theta[numTheta[p]], new_thetas[p]);
							}
						}
					}
					for (int k = 0; k < 4; k++)
					{
						for (int w = 0; w < nWells; w++)
						{
							b[q] -= vAlphas[k] * deriv(vec_func[w][k], vecFuncW[q][w][k], ourWells[numWell[q]].theta[numTheta[q]], new_thetas[q]) * (vec_func[w][k] - 0);
						}
					}
				}

				if (omega_flag)
				{
					for (int q = 0; q < a_dim; q++)
					{
						for (int k = 0; k < 4; k++)
						{
							for (int w = 0; w < nWells; w++)
							{
								A[q][a_dim] += vAlphas[k] * deriv(vec_func[w][k], vecFuncW[q][w][k], ourWells[numWell[q]].theta[numTheta[q]], new_thetas[q]) *
									deriv(vec_func[w][k], vecFuncW[q][w][k], ourOmega.paramPoly, new_thetas[a_dim]);
							}
						}
						A[a_dim][q] = A[q][a_dim];
					}

					for (int k = 0; k < 4; k++)
					{
						for (int w = 0; w < nWells; w++)
						{
							b[a_dim] -= vAlphas[k] * deriv(vec_func[w][k], vecFuncW[a_dim][w][k], ourOmega.paramPoly, new_thetas[a_dim]) * (vec_func[w][k] - 0);

						}
					}

					for (int k = 0; k < 4; k++)
					{
						for (int w = 0; w < nWells; w++)
						{
							A[a_dim][a_dim] += vAlphas[k] * deriv(vec_func[w][k], vecFuncW[a_dim][w][k], ourOmega.paramPoly, new_thetas[a_dim]) *
								deriv(vec_func[w][k], vecFuncW[a_dim][w][k], ourOmega.paramPoly, new_thetas[a_dim]);
						}
					}
				}

				normMatVec(A, b);

				// Добавляем регуляризацию
				for (int i = 0; i < new_dim; i++)
					A[i][i] += gamma_reg[i] * A[i][i];

				cout << "----------A[q][p]-----------" << endl;
				writeMatrix(A);
				//writeMatrix(A, gamma_reg, d_iter, "./debug_data/deriv_omega/matrix_A.txt");
				writeMatrix(A, gamma_reg, d_iter, iters_of_check, "./debug_data/deriv_omega/matrix_A.txt");
				cout << "----------------------------" << endl;
				cout << endl;

				cout << "------------b[q]------------" << endl;
				writeVector(b);
				cout << "----------------------------" << endl;
				cout << endl;

				// решить СЛАУ
				gauss_plus(A, b, dP, new_dim);

				cout << "------------dP[q]------------" << endl;
				writeVector(dP);
				cout << "----------------------------" << endl;
				cout << endl;
				cout << endl;

				dif = Get_norm(dP);

				for (int i = 0; i < a_dim; i++)
					ourWells[numWell[i]].prev_theta[numTheta[i]] = ourWells[numWell[i]].theta[numTheta[i]];
				if (omega_flag)
				{
					ourOmega.prev_paramPoly = ourOmega.paramPoly;
					ourOmega.prev_paramH2O = ourOmega.paramH2O;
				}

				res_check_th = checkMaxMinTheta(ourWells, dP, gamma_reg, error_gamma, fix_th_iter, d_iter, iters_of_check);
				res_check_P  = checkMaxMinPress(ourPressFiles, ourPress, gamma_reg, error_gamma, nWells, fix_th_iter, d_iter, iters_of_check);
				wdebugalphas(gamma_reg, fix_th_iter, error_gamma, d_iter, iters_of_check);
				// очистить СЛАУ
				clear_SLAU(A, b);
			}
			if (iters_of_check >= 200)
			{
				cout << "!!! --- iters_of_check >= 100 --- !!!" << endl;
				ofstream err_while;
				err_while.open("./debug_data/debug_iters.txt", std::ios::app);

				err_while << "iter = " << d_iter << "; sub_iter = " << iters_of_check << "; res_check_th = " << res_check_th << "; res_check_P = " << res_check_P << endl;

				err_while.close();
			}
				
			//cout << "!!! --- iters_of_check >= 20 --- !!!" << endl;
			else
				cout << "iters_of_check = " << iters_of_check << endl;

			iters_of_check = 0;
			res_check_th = 10;
			res_check_P = 10;

			write_well_param(ourWells, nWells, "./properties/");
			write_phase_param(ourOmega, "./properties/");

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
			vec_func = parts_func(fWellsAll, TWellsAll, vAlphas, tWOFiles, tPolyFiles, d_iter);
			func = functional(vec_func, vAlphas);

			sum_data_f.open("./debug_data/sum_val.txt", std::ios::app);
			if (sum_data_f.is_open())
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
			//alpha_loc *= 10;

			//debug info ------------------------------>
			dthetcount = 0;
			dfunccount = 0;
			dvalscount = 0;
			d_iter++;

			func_tread.clear(); // clear vector threads

			wdebugfunc(func, d_iter);
			wdebugvals(gamma_reg, d_iter);
			wdebugsumvals(tPolyFiles, d_iter);
			wdedugoutprop(d_iter);
			wdebugsummacc(d_iter);

			for (int i = 0; i < a_dim; i++)
			{
				gamma_reg[i] = alpha;
			}
			if (omega_flag)
				gamma_reg[a_dim] = alpha_gamma;
		}
	}
	else
		cout << "nWells != MFvec.size()\n";
	return 0;
}

void solvn::model_init()
{
	vector<string> tFileList;

	fLineParams("");
	read_well_param(ourWells, nWells, "./properties/");
	read_phase_param(ourOmega, "./properties/");
	read_press_param(ourPress, "./properties/");

	load_data(fWellsAll, TWellsAll, pathComlex, pathDatas);
	fChkFxThetas("");
}

void solvn::model_solve()
{
	string path = "";
	opti_task(path);
}