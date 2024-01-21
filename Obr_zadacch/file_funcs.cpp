#include "solver.h"

void solvn::prep_deriv_folder()
{
	makeDir("./debug_data/");
	for (int i = 0; i < a_dim; i++)
	{
		makeDir("./debug_data/deriv" + to_string(i) + "/");
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
		std::filesystem::copy("./time_points.txt", "./debug_data/deriv" + to_string(i));
	}

	// ============================================Для омеги:==========================================================================
	makeDir("./debug_data/deriv_omega/");
	std::filesystem::copy("./graph", "./debug_data/deriv_omega/graph", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./mesh", "./debug_data/deriv_omega/mesh", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./output_balance_Q", "./debug_data/deriv_omega/output_balance_Q", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./output3D", "./debug_data/deriv_omega/output3D", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./properties", "./debug_data/deriv_omega/properties", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./tables", "./debug_data/deriv_omega/tables", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./alignment.txt", "./debug_data/deriv_omega");
	std::filesystem::copy("./concrt140.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./DllInitAndDrawText.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./FilesConfig.cfg", "./debug_data/deriv_omega");
	std::filesystem::copy("./filtr.cfg", "./debug_data/deriv_omega");
	std::filesystem::copy("./filtr.log", "./debug_data/deriv_omega");
	std::filesystem::copy("./Filtr3D.exe", "./debug_data/deriv_omega");
	std::filesystem::copy("./libiomp5md.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./lines.txt", "./debug_data/deriv_omega");
	std::filesystem::copy("./mkl_avx.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./mkl_core.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./mkl_def.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./mkl_intel_thread.dll", "./debug_data/deriv_omega");
	std::filesystem::copy("./time_points.txt", "./debug_data/deriv_omega");
	//=============================================================================================================================================

	// ============================================Для проверки на давление в скважине:++++========================================================
	makeDir("./debug_data/check_press/");
	std::filesystem::copy("./graph", "./debug_data/check_press/graph", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./mesh", "./debug_data/check_press/mesh", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./output_balance_Q", "./debug_data/check_press/output_balance_Q", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./output3D", "./debug_data/check_press/output3D", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./properties", "./debug_data/check_press/properties", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./tables", "./debug_data/check_press/tables", std::filesystem::copy_options::recursive);
	std::filesystem::copy("./alignment.txt", "./debug_data/check_press");
	std::filesystem::copy("./concrt140.dll", "./debug_data/check_press");
	std::filesystem::copy("./DllInitAndDrawText.dll", "./debug_data/check_press");
	std::filesystem::copy("./FilesConfig.cfg", "./debug_data/check_press");
	std::filesystem::copy("./filtr.cfg", "./debug_data/check_press");
	std::filesystem::copy("./filtr.log", "./debug_data/check_press");
	std::filesystem::copy("./Filtr3D.exe", "./debug_data/check_press");
	std::filesystem::copy("./libiomp5md.dll", "./debug_data/check_press");
	std::filesystem::copy("./lines.txt", "./debug_data/check_press");
	std::filesystem::copy("./mkl_avx.dll", "./debug_data/check_press");
	std::filesystem::copy("./mkl_core.dll", "./debug_data/check_press");
	std::filesystem::copy("./mkl_def.dll", "./debug_data/check_press");
	std::filesystem::copy("./mkl_intel_thread.dll", "./debug_data/check_press");
	std::filesystem::copy("./time_points.txt", "./debug_data/check_press");
	//=============================================================================================================================================
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
		in >> is_polymer_flag;
		in >> omega_flag;

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
	if (omega_flag)
	{
		vFxThetas.resize(a_dim + 1);
	}
	else
	{
		vFxThetas.resize(a_dim);
	}

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
			if (omega_flag)
			{
				in >> vFxThetas[a_dim];
			}
		}

		in.close();
		in.clear();

		for (int i = 0; i < vFxThetas.size(); i++)
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
	ifstream in(path + "/" + tFile);

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
	vector<string> fWOList, fPolyList, fPressList;

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

	vAlphas.resize(4);
	ifstream in;
	in.open("alphas.txt", ios_base::in);
	if (in.is_open())
	{
		in >> alphas_flag;
		switch (alphas_flag)
		{
		case 0:
			for (int i = 0; i < 4; i++)
				vAlphas[i] = 1;
			break;
		case 1:
			for (int i = 0; i < 4; i++)
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

int solvn::genFList(vector<string>& fWOList, vector<string>& fPolyList, vector<string>& fPressList)
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

		//---pressure---
		filename = ourWells[i].name + "_P_small.txt";
		fPressList.push_back(filename);
		//--------------
	}

	return 0;
}