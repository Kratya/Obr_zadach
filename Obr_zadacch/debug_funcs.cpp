#include "solver.h"

void solvn::checkGauss()
{
	vector<vector<double>> A = { {0, 2, 1},
								 {1, 1, 2},
								  {2, 1, 1} };
	vector<double> b = { 4, 6, 7 };
	vector<double> x = { 0, 0, 0 };

	gauss_plus(A, b, x, 3);

	for (int i = 0; i < 3; i++)
		cout << "x[" << i << "] = " << x[i] << endl;
}

int solvn::wdebugsummacc(int d_iter)
{
	ifstream in;
	ofstream out;

	int days_out = 0;
	double summ_ext = 0.0;

	in.open("./output/SumAccumulatedOilWells.txt");
	out.open("./debug_data/summ_acc_out.txt", std::ios::app);

	while (in >> days_out)
	{
		in >> summ_ext;
	}

	out << d_iter << ": " << summ_ext << endl;

	in.close();
	out.close();

	return 0;
}

int solvn::wdedugoutprop(int d_iter)
{
	makeDir("./debug_outprop/");

	std::filesystem::copy("./output", "./debug_outprop/debug_output_iter_" + std::to_string(d_iter), std::filesystem::copy_options::recursive);
	std::filesystem::copy("./properties", "./debug_outprop/debug_properties_iter_" + std::to_string(d_iter), std::filesystem::copy_options::recursive);

	return 0;
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

int solvn::wdebugalphas(vector<double> gamma_reg, vector<int> fix_th_iter, vector<int> error_gamma, int d_iter, int sub_iter)
{
	makeDir("./debug_data/");
	ofstream deb_alphas; //f_vals;

	deb_alphas.open("./debug_data/debug_alpas.txt", std::ios::app);
	deb_alphas << "iter = " << d_iter << "; sub_iter = " << sub_iter << "; [";
	for (int i = 0; i < gamma_reg.size(); i++)
	{
		deb_alphas << gamma_reg[i] << ";";
	}
	deb_alphas << "]; fix = [";

	for (int i = 0; i < fix_th_iter.size(); i++)
	{
		deb_alphas << fix_th_iter[i] << ";";
	}
	deb_alphas << "]; err_g = [";

	for (int i = 0; i < error_gamma.size(); i++)
	{
		deb_alphas << error_gamma[i] << ";";
	}
	deb_alphas << "];" << endl;

	deb_alphas.close();

	return 0;
}

int solvn::wdebugvals(vector<double> gamma_reg, int d_iter)
{
	makeDir("./debug_data/");
	ofstream deb_vals; //f_vals;
	ifstream f_vals;

	deb_vals.open("./debug_data/debug_vals.txt", std::ios::app);
	f_vals.open("./output/Well1_s2.txt");

	deb_vals << d_iter << ": ";
	double TVal, FVal;

	while (f_vals >> TVal >> FVal);
	deb_vals << FVal << " | ";
	deb_vals << "\n";

	f_vals.close();
	deb_vals.close();

	deb_vals.open("./debug_data/debug_alpas.txt", std::ios::app);
	deb_vals << "iter =" << d_iter << ":";
	for (int i = 0; i < gamma_reg.size(); i++)
	{
		deb_vals << gamma_reg[i] << ";";
	}
	deb_vals << endl;
	deb_vals.close();

	return 0;
}

int solvn::wdebugsumvals(vector<string> poly_files, int d_iter)
{
	makeDir("./debug_data/");

	string new_poly_files_name;
	ofstream deb_vals; //f_vals;
	ifstream inp_vals;

	double t_val, f_val, sum_val = 0.0;

	for (int i = 0; i < poly_files.size(); i++)
	{
		new_poly_files_name = "sum_" + poly_files[i];
		inp_vals.open("./output/" + poly_files[i]);

		sum_val = 0.0;
		while(inp_vals >> t_val >> f_val)
		{
			sum_val += f_val * 30;
			if(t_val >= 600)
				break;
		}

		deb_vals.open("./debug_data/" + new_poly_files_name, std::ios::app);
		deb_vals << d_iter << "   " << sum_val << endl;
		deb_vals.close();
		inp_vals.close();
	} 

	return 0;
}

int test_phase_rw(paramPhase& stfPhase)
{
	int result = 0;

	if (!read_phase_param(stfPhase, "./properties/"))
		result = write_phase_param(stfPhase, "./");
	else return -1;

	return result;
}