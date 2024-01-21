#include "solver.h"
#define GAMMA_VAL 6
#define NUM_CHECKS_SUB_ITERS 50


int solvn::checkMaxMinPress(vector<string> fPressList, paramBord tPress, vector<double>& gamma_reg, vector<int> &error_gamma, int num_wells, vector<int> &fix_th_iter, int d_iter, int sub_iter)
{
	ifstream in;
	ofstream out;
	std::filesystem::path cwd = std::filesystem::current_path();

	double wpress_min, wpress_max, prev_press, cur_press, read_press;
	int days;
	int error_code = 0;

	for (int indx_well = 0; indx_well < num_wells; indx_well++)
	{

		wpress_min = tPress.press_min[indx_well];
		wpress_max = tPress.press_max[indx_well];

		in.open("./output/" + fPressList[indx_well]);
		out.open("./debug_data/debug_check_press.txt", std::ios::app);

		in >> days >> cur_press;
		while (in >> days >> read_press)
		{
			prev_press = cur_press;
			cur_press = read_press;
			if (prev_press >= wpress_max)
			{
				out << "iter:" << d_iter << "; sub_iter = " << sub_iter <<  "; prev_press > wpress_max" << "; num_well = " << indx_well << "; prev_press = " << prev_press << endl;
				for (int i = 0; i < a_dim; i++)
				{
					if ((numWell[i] == indx_well) && !fix_th_iter[i])
					{

						gamma_reg[i] *= GAMMA_VAL;
						error_gamma[i] += 1;
						error_code = 1;
					}
				}
			}
			if (prev_press <= wpress_min)
			{
				out << "iter:" << d_iter << "; sub_iter = " << sub_iter << "; prev_press < wpress_min" << "; num_well = " << indx_well << "; prev_press = " << prev_press << endl;
				for (int i = 0; i < a_dim; i++)
				{
					if ((numWell[i] == indx_well) && !fix_th_iter[i])
					{
						gamma_reg[i] *= GAMMA_VAL;
						error_gamma[i] += 1;
						error_code = 1;
					}
				}
			}
			if (error_code)
				break;
		}

		if (error_code != 0)
		{
			for (int i = 0; i < a_dim; i++)
			{
				if (!fix_th_iter[i])
				{
					ourWells[numWell[i]].theta[numTheta[i]] = ourWells[numWell[i]].prev_theta[numTheta[i]];
					if (error_gamma[i] >= NUM_CHECKS_SUB_ITERS)
					{
						fix_th_iter[i] = 1;
					}
				}
			}
			/*
			if (omega_flag)
			{
				if (!fix_th_iter[a_dim]) 
				{
					ourOmega.paramPoly = ourOmega.prev_paramPoly;
					ourOmega.paramH2O = ourOmega.prev_paramH2O;
					if (NUM_CHECKS_SUB_ITERS)
					{
						fix_th_iter[a_dim] = 1;
					}
				}
			}
			*/
		}

		out.close();
		in.close();
	}

	return error_code;
}

int solvn::checkMaxMinTheta(vector<paramWell>& fourWells, vector<double> deltaP, vector<double>& gamma_reg, vector<int> &error_gamma, vector<int>& fix_th_iter, int d_iter, int sub_iter)
{
	double max_theta, min_theta;
	double error_code = 0;

	ofstream out;
	out.open("./debug_data/debug_check_dP.txt", std::ios::app);

	if (out.is_open())
	{
		out << "iter = " << d_iter << "; [";
		for (int i = 0; i < a_dim; i++)
		{
			out << deltaP[i] << ", ";
		}
		if (omega_flag)
		{
			out << deltaP[a_dim];
		}
		out << "]" << endl;
		out.close();
	}
	else
	{
		cout << "can't open debug_check_dP.txt" << endl;
		//terminate();
	}

	for (int i = 0; i < a_dim; i++)
	{
		if (!fix_th_iter[i])
		{
			max_theta = abs(ourWells[numWell[i]].theta[numTheta[i]] * dThetaProc);
			if (abs(deltaP[i]) < max_theta)
				ourWells[numWell[i]].theta[numTheta[i]] += deltaP[i];
			else
			{
				gamma_reg[i] *= GAMMA_VAL;
				error_gamma[i] += 1;
				error_code = 1;
				//ourWells[numWell[i]].theta[numTheta[i]] += max_theta * deltaP[i] / abs(deltaP[i]);
			}
		}
	}

	if (omega_flag)
	{
		ourOmega.prev_paramPoly = ourOmega.paramPoly;
		ourOmega.prev_paramH2O = ourOmega.paramH2O;

		if (!fix_th_iter[a_dim])
		{
			max_theta = ourOmega.paramPoly + 0.10;
			min_theta = ourOmega.paramPoly - 0.10;

			if ((deltaP[a_dim] + ourOmega.paramPoly) < max_theta && (deltaP[a_dim] + ourOmega.paramPoly) > min_theta
				&& (deltaP[a_dim] + ourOmega.paramPoly) < 1 && (deltaP[a_dim] + ourOmega.paramPoly) > 0)
			{
				ourOmega.paramPoly += deltaP[a_dim];
				ourOmega.paramH2O = 1 - ourOmega.paramPoly;
			}
			else
			{
				gamma_reg[a_dim] *= GAMMA_VAL;
				error_gamma[a_dim] += 1;
				error_code = 1;
				//ourWells[numWell[i]].theta[numTheta[i]] += max_theta * deltaP[i] / abs(deltaP[i]);
			}
		}
	}

	if (error_code == 1)
	{
		for (int i = 0; i < a_dim; i++)
		{
			if (!fix_th_iter[i])
			{
				ourWells[numWell[i]].theta[numTheta[i]] = ourWells[numWell[i]].prev_theta[numTheta[i]];
				if (error_gamma[i] >= NUM_CHECKS_SUB_ITERS)
				{
					fix_th_iter[i] = 1;
				}
			}
		}

		if (omega_flag)
		{
			if (!fix_th_iter[a_dim]) 
			{
				ourOmega.paramPoly = ourOmega.prev_paramPoly;
				ourOmega.paramH2O = ourOmega.prev_paramH2O;

				if (error_gamma[a_dim] >= NUM_CHECKS_SUB_ITERS)
				{
					fix_th_iter[a_dim] = 1;
				}
			}
		}

		return error_code;
	}

	// Проверка на выход из диапазона (желательно доработать)
	for (int i = 0; i < a_dim; i++)
	{
		if (fourWells[numWell[i]].theta[numTheta[i]] > absThetaMax && fourWells[numWell[i]].stat == 0)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 2;
		}
		if (fourWells[numWell[i]].theta[numTheta[i]] < -absThetaMax && fourWells[numWell[i]].stat == 1)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 3;
		}
		if (fourWells[numWell[i]].theta[numTheta[i]] < 0 && fourWells[numWell[i]].stat == 0)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 4;
		}
		if (fourWells[numWell[i]].theta[numTheta[i]] > 0 && fourWells[numWell[i]].stat == 1)
		{
			fourWells[numWell[i]].theta[numTheta[i]] = fourWells[numWell[i]].prev_theta[numTheta[i]];
			return 5;
		}
	}

	return error_code;
}