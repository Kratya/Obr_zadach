#include "params.h"

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

int read_phase_param(paramPhase& stfPhase, string path)
{
	ifstream in(path + "/phasecomplist.txt");
	vector<string> buffstring(10);
	vector<double> buffnum(10);
	double paramPoly, paramH2O;
	int num;

	in >> num;

	getline(in, buffstring[0], '\t');
	getline(in, buffstring[0], '\t');
	getline(in, buffstring[1], '\t');
	getline(in, buffstring[2], '\n');

	getline(in, buffstring[3], '\t');
	in >> buffnum[0];
	in >> buffnum[1];
	in >> buffnum[2];

	getline(in, buffstring[4], '\n');
	getline(in, buffstring[4], '\t');
	in >> buffnum[3];
	in >> buffnum[4];
	in >> buffnum[5];

	getline(in, buffstring[5], '\t');
	getline(in, buffstring[5], '\t');
	getline(in, buffstring[6], '\t');
	getline(in, buffstring[7], '\n');

	getline(in, buffstring[8], '\t');
	in >> paramH2O;
	in >> buffnum[6];
	in >> paramPoly;

	getline(in, buffstring[9], '\n');
	getline(in, buffstring[9], '\t');
	in >> buffnum[7];
	in >> buffnum[8];
	in >> buffnum[9];

	paramPhase new_param_phase(num, buffstring, buffnum, paramPoly, paramH2O);
	stfPhase = new_param_phase;

	return 0;
}

int write_phase_param(paramPhase& stfPhase, string path)
{
	ofstream out(path + "/phasecomplist.txt");
	if (out.is_open())
	{
		out << stfPhase.num << '\n';
		out << '\t' << stfPhase.buffstring[0] << '\t' << stfPhase.buffstring[1] << '\t' << stfPhase.buffstring[2] << '\n';
		out << stfPhase.buffstring[3] << '\t' << stfPhase.buffnum[0] << '\t' << stfPhase.buffnum[1] << '\t' << stfPhase.buffnum[2] << '\n';
		out << stfPhase.buffstring[4] << '\t' << stfPhase.buffnum[3] << '\t' << stfPhase.buffnum[4] << '\t' << stfPhase.buffnum[5] << '\n';
		out << '\t' << stfPhase.buffstring[5] << '\t' << stfPhase.buffstring[6] << '\t' << stfPhase.buffstring[7] << '\n';
		out << stfPhase.buffstring[8] << '\t' << std::fixed << std::setprecision(8) << stfPhase.paramH2O << '\t' << stfPhase.buffnum[6] << '\t' << std::fixed << std::setprecision(8) << stfPhase.paramPoly << '\n';
		out << stfPhase.buffstring[9] << '\t' << stfPhase.buffnum[7] << '\t' << stfPhase.buffnum[8] << '\t' << stfPhase.buffnum[9] << '\n';
	}
	else
		return -1;

	return 0;
}

int read_press_param(paramBord &ourPressParam, string path)
{
	ifstream in(path + "/wellsPressureRestrictions.txt");
	string buffstring;
	int num_wells = 0;

	in >> num_wells;

	ourPressParam.num_wells = num_wells;
	ourPressParam.press_min.resize(num_wells);
	ourPressParam.press_max.resize(num_wells);
	
	for (int i = 0; i < num_wells; i++)
	{
		in >> buffstring >> ourPressParam.press_min[i] >> ourPressParam.press_max[i];
	}

	return 0;
}