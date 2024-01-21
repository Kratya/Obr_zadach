#pragma once

#include "basic_op.h"

class paramWell
{
public:

	int num_well; // Номер скважины
	string name;
	int n_pump;
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

class paramPhase
{
public:

	int num;
	vector<string> buffstring;
	vector<double> buffnum;
	double paramPoly, paramH2O;
	double prev_paramPoly, prev_paramH2O;

	paramPhase() {};
	paramPhase(int _num, vector<string> _buffstring, vector<double> _buffnum, double _paramPoly, double _paramH2O)
	{
		this->num = _num;
		this->buffstring = _buffstring;
		this->buffnum = _buffnum;
		this->paramPoly = _paramPoly;
		this->paramH2O = _paramH2O;
	};

	~paramPhase() {};
};

class paramBord
{

public:
	int num_wells;
	vector<double> press_min, press_max;

	paramBord() {};
	paramBord(int _num_wells)
	{
		this->num_wells = _num_wells;
		this->press_min.resize(_num_wells);
		this->press_max.resize(_num_wells);
	};

	~paramBord() {};
};

int read_well_param(vector<paramWell>& ourWells, int& num_all_well, string path);
int write_well_param(vector<paramWell> ourWells, int num_all_well, string path);

int read_phase_param(paramPhase& stfPhase, string path);
int write_phase_param(paramPhase& stfPhase, string path);

int read_press_param(paramBord& ourPressParam, string path);