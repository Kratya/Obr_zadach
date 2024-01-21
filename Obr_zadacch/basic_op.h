#pragma once

#include <vector>
#include <algorithm>
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

void normMatVec(vector<vector<double>>& A, vector<double>& b);
void writeMatrix(vector<vector<double>> A, vector<double> gamma_reg, int num_iter, string fileName);
void writeMatrix(vector<vector<double>> A, vector<double> gamma_reg, int num_iter, int num_subiter, string fileName);
void writeMatrix(vector<vector<double>> A);
void writeVector(vector<double> v, string fileName);
void writeVector(vector<double> v);
vector<double> operator+(vector<double> c1, vector<double> c2);
vector<double> operator+(vector<double> c1, double aC);
vector<double> operator*(vector<double> c1, double aC);
int makeDir(string dir);