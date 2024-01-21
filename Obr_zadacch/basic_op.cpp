#include "basic_op.h"

void writeMatrix(vector<vector<double>> A, vector<double> gamma_reg, int num_iter, string fileName)
{
	ofstream stream;
	stream.open(fileName, std::ios::app);
	int n = A.size();

	for (int i = 0; i < n; i++)
	{
		stream << "-";
	}
	stream << endl;
	stream << "iter = " << num_iter << "; ";

	stream << "gamma =[";
	for (int i = 0; i < n; i++)
	{
		stream << gamma_reg[i] << ";";
	}
	stream << "]" << endl;

	for (int q = 0; q < n; q++)
	{
		for (int p = 0; p < n; p++)
		{
			stream << A[q][p] << " ";
		}
		stream << endl;
	}

	for (int i = 0; i < n; i++)
	{
		stream << "-";
	}
	stream << endl;

	stream.close();
}

void writeMatrix(vector<vector<double>> A, vector<double> gamma_reg, int num_iter, int num_subiter, string fileName)
{
	ofstream stream;
	stream.open(fileName, std::ios::app);
	int n = A.size();

	for (int i = 0; i < n; i++)
	{
		stream << "-";
	}
	stream << endl;
	stream << "iter = " << num_iter << "; subiter = " << num_subiter << "; ";

	stream << "gamma =[";
	for (int i = 0; i < n; i++)
	{
		stream << gamma_reg[i] << ";";
	}
	stream << "]" << endl;

	for (int q = 0; q < n; q++)
	{
		for (int p = 0; p < n; p++)
		{
			stream << A[q][p] << " ";
		}
		stream << endl;
	}

	for (int i = 0; i < n; i++)
	{
		stream << "-";
	}
	stream << endl;

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

void normMatVec(vector<vector<double>> &A, vector<double> &b)
{
	int n = b.size();
	double min_str_val;
	vector<double> t_vec;

	t_vec.resize(n + 1);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			t_vec[j] = A[j][i];
		}
		t_vec[n] = b[i];

		min_str_val = *min_element(begin(t_vec), end(t_vec));
		if (min_str_val > 1)
		{
			for (int k = 0; k < n; k++)
			{
				A[k][i] = A[k][i] / min_str_val;
			}
			b[i] = b[i] / min_str_val;
		}
	}

}

int makeDir(string dir)
{
	error_code err;
	if (std::filesystem::exists(dir))
		return std::filesystem::is_directory(std::filesystem::status(dir));
	else
		return std::filesystem::create_directories(dir, err);
}