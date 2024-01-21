#define _USE_MATH_DEFINES

#include "solver.h"

/// <summary>
/// Вывод матрицы в файл
/// </summary>
/// <param name="A"> - Матрица для вывода</param>
/// <param name="fileName"> - Имя файла, куда выводить (Относительное или абсолютное)</param>

int main()
{
	//paramPhase stfPhase;
	solvn Solver;
	system("chcp 1251");
	Solver.checkGauss();

	Solver.model_init();
	Solver.model_solve();
	
	//cout << "---***Done!!!***---" << endl;

	//int test_res = 0;
	//test_res = test_phase_rw(stfPhase);
	//cout << "test_res = " << test_res << endl;

	return 0;
}



