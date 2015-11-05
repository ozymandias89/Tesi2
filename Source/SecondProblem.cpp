/*
 * SecondProblem.cpp
 *
 *  Created on: 04 nov 2015
 *      Author: riccardo
 */

#include "SecondProblem.h"

SecondProblem::SecondProblem() {
	// TODO Auto-generated constructor stub
}

SecondProblem::~SecondProblem() {
	// TODO Auto-generated destructor stub
}


void SecondProblem::setupSP(CEnv env, Prob lp, int num_rows) {

	{
		// variables
		static const char* varType = NULL;
		double obj = 1.0;
		double lb = -CPX_INFBOUND;
		double ub = 0.0;

		// variable u_0
		snprintf(name, NAME_SIZE, "u_%i", 0);
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		ub = CPX_INFBOUND;

		// variable u
		for (int i = 1; i <= num_rows; i++) {
			snprintf(name, NAME_SIZE, "u_%i", i);
			varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
					&varName);
		}

		// variables a
		for (int i = 0; i < N; i++) {
			snprintf(name, NAME_SIZE, "a_%i", i);
			varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
					&varName);
		}

		// variables b
		snprintf(name, NAME_SIZE, "b");
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		// variables v_0
		lb = 0.0;
		snprintf(name, NAME_SIZE, "v_%i", 0);
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		lb = -CPX_INFBOUND;

		// variables v
		for (int i = 1; i <= num_rows; i++) {
			snprintf(name, NAME_SIZE, "v_%i", i);
			varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
					&varName);
		}

	}

	// constraints A_T * u + e_k * u_0

	{
		std::vector<int> idx;
		std::vector<double> coef;

		// --------------------------------------------------
		//  A_T * u
		// --------------------------------------------------
		for (int i = 0; i < N; i++) {
			char sense = 'L';
			int matbeg = 0;
			double rhs = 0;
			int nzcnt = 0;

			int iter = 0;
			int u = 1;

			while (iter < num_constraint) {

				if (A[iter][i] != 0) {
					idx.push_back(u);
					coef.push_back(A[iter][i]);
					nzcnt++;
				}
				iter++;
				u++;
			}


			// --------------------------------------------------
			//  e_k * u0
			// --------------------------------------------------
			if (i == k) {
				idx.push_back(0);
				coef.push_back(1);
				nzcnt++;
			}

			// --------------------------------------------------
			//  -a_i
			// --------------------------------------------------
			idx.push_back(num_rows + 1 + i);
			coef.push_back(-1);
			nzcnt++;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

	int v_0;

	//constraint b_T * u + u_0 * gamma

	{
		std::vector<int> idx;
		std::vector<double> coef;

		char sense = 'G';
		int matbeg = 0;
		double rhs = 0;
		int nzcnt = 0;
		int u = 1;

		for (int i = 0; i < num_constraint; i++) {

			if (b[i] != 0) {
				idx.push_back(u);
				coef.push_back(b[i]);
				nzcnt++;
			}
			u++;

		}


		if (gam != 0) {
			idx.push_back(0);
			coef.push_back(gam);
			nzcnt++;
		}

		// --------------------------------------------------
		//  -b
		// --------------------------------------------------
		idx.push_back(num_rows + 1 + N);
		coef.push_back(-1);
		nzcnt++;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		idx.clear();
		coef.clear();
		v_0 = num_rows + N + 2;

	}

	// constraints A_T * v + e_k * v_0

	{
		std::vector<int> idx;
		std::vector<double> coef;

		// --------------------------------------------------
		//  A_T * v
		// --------------------------------------------------
		for (int i = 0; i < N; i++) {
			char sense = 'L';
			int matbeg = 0;
			double rhs = 0;
			int nzcnt = 0;

			int iter = 0;
			int v = v_0;
			v++;

			while (iter < num_constraint) {

				if (A[iter][i] != 0) {
					idx.push_back(v);
					coef.push_back(A[iter][i]);
					nzcnt++;
				}
				iter++;
				v++;
			}


			// --------------------------------------------------
			//  e_k * v0
			// --------------------------------------------------
			if (i == k) {
				idx.push_back(v_0);
				coef.push_back(1);
				nzcnt++;
			}

			// --------------------------------------------------
			//  -a_i
			// --------------------------------------------------
			idx.push_back(num_rows + 1 + i);
			coef.push_back(-1);
			nzcnt++;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

	//constraint b_T * v + v_0 * (gamma+1)

	{
		std::vector<int> idx;
		std::vector<double> coef;

		char sense = 'G';
		int matbeg = 0;
		double rhs = 0;
		int nzcnt = 0;
		int v = v_0;
		v++;


		// --------------------------------------------------
		//  b_T * v
		// --------------------------------------------------
		for (int i = 0; i < num_constraint; i++) {

			if (b[i] != 0) {
				idx.push_back(v);
				coef.push_back(b[i]);
				nzcnt++;
			}
			v++;

		}


		// --------------------------------------------------
		//  (gamma+1) * v0
		// --------------------------------------------------
		if (gam != 0) {
			idx.push_back(v_0);
			coef.push_back(gam + 1);
			nzcnt++;
		}

		// --------------------------------------------------
		//  -b
		// --------------------------------------------------
		idx.push_back(num_rows + N + 1);
		coef.push_back(-1);
		nzcnt++;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		idx.clear();
		coef.clear();

	}

	// A_i* a + b_i* beta + u_i + v_i  =  A_i*c + b_i*z + u_i + v_i

	{

		std::vector<int> idx;
		std::vector<double> coef;

		//A_i*a
		for (int i = 0; i < num_constraint; i++) {

			double r = 0;
			int nzcnt = 0;

			for (int iter = 0; iter < N; iter++) {
				r += A[i][iter] * c[iter];
				if (A[i][iter] != 0) {
					idx.push_back(num_rows + 1 + iter);
					coef.push_back(A[i][iter]);
					nzcnt++;
				}

			}
			r += b[i] * min_sol + dual_varVals_P1[i] + dual_varVals_P2[i];

			//b[i]*beta
			if (b[i] != 0) {
				idx.push_back(num_rows + 1 + N);
				coef.push_back(b[i]);
				nzcnt++;
			}

			//u_i
			idx.push_back(1 + i);
			coef.push_back(1);
			nzcnt++;

			//v_i
			idx.push_back(num_rows + 3 + N + i);
			coef.push_back(1);
			nzcnt++;

			char sense = 'E';
			int matbeg = 0;
			double rhs = r;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}



	}

}


std::vector<double> SecondProblem::evaluate_rT() {

	std::vector<double> rt;
	double sum = 0;

	// --------------------------------------------------
	//evaluate A_T + e_k -1;
	// --------------------------------------------------
	for (int i = 0; i < N; i++) {

		sum = 0;
		for (int j = 0; j < num_constraint; j++)
			sum += A[j][i];



		// --------------------------------------------------
		//  e_k
		// --------------------------------------------------
		if (i == k)
			sum++;

		// --------------------------------------------------
		//  -1 (coeff a)
		// --------------------------------------------------
		sum--;

		rt.push_back(sum);

	}

	// --------------------------------------------------
	//evaluate b_T + gamma -1;
	// --------------------------------------------------
	sum = 0;
	for (int i = 0; i < num_constraint; i++)
		sum += b[i];


	// --------------------------------------------------
	//  gamma
	// --------------------------------------------------
	sum += gam;

	// --------------------------------------------------
	//  - 1 (coeff b)
	// --------------------------------------------------
	sum--;

	rt.push_back(sum);

	// --------------------------------------------------
	//duplicate vector for other constraints with gamma+1
	// --------------------------------------------------
	std::vector<double> temp = rt;
	temp[temp.size() - 1]++;

	rt.insert(rt.end(), temp.begin(), temp.end());

	cout << "vettore r duplicato" << endl;
	for (std::vector<double>::const_iterator j = rt.begin(); j != rt.end(); ++j)
		cout << *j << " ";
	cout << endl;

	// --------------------------------------------------
	// A + b + 1 + 1 (the last constraints)
	// --------------------------------------------------
	for (int i = 0; i < num_constraint; i++) {
		sum = 0;

		for (int j = 0; j < N; j++)
			sum += A[i][j];

		sum += b[i];

		// u
		sum++;

		// v
		sum++;

		rt.push_back(sum);

	}


	cout << "vettore r " << endl;
	for (std::vector<double>::const_iterator j = rt.begin(); j != rt.end(); ++j)
		cout << *j << " ";
	cout << endl;

	return rt;
}

void SecondProblem::set_solution(CEnv env, Prob lp, int num_rows){

		vector<double> varibles;

		cout << "!!!!!!!  VARIABILI DEL PROBLEMA DUALE: " << endl;
		int cur_numcols = CPXgetnumcols(env, lp);

		varibles.resize(cur_numcols);
		CHECKED_CPX_CALL(CPXgetx, env, lp, &varibles[0], 0, cur_numcols - 1);

		int surplus;
		status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0,
				cur_numcols - 1);
		int cur_colnamespace = -surplus; // the space needed to save the names

		// allocate memory
		char** cur_colname = (char **) malloc(sizeof(char *) * cur_numcols);
		char* cur_colnamestore = (char *) malloc(cur_colnamespace);

		// get the names
		CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace,
				&surplus, 0, cur_numcols - 1);


		//  print index, name and value of each column
		for (int i = 0; i < cur_numcols; i++) {
			cout << cur_colname[i] << " = " << varibles[i] << endl;
		}

		cout << endl;
		cout << endl;
		//  set variables
		u0 = varibles[0];
		cout << cur_colname[0] << " = " << u0 << endl;

		for (int i = 1; i <= num_rows; i++) {
			cout << cur_colname[i] << " = " << varibles[i] << endl;
			u.push_back(varibles[i]);
		}

		for (int i = num_rows + 1; i < num_rows + 1+N; i++) {
			cout << cur_colname[i] << " = " << varibles[i] << endl;
			a.push_back(varibles[i]);
		}

		beta = varibles[num_rows + N + 1];
		cout << cur_colname[num_rows + N + 1] << " = " << beta << endl;


		v0 = varibles[num_rows + N + 2];
		cout << cur_colname[num_rows + N + 2] << " = " << v0 << endl;

		for (int i = num_rows + N + 3; i < 2*num_rows + N + 3; i++) {
					cout << cur_colname[i] << " = " << varibles[i] << endl;
					v.push_back(varibles[i]);
		}


		// free
		free(cur_colname);
		free(cur_colnamestore);


}

void SecondProblem::solve(CEnv env, Prob lp, int num_rows){

	CHECKED_CPX_CALL(CPXlpopt, env, lp);
	print_objval(env, lp);

	set_solution(env, lp, num_rows);


}
