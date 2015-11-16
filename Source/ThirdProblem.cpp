/*
 * Thirdproblem.cpp
 *
 *  Created on: 14 nov 2015
 *      Author: riccardo
 */

#include "ThirdProblem.h"

ThirdProblem::ThirdProblem(vector<double> y_til) {
	// TODO Auto-generated constructor stub
	this->y_tilde = y_til;
	y_tilde_MIN_y_bar();

}

ThirdProblem::~ThirdProblem() {
	// TODO Auto-generated destructor stub
}

void ThirdProblem::y_tilde_MIN_y_bar() {

	double difference;

	// --------------------------------------------------
	// 1. c-a
	// --------------------------------------------------
	cout << endl;

	cout << "Comparison y bar and y til: " << endl;

	cout << "vector c = ";
	for (unsigned int i = 0; i < c.size(); i++)
		cout << c[i] << " ";

	cout << endl;

	cout << "vector a = ";
	for (unsigned int i = 0; i < c.size(); i++)
		cout << y_tilde[i] << " ";

	cout << endl << "\t     ";

	int j = 0;
	for (unsigned int i = 0; i < c.size(); i++) {
		difference = c[i] - y_tilde[j];

		cout << difference << " ";
		//tolerance error
		if (difference < epsilon && difference > -epsilon)
			difference = 0.0;

		t.push_back(difference);
		j++;
	}

	cout << endl;

	// --------------------------------------------------
	// 2. z-b
	// --------------------------------------------------

	cout << "z = " << min_sol << endl;
	cout << "y_tilde = " << y_tilde[j];
	cout << endl;

	difference = min_sol - y_tilde[j];

	cout << "difference = " << difference << endl;

	//tolerance error
	if (difference < epsilon && difference > -epsilon)
		difference = 0.0;

	t.push_back(difference);
	j++;

	// --------------------------------------------------
	// 3. u-u && u_0-u_0
	// --------------------------------------------------
	cout << endl;
	cout << "vector u (y bar)= ";
	for (unsigned int i = 0; i < dual_varVals_P1.size(); i++)
		cout << dual_varVals_P1[i] << " ";

	cout << endl;

	cout << "vector u (y til)= ";
	for (unsigned int i = j; i < dual_varVals_P1.size() + j; i++)
		cout << y_tilde[i] << " ";

	cout << endl << "\t\t";

	for (unsigned int i = 0; i < dual_varVals_P1.size(); i++) {
		difference = dual_varVals_P1[i] - y_tilde[j];

		cout << difference << " ";
		//tolerance error
		if (difference < epsilon && difference > -epsilon)
			difference = 0.0;

		t.push_back(difference);
		j++;
	}

	// --------------------------------------------------
	// 5. v-v && v_0 - v_0
	// --------------------------------------------------
	cout << endl;
	cout << "vector v (y bar)= ";
	for (unsigned int i = 0; i < dual_varVals_P2.size(); i++)
		cout << dual_varVals_P2[i] << " ";

	cout << endl;

	cout << "vector v (y til)= ";
	for (unsigned int i = j; i < dual_varVals_P2.size() + j; i++)
		cout << y_tilde[i] << " ";

	cout << endl << "\t\t";

	for (unsigned int i = 0; i < dual_varVals_P2.size(); i++) {
		difference = dual_varVals_P2[i] - y_tilde[j];
		cout << difference << " ";

		//tolerance error
		if (difference < epsilon && difference > -epsilon)
			difference = 0.0;

		t.push_back(difference);
		j++;
	}

}

void ThirdProblem::setup(CEnv env, Prob lp) {

	{
		static const char* varType = NULL;
		double obj = -1.0;
		double lb = -CPX_INFBOUND;
		double ub = CPX_INFBOUND;

		snprintf(name, NAME_SIZE, "lambda");
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

	}

	{

		// constraints

		std::vector<int> idx;
		std::vector<double> coef;

		char sense = 'G';
		int matbeg = 0;
		double rhs;
		int nzcnt = 1;
		idx.push_back(0);

		double cof;

		{
			//per ogni riga della matrice (1 vincolo)
			for (int i = 0; i < N; i++) {
				cof = 0;
				rhs = 0;
				//1* y_tilde[a] && 1*t[a]
				int j = 0;
				while (j < N) {
					cof += t[j];
					rhs += y_tilde[j];
					j++;
				}

				//beta are 0... skip
				j++;

				//A_T * y_tilde[u] && A_T * t[u]
				for (int iter = 0; iter < num_constraint; iter++) {
					cof += A[iter][i] * t[j];
					rhs += A[iter][i] * y_tilde[j];
					j++;
				}

				if (i == k) {
					cof -= t[j];
					rhs -= y_tilde[j];
				}
				//v and v_0 are 0.. skip

				rhs = -rhs;

				coef.push_back(cof);

				//add constraints
				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);

				idx.clear();
				coef.clear();
			}
		}
		{	//2 vincolo
			cof = 0;
			rhs = 0;

			//a is 0...skip
			int j = N;

			//-beta
			cof -= t[j];
			rhs -= y_tilde[j];

			j++;

			//u
			for (int i = 0; i < num_constraint; i++) {
				cof += b[i] * t[j];
				rhs += b[i] * y_tilde[j];

				j++;
			}

			//u_0
			cof += gam * t[j];
			rhs += gam * y_tilde[j];

			//v and v_0 are 0... skip
			rhs = -rhs;

			coef.push_back(cof);

			//add constraints
			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}

		{
			//per ogni riga della matrice (3 vincolo)
			for (int i = 0; i < N; i++) {
				cof = 0;
				rhs = 0;
				//1* y_tilde[a] && 1*t[a]
				int j = 0;
				while (j < N) {
					cof += t[j];
					rhs += y_tilde[j];
					j++;
				}

				//beta are 0... skip
				j++;

				//u && u_0 are 0... skip
				j += num_constraint;

				//A_T * y_tilde[v] && A_T * t[v]
				for (int iter = 0; iter < num_constraint; iter++) {
					cof += A[iter][i] * t[j];
					rhs += A[iter][i] * y_tilde[j];
					j++;
				}

				//v_0
				if (i == k) {
					cof -= t[j];
					rhs -= y_tilde[j];
				}

				rhs = -rhs;

				coef.push_back(cof);
				//add constraints
				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);

				idx.clear();
				coef.clear();
			}

		}

		{
			//4 vincolo
			cof = 0;
			rhs = 0;

			//a is 0...skip
			int j = N;

			//-beta
			cof -= t[j];
			rhs -= y_tilde[j];

			j++;

			//u and u_0 are 0... skip
			j += num_constraint;

			//b_t * t[v]
			for (int i = 0; i < num_constraint; i++) {
				cof += b[i] * t[j];
				rhs += b[i] * y_tilde[j];

				j++;
			}

			//v_0
			cof += gam * t[j];
			rhs += gam * y_tilde[j];

			rhs = -rhs;

			coef.push_back(cof);
			//add constraints
			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();

		}
	}

}
