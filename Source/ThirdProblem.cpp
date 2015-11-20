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

	int j = 0;
	for (unsigned int i = 0; i < c.size(); i++) {
		difference = c[i] - y_tilde[j];
		t.push_back(difference);
		j++;
	}

	// --------------------------------------------------
	// 2. z-b
	// --------------------------------------------------

	difference = min_sol - y_tilde[j];

	t.push_back(difference);
	j++;

	// --------------------------------------------------
	// 3. u-u && u_0-u_0
	// --------------------------------------------------

	for (unsigned int i = 0; i < dual_varVals_P1.size(); i++) {
		difference = dual_varVals_P1[i] - y_tilde[j];
		t.push_back(difference);
		j++;
	}

	// --------------------------------------------------
	// 5. v-v && v_0 - v_0
	// --------------------------------------------------

	for (unsigned int i = 0; i < dual_varVals_P2.size(); i++) {
		difference = dual_varVals_P2[i] - y_tilde[j];
		t.push_back(difference);
		j++;
	}

}

void ThirdProblem::setup(CEnv env, Prob lp) {

	{
		// lambda variable
		static const char* varType = NULL;
		double obj = 1.0;
		double lb = -CPX_INFBOUND;
		double ub = CPX_INFBOUND;

		snprintf(name, NAME_SIZE, "lambda");
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		CPXchgobjsen (env, lp, CPX_MAX);
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
			//first constraint (second problem)
			for (int i = 0; i < N; i++) {
				cof = 0;
				rhs = 0;

				//1* y_tilde[a] && 1*t[a]
				cof += t[i];
				rhs += y_tilde[i];

				int j = N;
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

				//tolerance error
				if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
					rhs = 0.0;

				if (cof < epsilon_8_4 && cof > -epsilon_8_4)
					cof = 0.0;

				coef.push_back(cof);

				//add constraints
				//if (cof!=0)
				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);

				idx.clear();
				coef.clear();
			}
		}
		{	//Second constraint
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

//			//tolerance error
			if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
				rhs = 0.0;

			if (cof < epsilon_8_4 && cof > -epsilon_8_4)
				cof = 0.0;

			coef.push_back(cof);

			//add constraints
//			if (cof!=0)
			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}

		{
			//for each lines of matrix (third constraint)
			for (int i = 0; i < N; i++) {

				cof = 0;
				rhs = 0;
				//1* y_tilde[a] && 1*t[a]

				cof += t[i];
				rhs += y_tilde[i];


				int j = N;
				//beta are 0... skip
				j++;

				//u && u_0 are 0... skip
				j += num_constraint + 1;

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

//				//tolerance error
				if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
					rhs = 0.0;

				if (cof < epsilon_8_4 && cof > -epsilon_8_4)
					cof = 0.0;

				coef.push_back(cof);

				//add constraints
//				if (cof != 0)
				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);

				idx.clear();
				coef.clear();
			}

		}

		{
			//fourth constraint
			cof = 0;
			rhs = 0;

			//a is 0...skip
			int j = N;

			//-beta
			cof -= t[j];
			rhs -= y_tilde[j];
			j++;

			//u and u_0 are 0... skip
			j += num_constraint + 1;

			//b_t * t[v]
			for (int i = 0; i < num_constraint; i++) {
				cof += b[i] * t[j];
				rhs += b[i] * y_tilde[j];
				j++;
			}

			//v_0
			cof += (gam + 1) * t[j];
			rhs += (gam + 1) * y_tilde[j];
			rhs = -rhs;

//			//tolerance error
			if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
				rhs = 0.0;

			if (cof < epsilon_8_4 && cof > -epsilon_8_4)
				cof = 0.0;

			coef.push_back(cof);

			//add constraints
			//		if (cof != 0)
			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();

		}

		{
			// constraint -u_0 >= 0
			cof = 0;
			rhs = 0;

			int j = N + 1 + num_constraint;
			//-u
			cof -= t[j];
			rhs -= y_tilde[j];

			rhs = -rhs;
			//tolerance error
			if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
				rhs = 0.0;

			if (cof < epsilon_8_4 && cof > -epsilon_8_4)
				cof = 0.0;

			coef.push_back(cof);

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();

		}

		{
			// constraint v_0 >= 0
			cof = 0;
			rhs = 0;

			//v
			cof += t.back();
			rhs += y_tilde.back();

			rhs = -rhs;
			//tolerance error
			if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
				rhs = 0.0;

			if (cof < epsilon_8_4 && cof > -epsilon_8_4)
				cof = 0.0;
			coef.push_back(cof);

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();

		}
	}

	CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/third_problem.lp", 0);

}

void ThirdProblem::update_y_bar(CEnv env, Prob lp) {

	CHECKED_CPX_CALL(CPXrefineconflict, env, lp, NULL, NULL);
	int stat = CPXgetstat(env, lp);

	if (stat == CPX_STAT_CONFLICT_FEASIBLE) {
		CHECKED_CPX_CALL(CPXprimopt, env, lp);

//		print_objval(env, lp);

		vector<double> varibles;

		cout << "VARIABLES THIRD PROBLEM: " << endl;
		int cur_numcols = 1;

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

		//  print index, name and value of lambda
		cout << cur_colname[0] << " = " << varibles[0] << endl;

		double lambda = varibles[0];

		// free
		free(cur_colname);
		free(cur_colnamestore);


		vector <double> r_mul_lamb;

		cout << endl;

		//-------------------------------------------------
		// lambda * (y_bar-y_tilde)
		//-------------------------------------------------

		for (vector<double>::iterator it = t.begin() ; it != t.end(); ++it){
			r_mul_lamb.push_back((*it * lambda));
		}

		//-------------------------------------------------
		// y_tilde + r_mul_lamb
		//-------------------------------------------------

		int j=0;

		//c
		for (unsigned int i=0; i < c.size(); ++i) {
			c[j] = y_tilde[j] + r_mul_lamb[j];
			j++;
		}

		//z
		min_sol= y_tilde[j] + r_mul_lamb[j];
		j++;


		//u and u_0
		for (unsigned int i=0; i < dual_varVals_P1.size(); ++i){
			dual_varVals_P1[i] = y_tilde[j] + r_mul_lamb[j];
			j++;
		}

		cout << endl;

		//v and v_0
		for (unsigned int i = 0; i < dual_varVals_P2.size(); ++i) {
			dual_varVals_P2[i] = y_tilde[j] + r_mul_lamb[j];
			j++;
		}

		print_vect_c();
		cout << " beta " << min_sol <<endl;
		print_u_variables();
		print_v_variables();

	} else {
		cerr << "Third  problem has conflict!!!!!" << endl;
		exit(1);
	}

}
