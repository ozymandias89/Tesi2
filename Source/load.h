/**
 * @file load.h
 * library for support
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 2.0
 */

#ifndef SOURCE_LOAD_H_
#define SOURCE_LOAD_H_

// Includes:
#include <cstdio>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "cpxmacro.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

// error status and message buffer
int status;
char errmsg[BUF_SIZE];

const int NAME_SIZE = 512;
char name[NAME_SIZE];

//number of variable
int N;

//number of variable original
int Num_original_variables;

//number of variable original
int Num_original_constraints;

//index for e_k
int k;

//number of constraint
int num_constraint;

//Coefficient cost
std::vector<double> c;

//price matrix
std::vector< std::vector<double> > A;

//known terms
std::vector<double> b;

//primary variables
std::vector<double> varVals;

//dual variables P_1 problem
std::vector<double> dual_varVals_P1;

//dual variables P_2 problem
std::vector<double> dual_varVals_P2;

//gamma
int gam;
//min_sol from P_1/P_2
double min_sol;

/**
 Method that load from file the problem (for example of format file see folder data)
 @param  (ifstream &) , object ifstream
 @return void
 */
void load_problem(ifstream &myfile) {

	bool flag_obb = true;
	bool flag_known = true;

	std::string line;
	std::string number;

	if (!myfile.is_open()) {
		cerr << "Unable to open file!!!!!" << endl;
		exit(1);
	}
	do {
		getline(myfile, line);
		//if aren't a comment or empty line
		if ((line[0] != '/' && line[1] != '/') && line.length() != 0) {

			std::istringstream is(line);

			int count = 0;
			while (is >> number) {
				count++;
			}

			N = count;
			Num_original_variables=N;

			c.resize(N);

			std::istringstream it(line);

			int i = 0;
			while (it >> number) {
				c[i] = atof(number.c_str());
				i++;
			}
			flag_obb = false;
		}

	} while (flag_obb);

	do {

		getline(myfile, line);
		//if aren't a comment or empty line
		if ((line[0] != '/' && line[1] != '/') && line.length() != 0) {

			std::istringstream is(line);
			int count = 0;
			while (is >> number) {
				count++;
			}

			num_constraint = count;
			Num_original_constraints=num_constraint;

			b.resize(num_constraint);

			std::istringstream it(line);

			int i = 0;
			while (it >> number) {

				b[i] = atof(number.c_str());
				i++;
			}

			flag_known = false;
		}
	} while (flag_known);


	//ALLOCATE MATRIX A
	 A.resize(num_constraint);
	for (int i = 0; i < num_constraint; ++i) {
		A[i].resize(N);
	}

	int i = 0;
	int j = 0;
	while (getline(myfile, line)) {
		//if aren't a comment or empty line
		if ((line[0] != '/' && line[1] != '/') && line.length() != 0) {
			std::istringstream it(line);
			j = 0;
			if (i >= num_constraint) {
				cerr << "Matrix don't respect standard line";
				exit(1);
			}
			while (it >> number) {
				if (j >= N) {
					cerr << "Matrix don't respect standard columns";
					exit(1);
				}
				A[i][j] = atof(number.c_str());
				j++;

			}
			i++;

		}
	}

}

/**
 Method that set the primary problem
 @param  (CEnv env, Prob lp)
 @return void
 */
void setupLP(CEnv env, Prob lp) {

	{	// variables
		static const char* varType = NULL;
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;

		for (int i = 0; i < N; i++) {
			obj = c[i];
			snprintf(name, NAME_SIZE, "x_%i", i);
			char* varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
					&varName);

		}

	}

	// constraints

	{
		std::vector<int> idx;
		std::vector<double> coef;

		for (int i = 0; i < num_constraint; i++) {
			char sense = 'E';
			int matbeg = 0;
			double rhs = b[i];
			int nzcnt=0;

			for (int iter = 0; iter < N; iter++) {

				if (A[i][iter] != 0) {
					idx.push_back(iter);
					coef.push_back(A[i][iter]);
					nzcnt++;
				}

			}

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

}


 /**
 Method that chooses the best x fractional variable
 @param  (vector<double>)
 @return int, return index of higher variable functional (-1 if no variable is fractional)
 */
int select_fractionar_var(std::vector<double> varVals) {

	double threesold = 0.5;
	double temp;
	double temp2;
	int index = -1;

	for (int i = 0; i < Num_original_variables; ++i) {

		temp = (varVals[i] - (int) varVals[i]);
		temp2 = fabs((temp - 0.5));

		if (temp2 < threesold) {
			threesold = temp2;
			index = i;
		}

	}
	return index;
}

/**
 Method that print matrix A
 @param  none
 @return void
 */
void print_matrix() {

	cout << endl;
	cout << "Matrix A " << endl;

	for (unsigned int i = 0; i < A.size(); i++) {
		for (unsigned int j = 0; j < A[i].size(); j++) {

			cout << A[i][j] << " ";
		}
		cout << endl;
	}

}

/**
 Method that print vector_c
 @param  none
 @return void
 */
void print_vect_c() {

	cout << "Vector c" << endl;

	for (int i = 0; i < Num_original_variables; i++)
		cout << c[i] << " ";

	cout << endl;
}

/**
 Method that print vector_b
 @param  none
 @return void
 */
void print_vect_b() {

	cout << "Vector b" << endl;
	for (int i = 0; i < num_constraint; i++)
		cout << b[i] << " ";

	cout << endl;
}
/**
 Method that print object function
 @param  (CEnv env, Prob lp), environment of the problem and problem
 @return void
 */
void print_objval(CEnv env, Prob lp) {
	cout << endl;
	double objval;
	CHECKED_CPX_CALL(CPXgetobjval, env, lp, &objval);
	std::cout << "Obj val: " << objval << std::endl;
	cout << endl;

}

/**
 Method that set and print primal variable
 @param  (CEnv env, Prob lp), environmant of the problem and problem
 @return void
 */
void set_and_print_var_P(CEnv env, Prob lp) {

	cout << "PRIMAL VARIABLES: " << endl;
	int cur_numcols = CPXgetnumcols(env, lp);

	varVals.resize(cur_numcols);
	CHECKED_CPX_CALL(CPXgetx, env, lp, &varVals[0], 0, cur_numcols - 1);

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
		cout << cur_colname[i] << " = " << varVals[i] << endl;
	}
	// free
	free(cur_colname);
	free(cur_colnamestore);
}

/**
 Method that set and print dual variables
 @param  (CEnv env, Prob lp, bool prob), environment of the problem,
  problem and flag(true P1 problem, false P2 problem)
 @return void
 */
void set_and_print_var_D(CEnv env, Prob lp, bool prob) {

	cout << endl;
	cout << "DUAL VARIABLES (the last is u_0 or v_0): " << endl;
	int num_rows = CPXgetnumrows(env, lp);

	if (prob) {
		dual_varVals_P1.resize(num_rows);
		CHECKED_CPX_CALL(CPXgetpi, env, lp, &dual_varVals_P1[0], 0,
				num_rows - 1);

		for (int i = 0; i < num_rows; i++) {
			cout << dual_varVals_P1[i]<< " ";
		}
		cout << endl;

	} else {
		dual_varVals_P2.resize(num_rows);
		CHECKED_CPX_CALL(CPXgetpi, env, lp, &dual_varVals_P2[0], 0,
				num_rows - 1);

		for (int i = 0; i < num_rows; i++) {
			cout << dual_varVals_P2[i] << " ";
		}
		cout << endl;
	}

}

#endif /* SOURCE_LOAD_H_ */
