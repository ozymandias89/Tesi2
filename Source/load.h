/**
 * @file load.h
 * library for support
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 1.0
 */

#ifndef SOURCE_LOAD_H_
#define SOURCE_LOAD_H_

// Includes:
#include <cstdio>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <cmath>

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

//number of constraint
int num_constraint;

//Coefficient cost
double* c;

//price matrix
double** A;

//known terms
double* b;

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

	if (!myfile.is_open())
		cerr << "Unable to open file!!!!!";

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

			c = new double[N];

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

			b = new double[num_constraint];

			std::istringstream it(line);

			int i = 0;
			while (it >> number) {

				b[i] = atof(number.c_str());
				i++;
			}

			flag_known = false;
		}
	} while (flag_known);

	//cout << endl;
	//cout << "var e vincoli " << N << " " << num_constraint;

	A = new double*[num_constraint];
	for (int i = 0; i < num_constraint; ++i)
		A[i] = new double[N];

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
 Method that set the problem in this case only one constraint
 @param  (CEnv env, Prob lp)
 @return void
 */
void setupLP(CEnv env, Prob lp) {

	{
		char varType = 'C';
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;
		for (int i = 0; i < N; i++) {
			obj = c[i];
			snprintf(name, NAME_SIZE, "x_%i", i);
			char* varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, &varType,
					&varName);

		}

	}

	// constraint

	{
		std::vector<int> idx;
		std::vector<double> coef;

		for (int i = 0; i < num_constraint; i++) {
			char sense = 'E';
			int matbeg = 0;
			double rhs = b[i];

			for (int iter = 0; iter < N; iter++) {

				if (A[i][iter] != 0) {
					idx.push_back(iter);
					coef.push_back(A[i][iter]);
				}

			}

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 6, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

}
/**
 Method that choise the best fractional variable
 @param  (vector<double>)
 @return int, return index of highter variable functional (-1 if no variable is fractional)
 */
int fractionar_variable(std::vector<double> varVals) {

	double threesold = 0.5;
	double temp;
	double temp2;
	int index = -1;

	for (unsigned int i = 0; i < varVals.size(); ++i) {

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
 Method that print matrix
 @param  none
 @return void
 */
void print_matrix() {

	cout << "Matrix A" << endl;
	for (int i = 0; i < num_constraint; i++) {
		for (int j = 0; j < N; j++) {

			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

/**
 Method that print vector_c
 @param  none
 @return void
 */
void print_vect_c() {
	cout << "Vector c" << endl;

	for (int i = 0; i < N; i++)
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
 @param  (CEnv env, Prob lp), environmant of the problem and problem
 @return void
 */
void print_objval(CEnv env, Prob lp) {

	double objval;
	CHECKED_CPX_CALL(CPXgetobjval, env, lp, &objval);
	std::cout << "Obj val: " << objval << std::endl;
	cout << endl;

}

#endif /* SOURCE_LOAD_H_ */
