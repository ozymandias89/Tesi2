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

//index for e_k
int k;

//number of constraint
int num_constraint;

//Coefficient cost
double* c;

//price matrix
double** A;

//known terms
double* b;

//primary variables
std::vector<double> varVals;

//dual variables P_1 problem
std::vector<double> dual_varVals_P1;

//dual variables P_2 problem
std::vector<double> dual_varVals_P2;

//inequality insert in cut_A
std::vector<double*> cut_A;
//inequality insert in cut_b
std::vector<double> cut_b;

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
 Method that set the second problem
 @param  (CEnv env, Prob lp)
 @return void
 */
void setupSP(CEnv env, Prob lp, int num_rows, int num_cols) {

	{
		// variables
		static const char* varType = NULL;
		double obj = 1.0;
		double lb = -CPX_INFBOUND;
		double ub = 0.0;
		snprintf(name, NAME_SIZE, "u_%i", 0);
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,&varName);

		ub = CPX_INFBOUND;

		for (int i = 1; i < num_rows; i++) {
			snprintf(name, NAME_SIZE, "u_%i", i);
			varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType, &varName);
		}

		lb = 0.0;
		snprintf(name, NAME_SIZE, "v_%i", 0);
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		lb = -CPX_INFBOUND;

		for (int i = 1; i < num_rows; i++) {
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
			double rhs = c[i];
			int nzcnt=0;

			int iter = 0;
			int u=1;

			while ( iter < num_constraint) {

				//cout << A[iter][i];
				if (A[iter][i] != 0) {
					idx.push_back(u);
					coef.push_back(A[iter][i]);
					nzcnt++;
				}
				 iter++;
				 u++;
			}//cout << endl;

			for (std::vector<double*>::const_iterator j = cut_A.begin();
					j != cut_A.end(); ++j) {
				double*ptr = *j;
				//cout <<ptr[i];
				if (ptr[i] != 0) {
					idx.push_back(u);
					coef.push_back(ptr[i]);
					nzcnt++;
				}
				iter++;
				u++;
			}//cout << endl;

			// --------------------------------------------------
			//  e_k * u
			// --------------------------------------------------
			if (i == k) {
				idx.push_back(0);
				coef.push_back(1);
				nzcnt++;
			}


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
		double rhs = min_sol;
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

		for (std::vector<double>::const_iterator j = cut_b.begin();
				j != cut_b.end(); ++j) {
			//cout << *j;
			if (*j != 0) {
				idx.push_back(u);
				coef.push_back(*j);
				nzcnt++;
			}
			u++;
		}
		//cout << endl;

		if (gam != 0) {
			idx.push_back(0);
			coef.push_back(gam);
			nzcnt++;
		}

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		idx.clear();
		coef.clear();
		v_0=u;

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
				double rhs = c[i];
				int nzcnt = 0;

				int iter = 0;
				int v = v_0;
				v++;

				while ( iter < num_constraint) {

					//cout << A[iter][i];
					if (A[iter][i] != 0) {
						idx.push_back(v);
						coef.push_back(A[iter][i]);
						nzcnt++;
					}
					 iter++;
					 v++;
				}//cout << endl;


				for (std::vector<double*>::const_iterator j = cut_A.begin();
						j != cut_A.end(); ++j) {
					double*ptr = *j;
					//cout <<ptr[i];
					if (ptr[i] != 0) {
						idx.push_back(v);
						coef.push_back(ptr[i]);
						nzcnt++;
					}
					iter++;
					v++;
				}//cout << endl;


				// --------------------------------------------------
				//  e_k * u
				// --------------------------------------------------
				if (i == k) {
					idx.push_back(v_0);
					coef.push_back(1);
					nzcnt++;
				}


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
			double rhs = min_sol;
			int nzcnt = 0;
			int v = v_0;
			v++;


			for (int i = 0; i < num_constraint; i++) {

				if (b[i] != 0) {
					idx.push_back(v);
					coef.push_back(b[i]);
					nzcnt++;
				}
				v++;

			}

			for (std::vector<double>::const_iterator j = cut_b.begin();
					j != cut_b.end(); ++j) {
				//cout << *j;
				if (*j != 0) {
					idx.push_back(v);
					coef.push_back(*j);
					nzcnt++;
				}
				v++;
			}
			//cout << endl;

			if (gam != 0) {
				idx.push_back(0);
				coef.push_back(gam+1);
				nzcnt++;
			}

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();

		}


}

/**
 Method that chooses the best fractional variable
 @param  (vector<double>)
 @return int, return index of higher variable functional (-1 if no variable is fractional)
 */
int select_fractionar_var(std::vector<double> varVals) {

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
 Method that print matrix A
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

	cout << "DUAL VARIABLES: " << endl;
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
