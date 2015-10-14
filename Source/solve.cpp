/**
 * @file solve.cpp
 * solve file
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 1.0
 */

#include "solve.h"
#include <unistd.h>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;


/**
 Method that create P1 problem (make a branch of admissible region),
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return none
 */
void create_P1_prob(CEnv env, Prob lp, int index){

		cout << endl;
		cout << "PROBLEM P1" << endl;

		double rhs = floor(varVals[index]);

		char sense = 'L';
		int matbeg = 0;
		const int idx = index;
		const double coef = 1;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
				&coef, 0, 0);

		// print P1 problem to a file
		//CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp1", 0);

}
/**
 Method that create P2 problem (make a branch of admissible region),
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return none
 */
void create_P2_prob(CEnv env, Prob lp, int index) {

	cout << endl;
	cout << "PROBLEM P2" << endl;

	double rhs = floor(varVals[index]) + 1;
	char sense = 'G';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	//CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp2", 0);

}

/**
 Method that create P1 sub_problem (make a branch of admissible region),
 return solution if exist, otherwise return INFINITE value
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return double*, pointer to vector[2] where the first element is
 the result of P1 sub_problem, the second is the result of sub_P2 problem
 */
double* solve_P1_Problem(CEnv env, Prob lp, int index) {

	static double z[2];
	z[0] = CPX_INFBOUND;
	z[1] = CPX_INFBOUND;

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	CHECKED_CPX_CALL(CPXrefineconflict, env, lp, NULL, NULL);

	int stat = CPXgetstat(env, lp);
	int cur_numrows = CPXgetnumrows(env, lp);

// print and set solution and create and resolve P_2 problem"
	if (stat == CPX_STAT_CONFLICT_FEASIBLE) {
		cout << "FEASIBLE " << endl;
		gam = floor(varVals[index]);
		print_objval(env, lp);
		set_and_print_var_P(env, lp);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z[0]);
		set_and_print_var_D(env, lp, true);
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);

		create_P2_prob(env, lp, index);
		z[1] = solve_P2_Problem(env, lp, index);

	} else {
		cout << "No solution for P1 problem exists.. " << endl;

		double rhs = floor(varVals[index]) + 1;
		char sense = 'G';
		int matbeg = 0;
		const int idx = index;
		const double coef = 1;


		double* cut = new double[N];

		for (int i = 0; i < N; i++) {
			if (i == index) {
				cut[i] = 1;
			} else
				cut[i] = 0;
		}


		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg,
				&idx, &coef, 0, 0);

		cut_A.push_back(cut);
		cut_b.push_back(rhs);


		cout << "Resolve a new problem P1.. " << endl;
		solve(env, lp);

	}

	return z;

}

/**
 Method that create P2 problem (make a branch of admissible region),
 return solution if exist, otherwise return INFINITE value
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return double , result of P2 sub_problem
 */
double solve_P2_Problem(CEnv env, Prob lp, int index) {

	double z = CPX_INFBOUND;

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	CHECKED_CPX_CALL(CPXrefineconflict, env, lp, NULL, NULL);

	int stat = CPXgetstat(env, lp);

	int cur_numrows = CPXgetnumrows(env, lp);

	if (stat == CPX_STAT_CONFLICT_FEASIBLE) {
		cout << "FEASIBLE " << endl;
		print_objval(env, lp);
		set_and_print_var_P(env, lp);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z);
		set_and_print_var_D(env, lp, false);

	} else {

		cout << "No solution for P2 problem exists " << endl;
		double rhs = floor(varVals[index]);

		char sense = 'L';
		int matbeg = 0;
		const int idx = index;
		const double coef = 1;

		double* cut;
		cut = new double[N];

		for (int i = 0; i < N; i++) {
			if (i == index) {
				cut[i] = 1;
			} else
				cut[i] = 0;
		}

		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg,
				&idx, &coef, 0, 0);

		cut_A.push_back(cut);
		cut_b.push_back(rhs);


		cout << "Resolve a new problem P2.. " << endl;
		solve(env, lp);
	}

	return z;
}
/**
 Method solve
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return void
 */
void solve(CEnv env, Prob lp) {

	// --------------------------------------------------
	// 3. solve linear problem
	// --------------------------------------------------
	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	// --------------------------------------------------
	// 4. print solution
	// --------------------------------------------------
	print_objval(env, lp);

	// --------------------------------------------------
	// 5. set number and value of variable
	//    (cur_numcols,varVals) and print these
	// --------------------------------------------------
	set_and_print_var_P(env, lp);

	// --------------------------------------------------
	// 6. chose the best fractional variable
	// --------------------------------------------------
	int index = select_fractionar_var(varVals);


	//usleep(1000000);
	// --------------------------------------------------------
	// 7. if x solution aren't integer create P1 and P2 problem
	// --------------------------------------------------------
	if (index != -1) {

		cout << endl;
		cout << "Higher fractional variable choose " << varVals[index] << endl;

		cout << "Index " << index << endl;

		create_P1_prob(env, lp, index);

		double* z = solve_P1_Problem(env, lp, index);
		/////////////////////////////////////////////////////

		if (*z < CPX_INFBOUND && *(z + 1) < CPX_INFBOUND && flag_find) {
			flag_find = false;

			CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
			// --------------------------------------------------
			// print solution in standard format
			// --------------------------------------------------
			CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem.sol");

			min_sol = std::min(*z, *(z + 1));
			k=index;
			cout << "k = " << k << endl;
			cout << "objective function lesser = " << min_sol << endl;


		}

	}else{
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		cout << "The last solution is the best integer solution. STOP" << endl;
		CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem.sol");
	}

}

void print_cut_A() {
	cout << "cut: " << endl;
	for (std::vector<double*>::const_iterator i = cut_A.begin();
			i != cut_A.end(); ++i) {
		double*ptr = *i;
		for (int j = 0; j < N; ++j) {
			std::cout << ptr[j] << ' ';
		}
		cout << endl;
	}
}
void print_cut_b() {
	cout << "vector_b: " << endl;
	for (std::vector<double>::const_iterator i = cut_b.begin();
			i != cut_b.end(); ++i) {
		std::cout << *i << ' ';
		cout << endl;
	}
}

void print_u_variables() {
	cout << "u variable: " << endl;
	for (std::vector<double>::const_iterator i = dual_varVals_P1.begin();
			i != dual_varVals_P1.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}

void print_v_variables() {

	cout << "v variable: " << endl;
	for (std::vector<double>::const_iterator i =
			dual_varVals_P2.begin(); i != dual_varVals_P2.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}


