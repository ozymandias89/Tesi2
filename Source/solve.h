/**
 * @file solve.h
 * library for support
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 1.0
 */

#ifndef SOURCE_SOLVE_H_
#define SOURCE_SOLVE_H_

#include "load.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

/**
 Method that create P1 problem (make a branch of admissible region),
 return solution if exist, otherwise return INFINITE value
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return double
 */
double solve_P1_Problem(CEnv env, Prob lp, int index) {

	double z = CPX_INFBOUND;
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
	CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp1", 0);

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	CHECKED_CPX_CALL(CPXrefineconflict, env, lp, NULL, NULL);

	int stat = CPXgetstat(env, lp);
	cout << "Status of the problem: " << stat << endl;
	//CHECKED_CPX_CALL(CPXclpwrite, env, lp, "../data/conflict.lp1" );

	// print and set solution
	if (stat == 30) {
		print_objval(env, lp);
		print_var_P(env, lp);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z);
		print_var_D(env, lp, true);
	} else
		cout << "No solution for P1 problem exists " << endl;


	int cur_numrows = CPXgetnumrows(env, lp);


	CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);


	return z;

}

/**
 Method that create P2 problem (make a branch of admissible region),
 return solution if exist, otherwise return INFINITE value
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return void
 */
double solve_P2_Problem(CEnv env, Prob lp, int index) {

	double z = CPX_INFBOUND;
	cout << endl;
	cout << "PROBLEM P2" << endl;

	double rhs = floor(varVals[index]) + 1;
	char sense = 'G';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp2", 0);

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	CHECKED_CPX_CALL(CPXrefineconflict, env, lp, NULL, NULL);
	//CHECKED_CPX_CALL(CPXclpwrite, env, lp, "../data/conflict.lp2" );
	//CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem2.sol");

	int stat = CPXgetstat(env, lp);
	cout << "Status of the problem: " << stat << endl;

	if (stat == 30) {
		print_objval(env, lp);
		print_var_P(env, lp);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z);
		print_var_D(env, lp, false);
	} else
		cout << "No solution for P2 problem exists " << endl;

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
	// 2. write problem
	// --------------------------------------------------
	CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);

	// --------------------------------------------------
	// 3. Optimization
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
	print_var_P(env, lp);


	// --------------------------------------------------
	// 6. chose the best fractional variable
	// --------------------------------------------------
	int index = select_fractionar_var(varVals);

	// --------------------------------------------------------
	// 7. if x solution aren't integer create P1 and P2 problem
	// --------------------------------------------------------
	if (index != -1) {

		cout << endl;
		cout << "Higher fractional variable choose " << varVals[index] << endl;

		cout << "Index " << index << endl;

		double z1 = solve_P1_Problem(env, lp, index);
		/////////////////////////////////////////////////////
		double z2 = solve_P2_Problem(env, lp, index);

		if (z1<CPX_INFBOUND && z2<CPX_INFBOUND){
			double z = std::min (z1,z2);
			cout << "funzione obbiettivo minore = " << z << endl;
		}

	}

}


#endif /* SOURCE_SOLVE_H_ */
