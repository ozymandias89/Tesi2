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
		cout << "SUB_PROBLEM P1" << endl;

		double rhs = floor(varVals[index]);

		char sense = 'L';
		int matbeg = 0;
		const int idx = index;
		const double coef = 1;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
				&coef, 0, 0);

		cout << "insert inequality x_" << index << " <= " << rhs << endl;

}
/**
 Method that create P2 problem (make a branch of admissible region),
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return none
 */
void create_P2_prob(CEnv env, Prob lp, int index) {

	cout << endl;
	cout << "SUB_PROBLEM P2" << endl;

	double rhs = floor(varVals[index]) + 1;
	char sense = 'G';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	cout << "insert inequality x_" << index << " >= " << rhs << endl;
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
		cout << "delete last inequality " << endl;

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

		cout << "delete last inequality " << endl;
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg,
				&idx, &coef, 0, 0);


		cut_A.push_back(cut);
		cut_b.push_back(rhs);


		cout << "Resolve a new problem P1.. " << endl;
		cout << "add inequality x_" << index << " >= " << rhs  << endl;
		cout << "Now the new problem master is: "  << endl;
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

		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem.sol");

		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		cout << "delete last inequality " << endl;

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

		cout << "delete last inequality " << endl;
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg,
				&idx, &coef, 0, 0);
		cout << "Resolve a new problem P2.. " << endl;
		cout << "add inequality x_" << index << " <= " << rhs << endl;
		cout << "Now the new problem master is: "  << endl;

		cut_A.push_back(cut);
		cut_b.push_back(rhs);

		solve(env, lp);
	}

	return z;
}
/**
 Method solve, solve the original problem, then branch in P_1 and P_2 problem
 and resolve recursively the two sub_problems.
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variables selected
 @return void
 */
void solve(CEnv env, Prob lp) {

	// --------------------------------------------------
	// 3. solve linear problem
	// --------------------------------------------------
	cout << "PROBLEM MASTER:" << endl;
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

	//if you want to slow the iteration
	//usleep(1000000);

	// --------------------------------------------------------
	// 7. if x solution aren't integer create P1 and P2 problem
	// --------------------------------------------------------
	if (index != -1) {

		cout << endl;
		cout << "More fractional variable choose " << varVals[index] << endl;

		cout << "Index of variable choose: " << index << endl;

		//create problem P_1
		create_P1_prob(env, lp, index);

		// --------------------------------------------------------
		// 8. solve sub_problems (P_1 and P_2) return min solution
		// --------------------------------------------------------
		double* z = solve_P1_Problem(env, lp, index);

		/////////////////////////////////////////////////////


		// ------------------------------------------------
		// 9. only if both problems have solution else take
		//		the best solution and stop
		// ------------------------------------------------
		if (*z < CPX_INFBOUND && *(z + 1) < CPX_INFBOUND && flag_find) {
			flag_find = false;

			min_sol = std::min(*z, *(z + 1));
			k=index;
		}

	}else{
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		cout << "The last solution is the best integer solution. STOP" << endl;
		CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem.sol");
	}

}


/**
 Method that evaluate the vector r as sum of rows of C (matrix of dual problem (13-19) + 21)
 @param  none
 @return vector<double>
 */
std::vector<double> evaluate_rT() {

	std::vector<double> rt;
	double sum = 0;

	// --------------------------------------------------
	//evaluate A_T + e_k -1;
	// --------------------------------------------------
	for (int i = 0; i < N; i++) {

		sum = 0;
		for (int j = 0; j < num_constraint; j++)
			sum += A[j][i];

		for (std::vector<double*>::const_iterator j = cut_A.begin();
				j != cut_A.end(); ++j) {
			double*ptr = *j;
			sum += ptr[i];
		}

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

	for (std::vector<double>::const_iterator j = cut_b.begin();
			j != cut_b.end(); ++j)
		sum += *j;

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

	int iterator_b = 0;
	for (std::vector<double*>::const_iterator j = cut_A.begin();
			j != cut_A.end(); ++j) {

		sum = 0;

		for (int i = 0; i < N; i++) {

			double*ptr = *j;
			sum += ptr[i];
		}

		sum += cut_b[iterator_b];

		// u
		sum++;

		// v
		sum++;

		rt.push_back(sum);

		iterator_b++;

	}

	cout << "vettore r " << endl;
	for (std::vector<double>::const_iterator j = rt.begin(); j != rt.end(); ++j)
		cout << *j << " ";
	cout << endl;

	return rt;
}

void print_cut_A() {
	cout << "added cuts : " << endl;
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
	cout << "vector_b of added cuts : " << endl;
	for (std::vector<double>::const_iterator i = cut_b.begin();
			i != cut_b.end(); ++i) {
		std::cout << *i << ' ';
		cout << endl;
	}
}

void print_u_variables() {
	cout << "u variables (the last is u_0): " << endl;
	for (std::vector<double>::const_iterator i = dual_varVals_P1.begin();
			i != dual_varVals_P1.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}

void print_v_variables() {

	cout << "v variables (the last is v_0): " << endl;
	for (std::vector<double>::const_iterator i =
			dual_varVals_P2.begin(); i != dual_varVals_P2.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}


