/**
 * @file solve.cpp
 * solve file
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 2.0
 */

#include "solve.h"
#include <unistd.h>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

void create_P1_prob(CEnv env, Prob lp, int index, bool verbose) {

	if (verbose) {
		cout << endl;
		cout << "SUB_PROBLEM P1" << endl;
	}

	double rhs = floor(varVals[index]);

	char sense = 'L';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	if (verbose)
		cout << "insert inequality x_" << index << " <= " << rhs << endl;

}

void create_P2_prob(CEnv env, Prob lp, int index, bool verbose) {

	if (verbose) {
		cout << endl;
		cout << "SUB_PROBLEM P2" << endl;
	}

	double rhs = floor(varVals[index]) + 1;
	char sense = 'G';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	if (verbose)
		cout << "insert inequality x_" << index << " >= " << rhs << endl;

}

double* solve_P1_Problem(CEnv env, Prob lp, int index, bool verbose) {

	static double z[2];
	z[0] = CPX_INFBOUND;
	z[1] = CPX_INFBOUND;
	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	int stat = CPXgetstat(env, lp);
	if (verbose)
		cout << endl << "Status problem " << stat << endl;

	int cur_numrows = CPXgetnumrows(env, lp);

// print and set solution and create and resolve P_2 problem"
	if (stat != CPX_STAT_INFEASIBLE ) {
		if (verbose)
			cout << "FEASIBLE " << endl;
		gam = floor(varVals[index]);
		print_objval(env, lp, verbose);
		set_and_print_var_P(env, lp, verbose);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z[0]);
		set_and_print_var_D(env, lp, true, verbose);
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		if (verbose)
			cout << "delete last inequality " << endl;

		create_P2_prob(env, lp, index, verbose);
		z[1] = solve_P2_Problem(env, lp, index, verbose);

	} else {
		if (verbose)
			cout << "No solution for P1 problem exists.. " << endl;

		// add Slack variables
		static const char* varType = NULL;
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;
		snprintf(name, NAME_SIZE, "S_%i", slack);
		slack++;
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		//add cut
		std::vector<int> idx;
		std::vector<double> coef;

		double rhs = floor(varVals[index]) + 1;
		char sense = 'E';
		int matbeg = 0;

		// x_k
		idx.push_back(index);
		coef.push_back(1);

		// S
		idx.push_back(N);
		coef.push_back(-1);

		num_constraint++;
		N++;

		//add 0 to c
		c.push_back(0);

		//extend A matrix
		A.resize(num_constraint);
		for (int i = 0; i < num_constraint; i++)
			A[i].resize(N);

		for (int i = 0; i < N; i++) {
			if (i == index) {
				A[(num_constraint - 1)][i] = 1;
			} else if (i == N - 1) {
				A[(num_constraint - 1)][i] = -1;
			} else
				A[(num_constraint - 1)][i] = 0;
		}

		if (verbose)
			cout << "delete last inequality " << endl;
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 2, &rhs, &sense, &matbeg,
				&idx[0], &coef[0], 0, 0);

		b.push_back(rhs);

		if (verbose) {
			cout << "Resolve a new problem P1.. " << endl;
			cout << "add inequality x_" << index << " >= " << rhs << endl;
			cout << "Restart from step 1 with new problem master: " << endl;
		}

		flag_step1 = true;
		solve(env, lp, verbose);

	}

	return z;

}

double solve_P2_Problem(CEnv env, Prob lp, int index, bool verbose) {

	double z = CPX_INFBOUND;

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	//bool infeasible = test_problem_infeasible(env, lp, verbose);

	int stat = CPXgetstat(env, lp);
		if (verbose)
			cout << endl << "Status problem " << stat << endl;

	int cur_numrows = CPXgetnumrows(env, lp);

	if (stat != CPX_STAT_INFEASIBLE ) {
		if (verbose)
			cout << "FEASIBLE " << endl;
		print_objval(env, lp, verbose);
		set_and_print_var_P(env, lp, verbose);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z);
		set_and_print_var_D(env, lp, false, verbose);
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		if (verbose)
			cout << "delete last inequality " << endl;

	} else {
		if (verbose)
			cout << "No solution for P2 problem exists.. " << endl;

		// add Slack variables
		static const char* varType = NULL;
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;
		snprintf(name, NAME_SIZE, "S_%i", slack);
		slack++;
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		//add cut
		std::vector<int> idx;
		std::vector<double> coef;

		double rhs = floor(varVals[index]);
		char sense = 'E';
		int matbeg = 0;

		//x_k
		idx.push_back(index);
		coef.push_back(1);

		//S
		idx.push_back(N);
		coef.push_back(1);

		num_constraint++;
		N++;

		//add 0 to c
		c.push_back(0);

		//extend A matrix
		A.resize(num_constraint);
		for (int i = 0; i < num_constraint; i++)
			A[i].resize(N);

		for (int i = 0; i < N; i++) {
			if (i == index) {
				A[(num_constraint - 1)][i] = 1;
			} else if (i == N - 1) {
				A[(num_constraint - 1)][i] = 1;
			} else
				A[(num_constraint - 1)][i] = 0;
		}

		if (verbose) {
			for (int i = 0; i < N; i++)
				cout << A[(num_constraint - 1)][i] << " ";

			cout << endl;

			cout << "delete last inequality " << endl;
		}
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 2, &rhs, &sense, &matbeg,
				&idx[0], &coef[0], 0, 0);

		if (verbose) {
			cout << "Resolve a new problem P2.. " << endl;
			cout << "add inequality x_" << index << " <= " << rhs << endl;
			cout << "Restart from step 1 with new problem master: " << endl;
		}

		b.push_back(rhs);
		flag_step1 = true;
		solve(env, lp, verbose);
	}

	return z;
}

void solve_integer_problem(CEnv env, Prob lp, bool verbose) {

	//change problem to integer
	CHECKED_CPX_CALL(CPXchgprobtype, env, lp, CPXPROB_MILP);

	//set varibles to integer
	char ctype[N];
	for (int i = 0; i < N; i++) {
		if (i < Num_original_variables) {
			ctype[i] = CPX_INTEGER;
		} else
			ctype[i] = CPX_CONTINUOUS;
	}
	//const_cast<char *>(ctype);
	CHECKED_CPX_CALL(CPXcopyctype, env, lp, ctype);

	//check feasible
	//bool infeasible = test_problem_infeasible(env, lp, verbose);

	CHECKED_CPX_CALL(CPXmipopt, env, lp);
	int stat = CPXgetstat(env, lp);

	if (stat != CPXMIP_INFEASIBLE) {
		if (verbose) {
			cout << "Problem solved to integer: " << endl;
			print_objval(env, lp, verbose);
			double risult;
			CHECKED_CPX_CALL(CPXgetobjval, env, lp, &risult);
			if (integer == -CPX_INFBOUND)
				integer = risult;
			else if (fabs(integer - risult) > 10e+06 ) {
				cerr << "We have lose one integer solution: " << endl;
				cout << "Number of constraints added: " << CPXgetnumrows(env, lp) - Num_original_constraints << endl;
				exit(1);
			}
			set_and_print_var_P(env, lp, verbose);
		}

	} else {
		cout << endl;
		cout << " Integer problem not resolvable! " << endl;
		cout << " Iteration number: " << iter << endl;
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		// free allocate memory
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(const_cast<cpxenv **>(&env));
		exit(0);
	}

	CHECKED_CPX_CALL(CPXchgprobtype, env, lp, CPXPROB_LP);
}

void remove_constraint(CEnv env, Prob lp, int constraint, bool verbose) {

	if (verbose) {
		print_matrix();
		print_vect_c();
		print_vect_b();
		cout << "Delete redundant constraint number: " << constraint << endl;
	}

	CHECKED_CPX_CALL(CPXdelrows, env, lp, constraint, constraint);

	if (constraint < Num_original_constraints) {
		b.erase(b.begin() + constraint);
		A[constraint].clear();
		A.erase(A.begin() + constraint);

		Num_original_constraints--;
		num_constraint--;

	} else {

		int column = (Num_original_variables + constraint
				- Num_original_constraints);
		CHECKED_CPX_CALL(CPXdelcols, env, lp, column, column);

		b.erase(b.begin() + constraint);

		for (unsigned int i = 0; i < A.size(); ++i) {
			A[i].erase(A[i].begin() + column);
		}

		A[constraint].clear();
		A.erase(A.begin() + constraint);

		c.erase(c.begin() + column);

		num_constraint--;
		N--;
	}

	if (verbose) {
		print_matrix();
		print_vect_c();
		print_vect_b();
	}

}

void step1(CEnv env, Prob lp, bool verbose) {

	if (verbose) {
		cout << endl;
		cout << "STEP 1:" << endl;
	}

	CHECKED_CPX_CALL(CPXlpopt, env, lp);
	int stat = CPXgetstat(env, lp);

	if (verbose)
		cout << endl << "Status problem " << stat << endl;

	//check feasible
	//----------------------------------------------------------------
	// Remove redundant constraint (if exist) else STOP CONDITION 1
	//----------------------------------------------------------------
	if (stat != CPX_STAT_INFEASIBLE) {
		bool flag_redundant;
		do {
			flag_redundant = false;
			int num_cols = CPXgetnumcols(env, lp);
			int num_rows = CPXgetnumrows(env, lp);
			double redlb[num_cols], redub[num_cols];
			for (int i = 0; i < num_cols; ++i) {
				redlb[i] = 0;
				redub[i] = CPX_INFBOUND;
			}
			int rstat[num_rows];
			CHECKED_CPX_CALL(CPXbasicpresolve, env, lp, redlb, redub,
					&rstat[0]);

			int i = 0;
			while (i < num_rows) {
				if (rstat[i] == -1) {
					flag_redundant = true;
					break;
				}
				i++;
			}

			//remove constraint
			if (flag_redundant) {
				if (verbose)
					cout << "Constraint index " << i << " is redundant "
							<< endl;
				remove_constraint(env, lp, i, verbose);
			} else if (!flag_redundant && verbose) {
				cout << "No detect redundant constraints. " << endl;
				cout << endl;
			}

		} while (flag_redundant);

	} else {
		cout << endl;
		cout << " STOP CONDITION STEP 1 " << endl;
		cout << " Iteration number: " << iter << endl;
		cout << "Number of constraints added: " << CPXgetnumrows(env, lp) - Num_original_constraints;
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		// free allocate memory
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(const_cast<cpxenv **>(&env));
		exit(0);
	}

}

void solve(CEnv env, Prob lp, bool verbose) {

	// --------------------------------------------------
	// 1. STEP_1
	// --------------------------------------------------
	if (flag_step1)
		step1(env, lp, verbose);

	// --------------------------------------------------
	// 2. solve linear problem
	// --------------------------------------------------
	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	int stat = CPXgetstat(env, lp);

	// --------------------------------------------------
	// 3. STOP CONDITION
	// --------------------------------------------------
	if (stat==CPX_STAT_UNBOUNDED) {
		cout << endl;
		cout << " STOP CONDITION STEP 3 " << endl;
		cout << " Iteration number: " << iter << endl;
		cout << "Number of constraints added: " << CPXgetnumrows(env, lp) - Num_original_constraints << endl;
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		// free allocate memory
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(const_cast<cpxenv **>(&env));
		exit(0);
	}

	cout << endl << "PROBLEM MASTER:" << endl;

	// --------------------------------------------------
	// 4. print solution
	// --------------------------------------------------
	print_objval(env, lp, true);

	// --------------------------------------------------
	// 5. set number and value of variable
	//    (cur_numcols,varVals) and print these
	// --------------------------------------------------
	set_and_print_var_P(env, lp, true);

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

		if (verbose) {
			cout << endl;
			cout << "More fractional variable choose " << varVals[index]
					<< endl;

			cout << "Index of variable choose: " << index << endl;
		}

		//create problem P_1
		create_P1_prob(env, lp, index, verbose);

		// --------------------------------------------------------
		// 8. solve sub_problems (P_1 and P_2) return min solution
		// --------------------------------------------------------
		double* z = solve_P1_Problem(env, lp, index, verbose);

		/////////////////////////////////////////////////////

		// ------------------------------------------------
		// 9. only if both problems have solution else get
		//		the best solution and stop
		// ------------------------------------------------
		if (*z < CPX_INFBOUND && *(z + 1) < CPX_INFBOUND && flag_find) {
			flag_find = false;

			min_sol = std::min(*z, *(z + 1));
			k = index;
		}

	} else {
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		cout
				<< "The last solution is the best integer solution. STOP CONDITION STEP 4 "
				<< endl;
		cout << " Iteration number: " << iter << endl;
		cout << "Number of constraints added: " << CPXgetnumrows(env, lp) - Num_original_constraints << endl;
		CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem.sol");
		// free allocate memory
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(const_cast<cpxenv **>(&env));
		exit(0);
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
	for (std::vector<double>::const_iterator i = dual_varVals_P2.begin();
			i != dual_varVals_P2.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}
