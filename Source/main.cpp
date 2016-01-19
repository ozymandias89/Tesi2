/*
 @file    main.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 2.0
 */

/* --*- C++ -*------x-----------------------------------------------------------
 *
 *
 * Description: This main resolve a linear relaxation of the problem in the form
 * 				min cT x
 * 				s.t Ax=b
 * 				x€R+
 *
 * 				then print the higher fractional variable and create a branch
 * 				with P1 and P2 new problem to resolve.
 * 				In the end create problem step 7.
 *
 * -----------------x-----------------------------------------------------------
 */

#include "generator.h"
#include "solve.cpp"
#include <sstream>
#include <time.h>
#include "SecondProblem.cpp"
#include "ThirdProblem.cpp"

int main(int argc, char const *argv[]) {

	// --------------------------------------------------
	// 1. Parameters (generate problem)
	// --------------------------------------------------
	string file;
	int iteration;
	bool verbose;

	if (argc == 4) {
		file = argv[1];
		iteration = strtol(argv[2], NULL, 10);
		std::stringstream ss(argv[3]);
		if (!(ss >> std::boolalpha >> verbose)) {
			cerr << "Parameters error!!!!!" << endl;
			exit(1);
		}
	} else if (argc == 5) {
		int num_rows = strtol(argv[1], NULL, 10);
		int limit = strtol(argv[2], NULL, 10);
		iteration = strtol(argv[3], NULL, 10);
		file = generate_canonical_matrix(num_rows, limit);
		std::stringstream ss(argv[4]);
		if (!(ss >> std::boolalpha >> verbose)) {
			cerr << "Parameters error!!!!!" << endl;
			exit(1);
		}
	} else {
		cerr << "Parameters error!!!!!" << endl;
		cout << "The correct syntax is:" << endl
				<< "./main (string)'path_name_file' (int)max_iteration (bool)debug"
				<< endl << "or (to generate random matrix)" << endl
				<< "./main (int)number_constraints (int)limit_number_coefficients generate (int)max_iteration (bool)debug"
				<< endl;
		exit(1);
	}

	ifstream myfile(file.c_str(), std::ios::in);
	load_problem(myfile);

	myfile.close();

	// --------------------------------------------------
	// 2. Initialization problem
	// --------------------------------------------------

	iter = 0;

	//clock start
	clock_t t1;
	t1 = clock();

	// --------------------------------------------------
	// 1. Initialization problems
	// --------------------------------------------------
	DECL_ENV(env);
	DECL_PROB(env, lp, "resolve problem RL");

	try {

		setupLP(env, lp);

		CPXsetdblparam(env, CPXPARAM_Simplex_Tolerances_Feasibility, 1e-5);
		CPXsetdblparam(env, CPXPARAM_Simplex_Tolerances_Optimality, 1e-5);

		// --------------------------------------------------
		// 2. program (first part)
		// --------------------------------------------------

		do {
			//solve_integer_problem(env, lp, true);

			solve(env, lp, verbose);

			// --------------------------------------------------------------------
			// 3. if P_1 and P_2 have solution initialization of the second problem
			// --------------------------------------------------------------------
			DECL_ENV(env_dual);
			DECL_PROB(env_dual, lp_dual, "resolve second problem");

			SecondProblem* sec_prob = new SecondProblem(verbose);
			sec_prob->setupSP(env_dual, lp_dual);

			// --------------------------------------------------
			// 5. Evaluate vector r
			// --------------------------------------------------
			sec_prob->evaluate_rT();

			if (verbose) {
				print_vect_c();
				cout << "min sol:" << endl << min_sol << endl;
				print_u_variables();
				print_v_variables();
				cout << endl;
			}
			// --------------------------------------------------
			// 6. Cycle step 8
			// --------------------------------------------------
			bool flag;
			bool third_infeasiable;
			bool second_infeasible;

			CPXsetdblparam(env_dual, CPXPARAM_Simplex_Tolerances_Feasibility,
					1e-5);
			CPXsetdblparam(env_dual, CPXPARAM_Simplex_Tolerances_Optimality,
					1e-5);

			sec_prob->step8_1(env_dual, lp_dual);

			do {

				if (verbose)
					cout << endl
							<< "Number of constraints inserted in the cicle 8: "
							<< sec_prob->satisfy_constraint_list.size() << endl;

				clock_t t2;
				t2 = clock();
				double elapsed_secs = double(t2 - t1) / CLOCKS_PER_SEC;

				if (elapsed_secs > 120.0) {
					throw std::runtime_error("Timeout!");
				}

				sec_prob->step8_2(env_dual, lp_dual);

				CHECKED_CPX_CALL(CPXwriteprob, env_dual, lp_dual,
						"../data/second_problem.lp", 0);

				second_infeasible = sec_prob->solve(env_dual, lp_dual);

				if (!second_infeasible) {
					// --------------------------------------------------
					// 7. STOP condition
					// --------------------------------------------------
					flag = sec_prob->y_tilde_EQ_y_bar();

					if (!flag) {
						ThirdProblem* third_prob = new ThirdProblem(
								sec_prob->y_tilde, sec_prob->cost, verbose);
						if (verbose)
							print_vect_b();

						third_prob->setup();

						third_infeasiable = third_prob->solve(
								sec_prob->satisfy_constraint_list);

						if (third_prob->constraint_to_add != -1) {
							sec_prob->add_constraint(env_dual, lp_dual,
									third_prob->constraint_to_add);
						} else if (third_prob->constraint_to_add == -1
								&& !(third_infeasiable)) {
							third_infeasiable = true;
						}

						CHECKED_CPX_CALL(CPXwriteprob, env_dual, lp_dual,
								"../data/second_problem.lp", 0);

						if (!third_infeasiable) {
							third_prob->update_y_bar(sec_prob->cost);

							//delete last constraint ry=ry
							int num_constraint = CPXgetnumrows(env_dual,
									lp_dual);
							CHECKED_CPX_CALL(CPXdelrows, env_dual, lp_dual,
									num_constraint - 2, num_constraint - 2);
						}

						free(third_prob);

					}
				}
			} while (!flag && !(third_infeasiable) && !(second_infeasible));

			// --------------------------------------------------
			// 8. ADD constraint R in the first problem
			// --------------------------------------------------

			add_constraint_R(env, lp, sec_prob->R);

			CPXfreeprob(env_dual, &lp_dual);
			CPXcloseCPLEX(&env_dual);
			free(sec_prob);
			flag_find = true;

			CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);

			iter++;
			cout << "Number of iteration: " << iter << endl;
			clock_t t2;
			t2 = clock();
			double elapsed_secs = double(t2 - t1) / CLOCKS_PER_SEC;
			cout << "Elapsed time: " << elapsed_secs << endl;

		} while (iteration - iter != 0);

		cout << "Number of constraints added: "
				<< CPXgetnumrows(env, lp) - Num_original_constraints << endl;
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);

	} catch (std::exception& e) {
		cout << "Number of constraints added: "
				<< CPXgetnumrows(env, lp) - Num_original_constraints << endl;
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);

		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);

		clock_t t2;
		t2 = clock();
		double elapsed_secs = double(t2 - t1) / CLOCKS_PER_SEC;
		cout << "Elapsed time: " << elapsed_secs << endl;

	}

	return 0;
}
