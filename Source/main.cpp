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

#include "solve.cpp"
#include "SecondProblem.cpp"
#include "ThirdProblem.cpp"

int main(int argc, char const *argv[]) {

	if (argc < 2) { // you must insert file
		throw std::runtime_error("usage: ./main filename.txt");

	}

	ifstream myfile(argv[1], std::ios::in);

	load_problem(myfile);

	myfile.close();

	try {

		// --------------------------------------------------
		// 1. Initialization problem
		// --------------------------------------------------
		DECL_ENV(env);
		DECL_PROB(env, lp, "resolve problem RL");
		setupLP(env, lp);

		// --------------------------------------------------
		// 2. program (first part)
		// --------------------------------------------------
		do {

			solve(env, lp);

			// ---------------------------------------------------------
			// 3. if P_1 and P_2 have solution initialization of the second problem
			// ---------------------------------------------------------
			DECL_ENV(env_dual);
			DECL_PROB(env_dual, lp_dual, "resolve second problem");

			SecondProblem* sec_prob = new SecondProblem();
			sec_prob->setupSP(env_dual, lp_dual);

			// --------------------------------------------------
			// 5. Evaluate vector r
			// --------------------------------------------------
			sec_prob->evaluate_rT();

			// --------------------------------------------------
			// 6. Cycle step 8
			// --------------------------------------------------
			bool flag;
			int original_constraint = CPXgetnumrows(env_dual, lp_dual);

			do {
				sec_prob->print_c();

				sec_prob->step8_1(env_dual, lp_dual);

				sec_prob->step8_2(env_dual, lp_dual);

				int num_constraint = CPXgetnumrows(env_dual, lp_dual);

				sec_prob->solve(env_dual, lp_dual);

				// --------------------------------------------------
				// 7. STOP condition
				// --------------------------------------------------
				flag = sec_prob->y_tilde_EQ_y_bar();

				if (!flag) {
					DECL_ENV(env_third);
					DECL_PROB(env_third, lp_third, "resolve third problem");
					ThirdProblem* third_prob = new ThirdProblem(sec_prob->y_tilde, sec_prob->cost);
					third_prob->setup(env_third, lp_third);
					CHECKED_CPX_CALL(CPXwriteprob, env_dual, lp_dual, "../data/second_problem.lp", 0);
					CHECKED_CPX_CALL(CPXwriteprob, env_third, lp_third, "../data/third_problem.lp", 0);
					third_prob->solve(env_third, lp_third);
					third_prob->update_y_bar(env_third, lp_third, sec_prob->cost);

					free(third_prob);

				}

				CHECKED_CPX_CALL(CPXdelrows, env_dual, lp_dual,
						original_constraint, num_constraint - 1);

			} while (!flag);

			// --------------------------------------------------
			// 8. ADD constraint R in the first problem
			// --------------------------------------------------
			add_constraint_R(env, lp, sec_prob->R);

			CPXfreeprob(env_dual, &lp_dual);
			CPXcloseCPLEX(&env_dual);
			free(sec_prob);
			flag_find = true;

			CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);

		} while (1);

	} catch (std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}

	return 0;
}
