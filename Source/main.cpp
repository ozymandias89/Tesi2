/*
 @file    main.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
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
 * 				with P1 and P2 new problem to resolve
 *
 * -----------------x-----------------------------------------------------------
 */

#include "load.h"

int main(int argc, char const *argv[]) {

	if (argc < 2) { // per la ricerca locale servono 2 parametri
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
		CPXsetdblparam(env, CPX_PARAM_TILIM, 1800.0);


		// --------------------------------------------------
		// 2. write problem
		// --------------------------------------------------
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);


		// --------------------------------------------------
		// 3. Optimization
		// --------------------------------------------------
		CHECKED_CPX_CALL(CPXmipopt, env, lp);


		// --------------------------------------------------
		// 4. print solution
		// --------------------------------------------------
		print_objval(env, lp);


		// --------------------------------------------------
		// 5. set number and value of variable
		//    (cur_numcols,varVals) and print these
		// --------------------------------------------------
		print_variable(env, lp);


		// --------------------------------------------------
		// 6. chose the best fractional variable
		// --------------------------------------------------
		int index = fractionar_variable(varVals);

		// --------------------------------------------------------
		// 7. if x solution aren't integer create P1 and P2 problem
		// --------------------------------------------------------
		if (index != -1) {

			cout << endl;
			cout << "Higher fractional variable choose " << varVals[index]
					<< endl;

			cout << "Index " << index << endl;

			create_P1_Problem(env, lp, index);
			/////////////////////////////////////////////////////
			create_P2_Problem(env, lp, index);

		}


		// ---------------------------------------------------------
		// 8. free allocate memory
		// ---------------------------------------------------------
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);

		for (int i = 0; i < num_constraint; ++i) {
			delete[] A[i];
		}
		delete[] A;
		delete[] b;
		delete[] c;

	} catch (std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}

	return 0;
}
