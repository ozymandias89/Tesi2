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
 * 				then print the higher fractional variable
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

		// get the number and value of all variables
		int cur_numcols = CPXgetnumcols(env, lp);
		std::vector<double> varVals;
		varVals.resize(cur_numcols);
		CHECKED_CPX_CALL(CPXgetx, env, lp, &varVals[0], 0, cur_numcols - 1);

		int surplus; // will contain the space missing to save the column names

		// se a CPXgetcolname passiamo una array troppo piccolo per salvarci
		// tutti i nomi delle variabili, lui ritorna un codice d'errore e scrive
		// in surplus quanto spazio manca (un valore negativo)
		// NULL arguments to obtain the total memory required
		status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0,
				cur_numcols - 1);
		int cur_colnamespace = -surplus; // the space needed to save the names

		// allocate memory
		char** cur_colname = (char **) malloc(sizeof(char *) * cur_numcols);
		char* cur_colnamestore = (char *) malloc(cur_colnamespace);

		// get the names
		CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace,
				&surplus, 0, cur_numcols - 1);

		// print solution in standard format
		CHECKED_CPX_CALL(CPXsolwrite, env, lp, "../data/problem.sol");

		//  print index, name and value of each column
		for (int i = 0; i < cur_numcols; i++) {
			std::cout << cur_colname[i] << " = " << varVals[i] << std::endl;
		}

		//chose the best fractional variable
		int index = fractionar_variable(varVals);

		//if x solution aren't integer
		if (index != -1) {

			cout << endl;
			cout << "Higher fractional variable choose " << varVals[index]
					<< endl;

			cout << "Index " << index << endl;
			double rhs = floor(varVals[index]);

			char sense = 'L';
			int matbeg = 0;
			const int idx = index;
			const double coef = 1;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense,
					&matbeg, &idx, &coef, 0, 0);

			// print P1 problem to a file
			CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp1", 0);

			int cur_numrows = CPXgetnumrows(env, lp);

			CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1,
					cur_numrows - 1);

//			// print solution in standard format
//			           CHECKED_CPX_CALL(CPXsolwrite, env, lp, "problem1.sol");
			/////////////////////////////////////////////////////

			rhs = ceil(varVals[index]);
			sense = 'G';

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense,
					&matbeg, &idx, &coef, 0, 0);

			CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp2", 0);

//			// print solution in standard format
//						           CHECKED_CPX_CALL(CPXsolwrite, env, lp, "problem2.sol");
		}

		// free
		free(cur_colname);
		free(cur_colnamestore);
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
