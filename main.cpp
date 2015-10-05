/*
 @file    main.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

/* --*- C++ -*------x-----------------------------------------------------------
 *
 *
 * Description:
 *
 * -----------------x-----------------------------------------------------------
 */

// Includes:
#include <cstdio>
#include <iostream>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include "cpxmacro.h"


using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

// error status and message buffer
int status;
char errmsg[BUF_SIZE];

//if you want record time
timeval start, stop;
double elapsedTime;


const int NAME_SIZE = 512;
char name[NAME_SIZE];

//number of variable
int const N = 10;

//number of constraint
int const num_constraint = 5;

//Coefficient cost
double c[N];

//price matrix
double** A;

//known terms
double b[num_constraint];

void load_problem(ifstream &myfile) {

    if (!myfile.is_open())
        cerr << "Unable to open file!!!!!";

    for(int i = 0; i < N; ++i){
    	myfile >> c[i];
        cout << "numero letto funz. ob" << c[i] << endl;
    }


    A = new double* [N];
        for(int i = 0; i < N; ++i)
            A[i] = new double[N];


    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
           myfile >> A[i][j];
            cout << "numero letto mat A" << A[i][j] << endl;
        }
    }

    for(int i = 0; i < num_constraint; ++i){
       	myfile >> b[i];
           cout << "numero letto termin noto" << c[i] << endl;
       }



}



void setupLP(CEnv env, Prob lp) {

	{

		char varType = 'I';
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

	// vincolo

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

int main(int argc, char const *argv[]) {

    if (argc < 2) { // per la ricerca locale servono 2 parametri
            throw std::runtime_error("usage: ./main filename.dat");}

	ifstream myfile (argv[1], std::ios::in);

	load_problem(myfile);

	myfile.close();

//  try {
//        // init
//
//	  	//declaration of environment
//        DECL_ENV(env);
//        //call to name the problem
//        DECL_PROB(env, lp, "resolve problem RL");
//
//        setupLP(env, lp); // setup LP
//
//        CPXsetdblparam(env, CPX_PARAM_TILIM, 1800.0); //upper bound of time for Cplex resolution in seconds, return the best solution it has so far
//
//        // print problem description to a file
//        CHECKED_CPX_CALL(CPXwriteprob, env, lp, "quadro_elettrico.lp", 0);
//
//
//        struct timeval  tv1, tv2;
//		gettimeofday(&tv1, NULL);
//        gettimeofday(&tv1, NULL);
//        clock_t start = clock();
//        CHECKED_CPX_CALL(CPXmipopt, env, lp); // optimize
//        clock_t end = clock();
//        gettimeofday(&tv2, NULL);
//
//
//
//        // stampa il valore della funzione obiettivo
//        double objval;
//        CHECKED_CPX_CALL(CPXgetobjval, env, lp, &objval);
//        std::cout << "Obj val: " << objval << std::endl;
//
//        // get the number and value of all variables
//        int cur_numcols = CPXgetnumcols(env, lp);
//        std::vector<double> varVals;
//        varVals.resize(cur_numcols);
//        CHECKED_CPX_CALL(CPXgetx, env, lp, &varVals[0], 0, cur_numcols - 1);
//
//        int surplus; // will contain the space missing to save the column names
//
//        // se a CPXgetcolname passiamo una array troppo piccolo per salvarci
//        // tutti i nomi delle variabili, lui ritorna un codice d'errore e scrive
//        // in surplus quanto spazio manca (un valore negativo)
//        // NULL arguments to obtain the total memory required
//        status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0, cur_numcols - 1);
//        int cur_colnamespace = -surplus; // the space needed to save the names
//
//        // allocate memory
//        char** cur_colname = (char **) malloc(sizeof (char *)*cur_numcols);
//        char* cur_colnamestore = (char *) malloc(cur_colnamespace);
//
//        // get the names
//        CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace, &surplus, 0, cur_numcols - 1);
//
//        // print index, name and value of each column
//        //  for (int i = 0; i < cur_numcols; i++) {
//        //    std::cout << "Column " << i << ", " << cur_colname[i] << " = " << varVals[i] << std::endl;
//        //}
//
//        // stampa la soluzione in formato standard
//        CHECKED_CPX_CALL(CPXsolwrite, env, lp, "quadro_elettrico.sol");
//
//        std::cout << "TEMPO: " << (double)(tv2.tv_sec+tv2.tv_usec*1e-6 - (tv1.tv_sec+tv1.tv_usec*1e-6)) << " seconds (user time)\n";
//        cout << " TEMPO: " << ((float)(end - start))/ CLOCKS_PER_SEC << " seconds (CPU time)\n";
//
//        // free
//        free(cur_colname);
//        free(cur_colnamestore);
//        CPXfreeprob(env, &lp);
//        CPXcloseCPLEX(&env);
//
//        gettimeofday(&stop, NULL);
//
//        for(int i = 0; i < N; ++i) {
//                delete [] A[i];
//            }
//            delete [] A;
//
//
//    } catch (std::exception& e) {
//        std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
//    }

    return 0;
}
