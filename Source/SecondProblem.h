/*
 * SecondProblem.h
 *
 *  Created on: 04 nov 2015
 *      Author: riccardo
 */

#ifndef SOURCE_SECONDPROBLEM_H_
#define SOURCE_SECONDPROBLEM_H_


// Includes:
#include "load.h"

using namespace std;

class SecondProblem {
public:
	SecondProblem();
	virtual ~SecondProblem();


	// PREDICATES:
	/**
	 Method that set the second problem
	 @param  (CEnv env, Prob lp)
	 @return void
	 */
	void setupSP(CEnv env, Prob lp);


	/** Method that evaluate the vector r as sum of rows of C (matrix of dual problem (13-19) + 21)
	 @param  none
	 @return vector<double>
	 */
	std::vector<double> evaluate_rT();


	void set_solution(CEnv env, Prob lp);


	/**
	 Method that solve the second problem (13-19) + 21)
	 @param  none
	 @return vector<double>
	 */
	void solve(CEnv env, Prob lp);

	// ATTRIBUTES:
	std::set <double> R;
	vector<double> u;
	vector<double> v;
	vector<double> a;
	double u0;
	double v0;
	double beta;

};

#endif /* SOURCE_SECONDPROBLEM_H_ */
