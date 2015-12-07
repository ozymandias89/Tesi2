/*
 * Thirdproblem.h
 *
 *  Created on: 14 nov 2015
 *      Author: riccardo
 */

#ifndef SOURCE_THIRDPROBLEM_H_
#define SOURCE_THIRDPROBLEM_H_

// Includes:
#include "load.h"

class ThirdProblem {
public:
	ThirdProblem(vector<double> y_til, vector<double> cost,
			bool verbose = false);
	virtual ~ThirdProblem();

	// PREDICATES:

	/**
	 Print vector
	 \param  (vector<double> vector)
	 \return void
	 */
	void print_vector(vector<double> vector);

	/**
	 Calculate y bar minimum y tilde
	 \param  (vector<double> vector)
	 \return void
	 */
	void y_bar_MIN_y_tilde(vector<double> c);

	/**
	 Method that set the third problem
	 \param  (CEnv env, Prob lp)
	 \return void
	 */
	void setup(CEnv env, Prob lp);

	/**
	 solve third problem and calculate lambda
	 \param  (CEnv env, Prob lp, vector<double>& c)
	 \return void
	 */
	void solve(CEnv env, Prob lp);

	/**
	 update y_bar step 8.4
	 \param  (CEnv env, Prob lp, vector<double>& c)
	 \return void
	 */
	void update_y_bar(CEnv env, Prob lp, vector<double>& c);

	// ATTRIBUTES:
	vector<double> y_tilde;
	vector<double> t;
	double lambda;
	bool verbose;
};

#endif /* SOURCE_THIRDPROBLEM_H_ */
