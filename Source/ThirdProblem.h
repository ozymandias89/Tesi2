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
	ThirdProblem(vector <double> y_til);
	virtual ~ThirdProblem();

	// PREDICATES:

	void y_tilde_MIN_y_bar();

	/**
	 Method that set the third problem
	 @param  (CEnv env, Prob lp)
	 @return void
	 */
	void setup(CEnv env, Prob lp);

	/**
	 update y_bar step 8.4
	 @param  (CEnv env, Prob lp)
	 @return void
	 */
	void update_y_bar(CEnv env, Prob lp);

	// ATTRIBUTES:
	vector<double> y_tilde;
	vector<double> t;
};

#endif /* SOURCE_THIRDPROBLEM_H_ */
