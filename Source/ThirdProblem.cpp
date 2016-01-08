/*
 * Thirdproblem.cpp
 *
 *  Created on: 14 nov 2015
 *      Author: riccardo
 */

#include "ThirdProblem.h"

ThirdProblem::ThirdProblem(vector<double> y_til, vector<double> c,
		bool verbose) {
	this->verbose = verbose;
	this->y_tilde = y_til;
	y_bar_MIN_y_tilde(c);

}

ThirdProblem::~ThirdProblem() {
	// TODO Auto-generated destructor stub
}

void ThirdProblem::y_bar_MIN_y_tilde(vector<double> c) {

	double difference;

//	if (verbose){
//		cout << endl;
//		cout << "y_bar - y_til" << endl << endl;
//	}
	// --------------------------------------------------
	// 1. c-a
	// --------------------------------------------------
	int j = 0;
	for (unsigned int i = 0; i < c.size(); i++) {
		difference = c[i] - y_tilde[j];
		//tolerance error
		if (difference < epsilon_8_4 && difference > -epsilon_8_4)
			difference = 0.0;

		t.push_back(difference);
//		if (verbose)
//			cout << "a_" << i << "  =" << difference << endl;
		j++;
	}

	// --------------------------------------------------
	// 2. z-b
	// --------------------------------------------------

	difference = min_sol - y_tilde[j];
	//tolerance error
	if (difference < epsilon_8_4 && difference > -epsilon_8_4)
		difference = 0.0;

	t.push_back(difference);
//	if (verbose)
//		cout << "beta  = " << difference << endl;

	j++;

	// --------------------------------------------------
	// 3. u-u && u_0-u_0
	// --------------------------------------------------

	for (unsigned int i = 0; i < dual_varVals_P1.size(); i++) {
		difference = dual_varVals_P1[i] - y_tilde[j];
		//tolerance error
		if (difference < epsilon_8_4 && difference > -epsilon_8_4)
			difference = 0.0;
//		if (verbose)
//			cout << "u  =" << difference << endl;
		t.push_back(difference);
		j++;
	}

	// --------------------------------------------------
	// 5. v-v && v_0 - v_0
	// --------------------------------------------------

	for (unsigned int i = 0; i < dual_varVals_P2.size(); i++) {
		difference = dual_varVals_P2[i] - y_tilde[j];
		//tolerance error
		if (difference < epsilon_8_4 && difference > -epsilon_8_4)
			difference = 0.0;
//		if (verbose)
//			cout << "v  =" << difference << endl;
		t.push_back(difference);
		j++;
	}

	if (verbose) {
		cout << endl;
		cout << "Vector t (y_bar - y_tilde): " << endl;
		print_vector(t);
	}

}

template<typename T>
void ThirdProblem::print_vector(vector<T> vector) {

	for (unsigned int i = 0; i < vector.size(); i++)
		cout << vector[i] << " ";

	cout << endl;

}

void ThirdProblem::setup() {

	double cof, rhs;
	{
		//first constraint (second problem)
		for (int i = 0; i < N; i++) {
			cof = 0;
			rhs = 0;

			//1* y_tilde[a] && 1*t[a]
			cof += t[i];
			rhs += y_tilde[i];

			int j = N;
			//beta are 0... skip
			j++;

			//A_T * y_tilde[u] && A_T * t[u]
			for (int iter = 0; iter < num_constraint; iter++) {
				cof += A[iter][i] * t[j];
				rhs += A[iter][i] * y_tilde[j];
				j++;
			}

			if (i == k) {
				cof -= t[j];
				rhs -= y_tilde[j];
			}
			//v and v_0 are 0.. skip

			rhs = -rhs;

			//tolerance error
			if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
				rhs = 0.0;

			if (cof < epsilon_8_4 && cof > -epsilon_8_4)
				cof = 0.0;

			if (cof >= 0) {
				sense.push_back('g');
			} else
				sense.push_back('l');

			double ris;
			//add constraints
			if (cof != 0) {
				ris = rhs / cof;
			} else {
				ris = -CPX_INFBOUND;
			}

			result.push_back(ris);

		}
	}
	{	//Second constraint
		cof = 0;
		rhs = 0;

		//a is 0...skip
		int j = N;

		//-beta
		cof -= t[j];
		rhs -= y_tilde[j];
		j++;

		//u

		for (int i = 0; i < num_constraint; i++) {
			cof += b[i] * t[j];
			rhs += b[i] * y_tilde[j];
			j++;
		}

		//u_0
		cof += gam * t[j];
		rhs += gam * y_tilde[j];

		//v and v_0 are 0... skip
		rhs = -rhs;

//			//tolerance error
		if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
			rhs = 0.0;

		if (cof < epsilon_8_4 && cof > -epsilon_8_4)
			cof = 0.0;

		if (cof >= 0) {
			sense.push_back('g');
		} else
			sense.push_back('l');

		double ris;
		//add constraints
		if (cof != 0) {
			ris = rhs / cof;
		} else {
			ris = -CPX_INFBOUND;
		}

		result.push_back(ris);

	}

	{
		//for each lines of matrix (third constraint)
		for (int i = 0; i < N; i++) {

			cof = 0;
			rhs = 0;
			//1* y_tilde[a] && 1*t[a]

			cof += t[i];
			rhs += y_tilde[i];

			int j = N;
			//beta are 0... skip
			j++;

			//u && u_0 are 0... skip
			j += num_constraint + 1;

			//A_T * y_tilde[v] && A_T * t[v]
			for (int iter = 0; iter < num_constraint; iter++) {
				cof += A[iter][i] * t[j];
				rhs += A[iter][i] * y_tilde[j];
				j++;
			}

			//v_0
			if (i == k) {
				cof -= t[j];
				rhs -= y_tilde[j];
			}

			rhs = -rhs;

//				//tolerance error
			if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
				rhs = 0.0;

			if (cof < epsilon_8_4 && cof > -epsilon_8_4)
				cof = 0.0;

			if (cof >= 0) {
				sense.push_back('g');
			} else
				sense.push_back('l');

			double ris;
			//add constraints
			if (cof != 0) {
				ris = rhs / cof;
			} else {
				ris = -CPX_INFBOUND;
			}

			result.push_back(ris);

		}

	}

	{
		//fourth constraint
		cof = 0;
		rhs = 0;

		//a is 0...skip
		int j = N;

		//-beta
		cof -= t[j];
		rhs -= y_tilde[j];
		j++;

		//u and u_0 are 0... skip
		j += num_constraint + 1;

		//b_t * t[v]
		for (int i = 0; i < num_constraint; i++) {
			cof += b[i] * t[j];
			rhs += b[i] * y_tilde[j];
			j++;
		}

		//v_0
		cof += (gam + 1) * t[j];
		rhs += (gam + 1) * y_tilde[j];
		rhs = -rhs;

//			//tolerance error
		if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
			rhs = 0.0;

		if (cof < epsilon_8_4 && cof > -epsilon_8_4)
			cof = 0.0;

		if (cof >= 0) {
			sense.push_back('g');
		} else
			sense.push_back('l');

		double ris;
		//add constraints
		if (cof != 0) {
			ris = rhs / cof;
		} else {
			ris = -CPX_INFBOUND;
		}

		result.push_back(ris);

	}

	{
		// constraint -u_0 >= 0
		cof = 0;
		rhs = 0;

		int j = N + 1 + num_constraint;
		//-u
		cof -= t[j];
		rhs -= y_tilde[j];

		rhs = -rhs;
		//tolerance error
		if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
			rhs = 0.0;

		if (cof < epsilon_8_4 && cof > -epsilon_8_4)
			cof = 0.0;

		if (cof >= 0) {
			sense.push_back('g');
		} else
			sense.push_back('l');

		double ris;
		//add constraints
		if (cof != 0) {
			ris = rhs / cof;
		} else {
			ris = -CPX_INFBOUND;
		}

		result.push_back(ris);

	}

	{
		// constraint v_0 >= 0
		cof = 0;
		rhs = 0;

		//v
		cof += t.back();
		rhs += y_tilde.back();

		rhs = -rhs;
		//tolerance error
		if (rhs < epsilon_8_4 && rhs > -epsilon_8_4)
			rhs = 0.0;

		if (cof < epsilon_8_4 && cof > -epsilon_8_4)
			cof = 0.0;

		if (cof >= 0) {
			sense.push_back('g');
		} else
			sense.push_back('l');

		double ris;
		//add constraints
		if (cof != 0) {
			ris = rhs / cof;
		} else {
			ris = -CPX_INFBOUND;
		}

		result.push_back(ris);

	}
}

bool ThirdProblem::solve(set<int> constraints) {

	bool infeasible = false;
	ub = CPX_INFBOUND;
	lb = -CPX_INFBOUND;
	constraint_to_add = -1;

	for (unsigned int i = 0; i < result.size(); ++i) {

		if (sense[i] == 'g') {
			if ((result[i] - lb) > 1.e-6L) {
				lb = result[i];
			}
		} else {
			if ((result[i] - ub) < 1.e-6L) {
				ub = result[i];
			}
		}

	}

	//test if third problem is infeasible
	if ((ub - lb) < -1.e-6L) {
		if (verbose)
			cout << "Problem third is infeasible! " << endl;
		infeasible = true;
	}

	if (!infeasible) {
		for (unsigned int i = 0; i < result.size(); ++i) {
			if ((result[i] - ub) < 1.e-6L && sense[i] == 'l') {
				if (constraints.count(i) != 0)
					continue;
				else {
					constraint_to_add = i;
					break;
				}

			}
		}

		if (verbose) {
			cout << endl;
			cout << "lower bound " << lb << endl;
			cout << "upper bound " << ub << endl;
			cout << "index constraint to add: " << constraint_to_add << endl;
		}

		lambda = ub;

	}

	return infeasible;
}

void ThirdProblem::update_y_bar(vector<double>& c) {

	vector<double> r_mul_lamb;

	//-------------------------------------------------
	// lambda * (y_bar-y_tilde)
	//-------------------------------------------------

	for (vector<double>::iterator it = t.begin(); it != t.end(); ++it) {
		r_mul_lamb.push_back((*it * lambda));
	}

	//-------------------------------------------------
	// UPDATE y bar
	//-------------------------------------------------

	int j = 0;

	//c
	for (unsigned int i = 0; i < c.size(); ++i) {
		c[j] = y_tilde[j] + r_mul_lamb[j];
		j++;
	}

	//z
	min_sol = y_tilde[j] + r_mul_lamb[j];
	j++;

	//u and u_0
	for (unsigned int i = 0; i < dual_varVals_P1.size(); ++i) {
		dual_varVals_P1[i] = y_tilde[j] + r_mul_lamb[j];
		j++;
	}

	//v and v_0
	for (unsigned int i = 0; i < dual_varVals_P2.size(); ++i) {
		dual_varVals_P2[i] = y_tilde[j] + r_mul_lamb[j];
		j++;
	}

	if (verbose) {
		cout << endl << "Update y bar:" << endl;
		cout << "Vector c:" << endl;
		for (unsigned int i = 0; i < c.size(); ++i)
			cout << c[i] << " ";
		cout << "beta " << min_sol << endl;
		cout << "u variables (the last is u_0): " << endl;
		print_vector(dual_varVals_P1);
		cout << "v variables (the last is v_0): " << endl;
		print_vector(dual_varVals_P2);
		cout << endl;
	}

}
