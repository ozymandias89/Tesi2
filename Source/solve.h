/**
 * @file solve.h
 * library for support
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 1.0
 */

#ifndef SOURCE_SOLVE_H_
#define SOURCE_SOLVE_H_

#include "load.h"
static bool flag_find = true;
static bool flag_step1 = true;


void create_P1_prob(CEnv env, Prob lp, int index, bool verbose);
void create_P2_prob(CEnv env, Prob lp, int index, bool verbose);

double* solve_P1_Problem(CEnv env, Prob lp, int index, bool verbose);
double solve_P2_Problem(CEnv env, Prob lp, int index, bool verbose);

void solve_integer_problem(CEnv env, Prob lp, bool verbose=false);


void remove_constraint(CEnv env, Prob lp, int constraint,  bool verbose);


void step1(CEnv env, Prob lp, bool verbose);

void solve(CEnv env, Prob lp, bool verbose=false);

std::vector<double> evaluate_rT();


void print_u_variables();
void print_v_variables();

#endif /* SOURCE_SOLVE_H_ */
