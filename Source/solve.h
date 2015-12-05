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

/**
 Method that create P1 problem (make a branch of admissible region),
 @param  (CEnv env, Prob lp, int index, bool verbose)
 Environment of the problem, problem , index of fractional variable selected and verbose
 @return none
 */
void create_P1_prob(CEnv env, Prob lp, int index, bool verbose);

/**
 Method that create P2 problem (make a branch of admissible region),
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variable selected
 @return none
 */
void create_P2_prob(CEnv env, Prob lp, int index, bool verbose);

/**
 Method that create P1 sub_problem (make a branch of admissible region),
 return solution if exist, otherwise return INFINITE value
 @param  (CEnv env, Prob lp, index, bool verbose),
 Environment of the problem, problem , index of fractional variable selected and bool verbose
 @return double*, pointer to vector[2] where the first element is
 the result of P1 sub_problem, the second is the result of sub_P2 problem
 */
double* solve_P1_Problem(CEnv env, Prob lp, int index, bool verbose);

/**
 Method that create P2 problem (make a branch of admissible region),
 return solution if exist, otherwise return INFINITE value
 @param  (CEnv env, Prob lp, int index, bool verbose)
 Environment of the problem, problem , index of fractional variable selected and verbose
 @return double , result of P2 sub_problem
 */
double solve_P2_Problem(CEnv env, Prob lp, int index, bool verbose);

/**
 Method that solve the problem to integer
 @param  (CEnv env, Prob lp, int index, bool verbose)
 Environment of the problem, problem , index of fractional variable selected and verbose
 @return void
 */
void solve_integer_problem(CEnv env, Prob lp, bool verbose = false);

/**
 Ausiliar method that remove constraint and update data structures
 @param  (CEnv env, Prob lp, bool verbose),
 Environment of the problem, problem and verbose
 @return void
 */
void remove_constraint(CEnv env, Prob lp, int constraint, bool verbose);

/**
 Method that solve MIP problem, then remove redundant constraints from primary problem
 @param  (CEnv env, Prob lp, bool verbose),
 Environment of the problem, problem and verbose
 @return void
 */
void step1(CEnv env, Prob lp, bool verbose);

/**
 Method solve, solve the original problem, then branch in P_1 and P_2 problem
 and resolve recursively the two sub_problems.
 @param  (CEnv env, Prob lp, index),
 Environment of the problem, problem and index of fractional variables selected
 @return void
 */
void solve(CEnv env, Prob lp, bool verbose = false);

/**
 Print u variable primal problem
 @param none
 @return void
 */
void print_u_variables();

/**
 Print v variable primal problem
 @param none
 @return void
 */
void print_v_variables();

#endif /* SOURCE_SOLVE_H_ */
