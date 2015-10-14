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
std::vector<double*> cut_A;
std::vector<double> cut_b;
int gam;


void create_P1_prob(CEnv env, Prob lp, int index);
void create_P2_prob(CEnv env, Prob lp, int index);

double* solve_P1_Problem(CEnv env, Prob lp, int index);
double solve_P2_Problem(CEnv env, Prob lp, int index);

double solve(CEnv env, Prob lp);

void print_cut_A();
void print_cut_b();

#endif /* SOURCE_SOLVE_H_ */
