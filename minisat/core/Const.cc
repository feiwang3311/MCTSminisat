#include "minisat/core/Const.h"

const float Hyper_Const::c_act = 2.8f; // TODO: need a better value here for exploration

const gsl_rng* Hyper_Const::r = gsl_rng_alloc(gsl_rng_mt19937);

const double Hyper_Const::alpha[Hyper_Const::nact] = {2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0};

void Hyper_Const::generate_dirichlet(double* di) {
    gsl_ran_dirichlet(Hyper_Const::r, Hyper_Const::nact, Hyper_Const::alpha, di);
} 
