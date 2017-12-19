#ifndef HYPER_CONST
#define HYPER_CONST

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>


class Hyper_Const
{
public :
	static const int dim0 = 100;           // max_clause
    static const int dim1 = 20;            // max_var
    static const int dim2 = 2;             // nc
    static const int nact = 40;            // nact 
    static const float c_act;     	       // c_act is a hyperparameter for MCTS (decide the level of exploration) 
    static const int MCTS_size_lim = 400; // the size of MCT we want to achieve.

    static const gsl_rng *r;              // random generater
    static const double alpha[nact];      // alpha parameter
    static void generate_dirichlet(double*);     // function used to generate dirichlet noise
};

#endif