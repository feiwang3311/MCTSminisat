#ifndef Minisat_GymSolver_h
#define Minisat_GymSolver_h

#include "minisat/simp/SimpSolver.h"

namespace Minisat {

class GymSolver {
	
	SimpSolver S;

public:
	GymSolver(char*);                 // set up the basics for char*, which is the filename of the SAT problem
	bool   init(float* array, int n); // initialize the SAT problem and return the state. 
									  // If return false, the Solver is in finished state and array is empty
									  // one should call the constructor and the init() to reset on a SAT problem.

	// this function asks the solver to run a simulation
	// argument array is space to write the state to (once a simulation is run, normally we get a new state, whose pi and v needs evaluation)
	// argument pi and v are for the state of the last call of simulate()
	// no need to feed simulate an action or choice to simulate on, because simulation is based on MCTS rules within the solver and shadow class
	// if return = 0, or (00 in binary): no more simulation needed (meet the size constraint), no more state evaluation needed (leaf_shadow is NULL)
	// if return = 1, or (01 in binary): no more simulation needed (meet the size constraint), the leaf_shadow state still needs to be evaluated
	// if return = 2, or (10 in binary): more simulation needed (dismeet the size constraint), no more state evaluation needed (leaf shadow is NULL-probably due to finished state)
	// if return = 3, or (11 in binary): more simulation needed (dismeet the size constraint), the leaf_shadow state still needs to be evaluated
	// one should call simulate until the result is 0, to build a complete MCTS. 
	int    simulate(float* array, int n, float* pi, int m, float* v, int t); 
	void   get_visit_count(float* array, int n); // get the nn vector from the root of MCTS (for PI)

	void   set_decision(int decision);           // set the decision for the next step() call
	void   step(float* array, int n);            // the step call (real step, not simulation)
												 // one should call the set_decision() and step() to make a real step
	
	double getReward();                          // get the reward (most likely -1 for all intermediate steps)
	bool   getDone();                            // get if the state is done                     
	char*  getState();                           // get the pointer where state can be write to (NO LONGER FUNCTIONAL)
};

}

#endif