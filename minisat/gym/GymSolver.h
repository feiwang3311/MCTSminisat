#ifndef Minisat_GymSolver_h
#define Minisat_GymSolver_h

#include "minisat/simp/SimpSolver.h"

namespace Minisat {

class GymSolver {
	
	SimpSolver S;

public:
	GymSolver(char*);
	void step(int); 
	double getReward();
	bool getDone();
	char* getState();

};

}

#endif