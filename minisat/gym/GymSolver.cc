#include <errno.h>
#include <zlib.h>

#include "minisat/utils/System.h"
#include "minisat/utils/ParseUtils.h"
#include "minisat/utils/Options.h"
#include "minisat/core/Dimacs.h"
#include "minisat/simp/SimpSolver.h"
#include "minisat/gym/GymSolver.h"

using namespace Minisat;

//=================================================================================================
// Constructor/Destructor:

GymSolver::GymSolver(char* satProb) {
    
	IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 0, IntRange(0, 2));
	S.verbosity = verb;
	gzFile in = gzopen(satProb, "rb");
    if (in == NULL) {
    	printf("ERROR! Could not open file: %s\n", satProb); 
    	exit(1);
    }
    parse_DIMACS(in, S, true);
    gzclose(in);
    asprintf(&(S.snapTo), "%s%s", satProb, "snaps");

    S.eliminate(true);
    if (!S.okay()){
    	printf("ERROR! SAT problem from file: %s is UNSAT by simplification\n", satProb);
    	exit(1);    
    }    
}

bool GymSolver::init(float* array, int n) {
    // Comments by Fei: Now the solveLimited() function really just initialize the problem. It needs steps to finish up!
    vec<Lit> dummy;
    S.write_state_to = array;
    S.solveLimited(dummy);
    return S.env_hold; // return false if the problem is finished by simplification
}

int GymSolver::simulate(float* array, int n, float* pi, int m, float* v, int t) {
    return S.simulate(array, pi, v[0]);
}

void GymSolver::get_visit_count(float* array, int n){
    S.get_visit_count(array);
}

void GymSolver::set_decision(int decision) {
    if (decision < 0) {
        S.agent_decision = S.default_pickLit();
    } else {
        S.agent_decision = toLit(decision);
    }
}

void GymSolver::step(float* array, int n) {
    S.write_state_to = array;
    S.step();
    /*
    if (decision == 32767) {
        S.agent_decision = S.default_pickLit(); // if the decision is MaxInt number, let the minisat decide!
    } else {
    	S.agent_decision = mkLit(abs(decision)-1, decision < 0);
    }
	S.step();
    */
}

void GymSolver::step_forward(int decision) {
    set_decision(decision);
    S.step();
}

double GymSolver::getReward() {
	return S.env_reward;
}

bool GymSolver::getDone() {
	return !S.env_hold;
}

char* GymSolver::getState() {
	//return S.snapTo;
    return S.env_state;
}

