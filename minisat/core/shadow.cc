#include <math.h>
#include "minisat/mtl/Alg.h"
#include "minisat/mtl/Sort.h"
#include "minisat/utils/System.h"
#include "minisat/core/Solver.h"
#include "minisat/core/shadow.h"

using namespace Minisat;

shadow::shadow(Solver* from) : 
	verbosity                     (from -> verbosity),
    ccmin_mode                    (from -> ccmin_mode),
    phase_saving                  (from -> phase_saving),
    learntsize_factor             (from -> learntsize_factor),
    learntsize_inc                (from -> learntsize_inc),
    garbage_frac                  (from -> garbage_frac),
    max_learnts                   (from -> max_learnts),
    learntsize_adjust_start_confl (from -> learntsize_adjust_start_confl),
    learntsize_adjust_inc         (from -> learntsize_adjust_inc),
    learntsize_adjust_confl       (from -> learntsize_adjust_confl),
    learntsize_adjust_cnt         (from -> learntsize_adjust_cnt),
    cla_inc                       (from -> cla_inc),
    clause_decay                  (from -> clause_decay),
    trail_size                    (from -> trail.size()),
    trail_lim_size				  (from -> trail_lim.size()),
    qhead                         (from -> qhead),
    ca_size						  (from -> ca.size()),
    learnts_size 				  (from -> learnts.size())
	{ 
		ca_shadow.extra_clause_field = from -> ca.extra_clause_field; 
		learnts_copy_is_uninitialized = true;
		origin = from;
		parent = NULL;
		index_child_last_pick = -1;
		for (int i = 0; i < nact; i++) {
			childern[i] = NULL;
			pi[i] = 0.0;
			qu[i] = 0.0;
			uu[i] = 0.0;
			nn[i] = 0;
			done[i] = false;
		}
		sumN = 0;
        from -> seen.copyTo(seen); // Maybe optimized to use the same seen object instead of copying it..
	}

shadow::shadow(shadow* from) :
	verbosity                     (from -> verbosity),
    ccmin_mode                    (from -> ccmin_mode),
    phase_saving                  (from -> phase_saving),
    learntsize_factor             (from -> learntsize_factor),
    learntsize_inc                (from -> learntsize_inc),
    garbage_frac                  (from -> garbage_frac),
    max_learnts                   (from -> max_learnts),
    learntsize_adjust_start_confl (from -> learntsize_adjust_start_confl),
    learntsize_adjust_inc         (from -> learntsize_adjust_inc),
    learntsize_adjust_confl       (from -> learntsize_adjust_confl),
    learntsize_adjust_cnt         (from -> learntsize_adjust_cnt),
    cla_inc                       (from -> cla_inc),
    clause_decay                  (from -> clause_decay),
    trail_size                    (from -> trail_size),
    trail_lim_size				  (from -> trail_lim_size),
    qhead                         (from -> qhead),
    ca_size						  (from -> ca_size),
    learnts_size 				  (from -> learnts_size)
	{
		ca_shadow.extra_clause_field = from -> ca_shadow.extra_clause_field; 
		learnts_copy_is_uninitialized = true;
		origin = NULL;
		parent = from;
		index_child_last_pick = -1;
		for (int i = 0; i < nact; i++) {
			childern[i] = NULL;
			pi[i] = 0.0;
			qu[i] = 0.0;
			uu[i] = 0.0;
			nn[i] = 0;
			done[i] = false;
		}
		sumN = 0;
        from->seen.copyTo(seen);
	}

shadow::~shadow() {
	// destruct the data in std::unordered_map<int, vec<Solver::Watcher>* > watches_map;
	for (std::pair<int, vec<Solver::Watcher>* > element : watches_map) {
		element.second -> ~vec();
    }
}

// This function set the child at index action to be the new root of MCTS
// the connection between old root and new root is remove, and new root -> origin is set as Solver
// The pointer to the child at index action is returned (to assign to the root_shadow)
shadow* shadow::next_root(int action) {
    if (childern[action] == NULL) // the new root doesn't exist because it is in a finished state
        return NULL;
    childern[action] -> origin = origin;
    childern[action] -> parent = NULL;
    shadow* temp = childern[action];
    childern[action] = NULL;
    return temp;
}

// This function push forward the search within the MCTS forward
// The key logic is picking the best childern index, which is calculated based on nn, pi, qu and sumN
// The logic also prevent picking variables whose values are already assigned.
// it will throw exception if no good choice can be make (check the assert command)
// if a child index is pick, update the nn and sumN, but the qu (values) has to wait until next simulation call from Solver (need neural net evaluation)
// if the child to pick is marked done, it means that the child was visited before, and it stepped into finished state. Return NULL
// if the child to pick is not NULL, make recursive call from that child
// if the child to pick is NULL (given that done is not NULL), construct a new shadow copy for that child and ask the child to step on the index 
// step() write the new state to array argument, and returns whether the state is done (if done, return false)
// if the child stepped to "finished state", call its destructor, and set childern[index] as NULL (avoid dangling pointers)
// return childern[index] (could be NULL if the child stepped to finished state) (otherwise, the returned pointer is to the leaf_shadow whose pi needs evaluation)
// NOTE: index_child_last_pick is a field in this object, which tracks the most recent pick of childern IMPORTANT for assigning qu and pi later!!!
shadow* shadow::next_to_explore(float* array) {
	// pick a child to simulate by the score system
	index_child_last_pick = 0; float pick_val = -1e20;
	for (int i = 0; i < nact; i++) {
		uu[i] = c_act * pi[i] * sqrt(sumN + 1) / (1 + nn[i]);
		float val = nn[i] == 0? uu[i] : uu[i] + qu[i] / nn[i];
		if (val > pick_val && value(toLit(i)) == l_Undef ) { // IMPORTANT: only pick Lit that is not yet assigned
			index_child_last_pick = i; pick_val = val;
		}
	}
	assert (pick_val > (0.1 - 1e-20) && "failed to pick a good action for simulation");

	// found a child to simulate
    // printf("simulation steps on %d\n", index_child_last_pick);
	nn[index_child_last_pick] += 1; sumN++;
	if (done[index_child_last_pick]) { // the picked child is already visited before and the child is in a done state
		return NULL;
	}
	if (childern[index_child_last_pick] != NULL) {
		return childern[index_child_last_pick] -> next_to_explore(array);
	} else {
		childern[index_child_last_pick] = new shadow(this);
        done[index_child_last_pick] = !(childern[index_child_last_pick] -> step(toLit(index_child_last_pick), array));
		if (done[index_child_last_pick]) {
			childern[index_child_last_pick] -> ~shadow();
			childern[index_child_last_pick] = NULL;
		}
		return childern[index_child_last_pick];
	}
}

// this function writes nn (visit count) to array (some numpy array provided by RL algorithm)
void shadow::get_visit_count(float* array) {
	for (int i = 0; i < nact; i++)
		array[i] = nn[i];
}

// helper function for write_clause (return true if a Clause c is already satisfied)
bool shadow::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; 
}
// helper function for generate_state (write state to a 1D array and returns the next col to write to)
int shadow::write_clause(const Clause& c, int index_col, float* array) {
	if (satisfied(c)) return index_col;
    for (int i = 0; i < c.size(); i++) {
        if (value(c[i]) != l_False) {
        	int index_row = var(c[i]); int index_z = int(sign(c[i]));
        	int index = index_z + index_row * dim2 + index_col * dim1 * dim2;
        	array[index] = 1.0;
        }
    }
    return index_col + 1;
}
// write state in tensor "array" as side effect (if too many clauses to write, cut off by dim0)
// return true if state is not empty (not solved), false otherwise
bool shadow::generate_state(float* array) {
	int index_col = 0;
	// write clauses in array (optimization:: get hold of the Solver's clause and work from there)
	shadow* temp = this;
	while (temp -> parent != NULL) {
        temp = temp -> parent;
    }
	Solver* solver = temp -> origin;
	vec<CRef>& clauses = solver -> clauses;
	ClauseAllocator& ca = solver -> ca;
	for (int i = 0; i < clauses.size() && index_col < dim0; i++) 
        index_col = write_clause(ca[clauses[i]], index_col, array);
    // write learnts in array
    for (int i = 0; i < get_learnts_size() && index_col < dim0; i++)
    	index_col = write_clause(get_clause(get_learnts(i)), index_col, array);
    /* printf("clause %d, learnts %d\n", clauses.size(), get_learnts_size());
    for (int i = 0; i < trail_size; i++) {
        printf("%d_", get_trail(i).x);
    }
    printf("\n");
    for (int i = 0; i < 20; i++) {
        printf("%d_", get_assigns(i));
    }
    printf("\n"); */
    return index_col > 0;
}

// this function return true if state is not solved, false otherwise (note no parameters)
bool shadow::generate_state() {
	shadow* temp = this;
	while (temp -> parent != NULL) temp = temp -> parent;
	Solver* solver = temp -> origin;
	vec<CRef>& clauses = solver -> clauses;
	ClauseAllocator& ca = solver -> ca;
	for (int i = 0; i < clauses.size(); i++) {
		if (!satisfied(ca[clauses[i]])) return true;
	}
	for (int i = 0; i < get_learnts_size(); i++) {
		if (!satisfied(get_clause(get_learnts(i)))) return true;
	}
	return false;
}

// enque p as the next assignment. The qhead didn't increment, so the propagate() knows that this Lit p needs to be propagated.
void shadow::uncheckedEnqueue(Lit p, CRef from)
{
    assert(value(p) == l_Undef);
    set_assigns(var(p), lbool(!sign(p)));
    set_vardata(var(p), Solver::mkVarData(from, decisionLevel()));
    append_trail(p);
}

// propagate the newly assigned Lits
CRef shadow::propagate()
{
    CRef    confl     = CRef_Undef;
    while (qhead < trail_size){
        Lit                   p  = get_trail(qhead++); // 'p' is enqueued fact to propagate.
        vec<Solver::Watcher>& ws = get_watches_copied(p);

        if (get_dirty(p)) clean_watches(p); // Comments by Fei: the initial lookup function garantees that the watcher list is cleaned!!
        Solver::Watcher *i, *j, *end;

        for (i = j = (Solver::Watcher*)ws, end = i + ws.size();  i != end;){

            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; 
            }

            // Make sure the false literal is data[1]:
            CRef cr       = i->cref; 
            Lit false_lit = ~p;
            Lit first;
            const Clause& c = get_clause(cr);
            if (c[0] == false_lit) {
            	Clause& cc = get_clause_copied(cr);
            	first = cc[0] = cc[1]; cc[1] = false_lit;
            } else { 
  				assert(c[1] == false_lit);
  				first = c[0];
  			}
            i++; 

            // If 0th watch is true, then clause is already satisfied.
            Solver::Watcher w = Solver::Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; 
            }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                	Clause& cc = get_clause_copied(cr);
                	cc[1] = cc[k]; cc[k] = false_lit;
                    get_watches_copied(~cc[1]).push(w);
                    goto NextClause; 
                }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr; 
                qhead = trail_size;
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            } else
                uncheckedEnqueue(first, cr);

        NextClause:;
        }
        ws.shrink(i - j);
    }
    return confl;
}

struct reduceDB_lt {
	shadow* which_shadow;
	reduceDB_lt(shadow* this_shadow): which_shadow(this_shadow) {}
	bool operator() (CRef x, CRef y) {
		const Clause& a = which_shadow->get_clause(x);
		const Clause& b = which_shadow->get_clause(y);
		return a.size() > 2 && (b.size() == 2 || a.activity() < b.activity());
	} 
};
void shadow::reduceDB()
{
    int     i, j;
    double  extra_lim = cla_inc / get_learnts_size();    // Remove any clause below this activity
	
	// need a copy of learnts!
	get_copy_for_learnts();
	sort(learnts_copy, reduceDB_lt(this));

    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
	for (i = j = 0; i < get_learnts_size(); i++) {
		const Clause& c = get_clause(learnts_copy[i]);
		if (c.size() > 2 && !locked(c) && (i < get_learnts_size() / 2 || c.activity() < extra_lim))
			removeClause(learnts_copy[i]); 
		else 
			learnts_copy[j++] = learnts_copy[i];
	}

    learnts_copy.shrink(i - j); // learnts_size = learnts_copy.size(); 
    // checkGarbage(); // Comments by Fei: disabled for now for small SAT problems
}

void shadow::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel) {
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    out_learnt.push();      // (leave room for the asserting literal)
    int index = trail_size - 1;
    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        const Clause& c = get_clause(confl); 

        if (c.learnt()) {
        	Clause& cc = get_clause_copied(confl);
            claBumpActivity(cc);
        }
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
            Lit q = c[j];
            if (!seen[var(q)] && get_level(var(q)) > 0){ 
                // varBumpActivity(var(q));
                seen[var(q)] = 1;
                if (get_level(var(q)) >= decisionLevel())
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }
        // Select next clause to look at:
        while (!seen[var(get_trail(index--))]);
        p     = get_trail(index+1);
        confl = get_reason(var(p));
        seen[var(p)] = 0;
        pathC--;
    } while (pathC > 0);
    out_learnt[0] = ~p;
    
    // Simplify conflict clause:
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        for (i = j = 1; i < out_learnt.size(); i++)
            if (get_reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i])) 
                out_learnt[j++] = out_learnt[i];
        
    } else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (get_reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else {
                const Clause& c = get_clause(get_reason(var(out_learnt[i])));
                for (int k = 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && get_level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; 
                    }
            }
        }
    } else
        i = j = out_learnt.size();

    out_learnt.shrink(i - j);
    // max_literals += out_learnt.size(); remove code for stats
    // tot_literals += out_learnt.size(); remove code for stats
    
    // Find correct backtrack level:
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else {
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (get_level(var(out_learnt[i])) > get_level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = get_level(var(p));
    }

    for (int j = 0; j < analyze_toclear.size(); j++) 
    	seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}

// Check if 'p' can be removed from a conflict clause.
bool shadow::litRedundant(Lit p) {
    enum { seen_undef = 0, seen_source = 1, seen_removable = 2, seen_failed = 3 };
    assert(seen[var(p)] == seen_undef || seen[var(p)] == seen_source);
    assert(get_reason(var(p)) != CRef_Undef);

    const Clause* c = &get_clause(get_reason(var(p)));
    vec<ShrinkStackElem>& stack = analyze_stack;
    stack.clear();

    for (uint32_t i = 1; ; i++){
        if (i < (uint32_t)c->size()){
            // Checking 'p'-parents 'l':
            Lit l = (*c)[i];
            
            // Variable at level 0 or previously removable:
            if (get_level(var(l)) == 0 || seen[var(l)] == seen_source || seen[var(l)] == seen_removable) 
                continue; 
            
            // Check variable can not be removed for some local reason:
            if (get_reason(var(l)) == CRef_Undef || seen[var(l)] == seen_failed) {
                stack.push(ShrinkStackElem(0, p));
                for (int i = 0; i < stack.size(); i++)
                    if (seen[var(stack[i].l)] == seen_undef){
                        seen[var(stack[i].l)] = seen_failed;
                        analyze_toclear.push(stack[i].l);
                    }
                return false;
            }

            // Recursively check 'l':
            stack.push(ShrinkStackElem(i, p));
            i  = 0;
            p  = l;
            c  = &get_clause(get_reason(var(p)));
        } else {
            // Finished with current element 'p' and reason 'c':
            if (seen[var(p)] == seen_undef){
                seen[var(p)] = seen_removable;
                analyze_toclear.push(p);
            }

            // Terminate with success if stack is empty:
            if (stack.size() == 0) break;
            
            // Continue with top element on stack:
            i  = stack.last().i;
            p  = stack.last().l;
            c  = &get_clause(get_reason(var(p)));

            stack.pop();
        }
    }
    return true;
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
void shadow::cancelUntil(int level) {
    if (decisionLevel() > level) {
        for (int c = trail_size - 1; c >= get_trail_lim(level); c--) {
            Var x  = var(get_trail(c));
            set_assigns(x, l_Undef);
            if (phase_saving > 1 || (phase_saving == 1 && c > get_trail_lim(trail_lim_size - 1)))
                set_polarity(x, sign(get_trail(c)));
            // insertVarOrder(x);  remove code related with ordering
        }
        qhead = get_trail_lim(level);
        trail_map_clear_until(qhead);
        trail_lim_map_clear_until(level);
    } 
}

void shadow::attachClause(CRef cr){
    const Clause& c = get_clause(cr);
    assert(c.size() > 1);
    get_watches_copied(~c[0]).push(Solver::Watcher(cr, c[1]));
    get_watches_copied(~c[1]).push(Solver::Watcher(cr, c[0]));
    //if (c.learnt()) num_learnts++, learnts_literals += c.size();
    //else            num_clauses++, clauses_literals += c.size(); No need to update stats
}


void shadow::detachClause(CRef cr, bool strict){ // remove clause from watcher list
    const Clause& c = get_clause(cr);
    assert(c.size() > 1);
    
    // Strict or lazy detaching:
    if (strict){
    	remove(get_watches_copied(~c[0]), Solver::Watcher(cr, c[1]));
    	remove(get_watches_copied(~c[1]), Solver::Watcher(cr, c[0]));
    }else{
    	set_dirty(~c[0], 1);
    	set_dirty(~c[1], 1);
    }
//    if (c.learnt()) num_learnts--, learnts_literals -= c.size();
//    else            num_clauses--, clauses_literals -= c.size(); Comments: no need to track stats
}

void shadow::removeClause(CRef cr) {
	Clause& c = get_clause_copied(cr);
	detachClause(cr);

    // Don't leave pointers to free'd memory!
    if (locked(c)) set_vardata(var(c[0]), Solver::mkVarData(CRef_Undef, get_vardata(var(c[0])).level));
    c.mark(1);
    ca_shadow.free(cref_map.at(cr)); 
}


// this function tries to step on "action", and write the new state to array
// it returns true if the state is not yet finished and array has non-zero values in it.
// NOTE: restart and garbage collection are disabled in this function, and similar functionality has to be disabled in Solver::search as well to 
//       keep the simulation step consistent with the actual step
bool shadow::step(Lit action, float* array)
{ 
	newDecisionLevel();
    uncheckedEnqueue(action);

    while (true) {
        CRef confl = propagate(); // Comments by Fei: changed to field variable in Solver, USE local variable in shadow
        if (confl != CRef_Undef) {
            // CONFLICT
            // conflicts++; conflictCounts++; // Comments by Fei: replace usage! conflictC++; 
            if (decisionLevel() == 0) return false; // terminate with UNSAT, return false because nothing written in the array argument for evaluation 

            vec<Lit> learnt_clause; // Comments by Fei: make this variable local to loop (used and destroyed)
            learnt_clause.clear();
            int backtrack_level; // Comments by Fei: make this variable local to loop (used and destroyed)
            analyze(confl, learnt_clause, backtrack_level); // this function writes to learnts_clause and backtrack_level arguments
            cancelUntil(backtrack_level);

            if (learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{ 
            	CRef cr = get_alloc(learnt_clause, true);
            	append_learnts(cr);
                attachClause(cr);
                claBumpActivity(get_clause_copied(cr)); 
                uncheckedEnqueue(learnt_clause[0], cr);
            }

            // varDecayActivity();
            claDecayActivity(); 

            if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
                max_learnts             *= learntsize_inc;
            }

        } else {
            // NO CONFLICT (disabled restart) (disable simplify() when decisionLevel == 0)
            if (get_learnts_size() - nAssigns() >= max_learnts)
                reduceDB();

            // remove code about assumptions. 
			// save states and return true is state is not finished, false otherwise
            return generate_state(array); 
        }
    }
}

// Garbage Collection methods: Disabled for now
void shadow::relocAll(ClauseAllocator& to)
{
	/*
    // All watchers:
    watches.cleanAll();
    for (int v = 0; v < nVars(); v++)
        for (int s = 0; s < 2; s++) {
            Lit p = mkLit(v, s);
            vec<Watcher>& ws = watches[p];
            for (int j = 0; j < ws.size(); j++)
                ca.reloc(ws[j].cref, to);
        }

    // All reasons:
    for (int i = 0; i < trail.size(); i++) {
        Var v = var(trail[i]);
        // Note: it is not safe to call 'locked()' on a relocated clause. This is why we keep
        // 'dangling' reasons here. It is safe and does not hurt.
        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)]))){
            assert(!isRemoved(reason(v)));
            ca.reloc(vardata[v].reason, to);
        }
    }

    // All learnt:
    //
    int i, j;
    for (i = j = 0; i < learnts.size(); i++)
        if (!isRemoved(learnts[i])){
            ca.reloc(learnts[i], to);
            learnts[j++] = learnts[i];
        }
    learnts.shrink(i - j);

    // All original:
    //
    for (i = j = 0; i < clauses.size(); i++)
        if (!isRemoved(clauses[i])){
            ca.reloc(clauses[i], to);
            clauses[j++] = clauses[i];
        }
    clauses.shrink(i - j);*/
}

void shadow::garbageCollect()
{   /*
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted()); 

    relocAll(to);
    if (verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
               ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
    to.moveTo(ca);
    */
}