#ifndef Minisat_shadow_h
#define Minisat_shadow_h

#include "minisat/mtl/Vec.h"
#include "minisat/mtl/Heap.h"
#include "minisat/mtl/Alg.h"
#include "minisat/mtl/IntMap.h"
#include "minisat/utils/Options.h"
#include "minisat/core/SolverTypes.h"
#include "minisat/core/Const.h"
#include <unordered_map>

namespace Minisat {

class Solver; 
class shadow {
public:
    shadow( shadow* from); 
    shadow( Solver* from); 
    virtual ~shadow();

    const int nact = Hyper_Const::nact;
    const int dim0 = Hyper_Const::dim0;
    const int dim1 = Hyper_Const::dim1;
    const int dim2 = Hyper_Const::dim2;
    const float c_act = Hyper_Const::c_act;
    const int MCTS_size_lim = Hyper_Const::MCTS_size_lim;

    Solver* origin;                      // this points to the Solver instance that these shadows are cloned from (will be null if not root node)
    shadow* parent;                      // this points to the parent shadow node (will be null if this is the root node)
    int index_child_last_pick;           // this field remembers the child of choice during MCTS, for assigning Q after passing the state in neural net.
    shadow* childern[Hyper_Const::nact]; // this is an array of pointers to all the childern of this node (could be null if childern not visited yet)
    float pi[Hyper_Const::nact];         // this is an array of pi values (initial visiting probality for MCTS simulation)
    float qu[Hyper_Const::nact];         // this is an array of qu values (total expect score for MCTS simulation)
    int   nn[Hyper_Const::nact];         // this is an array of nn values (total visit count for each childern in MCTS simulation)
    float uu[Hyper_Const::nact];         // this is an array of U values (combine pi and nn values)
    bool done[Hyper_Const::nact];        // this is an array to label if a child branch leads to finished state
    int sumN;                            // this is the total number of MCTS simulations run from this node (sum of nn)
    
    // MCTS functions
    shadow* next_root(int action); // this function set child at index "action" to be the next root, it returns the pointer to the new root
    void    get_visit_count(float* count); // this function writes the nn array to array argument
    shadow* next_to_explore(float* state); // this function initiate simulation from this shadow, will write state to state argument, returns leaf shadow 

    bool generate_state(float*); // this function askes this node to write its state to the argument given by RL algorithm (no memory copy, inplace write)
    bool generate_state();       // overload of generate_state function. No float pointer to write to, so only returns true if state is not done.
    int  write_clause (const Clause& c, int index_col, float* array); // helper function of "generate_state" for a clause
    bool satisfied    (const Clause& c) const;                        // helper function of "write_clause"


    // Mode of operation: Directly copied from Solver class. Used during simulation 
    int       verbosity;          
    int       ccmin_mode;         // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
    int       phase_saving;       // Controls the level of phase saving (0=none, 1=limited, 2=full). 
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses. (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.  (default 1.1) 
    double    garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered. 
    double    max_learnts; 
    int       learntsize_adjust_start_confl;
    double    learntsize_adjust_inc;
    double    learntsize_adjust_confl; 
    int       learntsize_adjust_cnt; 
    double    cla_inc;            // Amount to bump next clause with. 
    double    clause_decay;     


    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is used, 
    // exept 'seen' wich is used in several places. Directly copied from Solver class. Used during simulation
    struct ShrinkStackElem {
        uint32_t i;
        Lit      l;
        ShrinkStackElem(uint32_t _i, Lit _l) : i(_i), l(_l){}
    };
    VMap<char>           seen;
    vec<ShrinkStackElem> analyze_stack;
    vec<Lit>             analyze_toclear;
    vec<Lit>             add_tmp;

    
    // fields and methods for caching the difference 
    // NOTE: data structures here are never directly accessed. Instead, use the getter and setter functions. 
    // The getter and setter functions take care of the logic of only storing the difference
    // Implementation of these functions are inlined in this file
    std::unordered_map<int, Lit> trail_map;
    int trail_size; 
    Lit get_trail(int x) const; 
    void append_trail(Lit y); 
    void trail_map_clear_until(int x);                    // corresponding to trail.shrink to x

    std::unordered_map<int, int> trail_lim_map;
    int trail_lim_size; 
    int get_trail_lim(int x) const; 
    void append_trail_lim(int y); 
    void trail_lim_map_clear_until(int level);            // corresponding to trail_lim.shrink to level

    int qhead;

    std::unordered_map<Var, lbool> assigns_map; 
    lbool get_assigns(Var x) const;        
    void set_assigns(Var x, lbool y); 

    std::unordered_map<Var, Solver::VarData> vardata_map; // Stores reason and level for each variable.
    Solver::VarData get_vardata(Var x) const;
    void set_vardata(Var x, Solver::VarData data); 

    std::unordered_map<Var, char> polarity_map;           // The preferred polarity of each variable.
    char get_polarity(Var x) const;
    void set_polarity(Var x, char y);

    std::unordered_map<int, vec<Solver::Watcher>* > watches_map; 
    std::unordered_map<int, char> dirty_map;
    vec<Solver::Watcher>& get_watches_copied(Lit p);
    char get_dirty(Lit p) const;
    void set_dirty(Lit p, char c);
    void clean_watches(Lit p);                            // remove deleted clause (mark is true) from watcher list 

    ClauseAllocator ca_shadow;                            // if clauses are changed, they are copied to ca_shadow, then changed from here
    std::unordered_map<CRef, CRef> cref_map;              // copied clauses often have new CRef. Use this as a mapping from old CRef to new CRef
    int ca_size;                                          // this tracks the size of ca of the parent (used as CRef when adding learnt clauses)
    const Clause& get_clause(CRef cr) const;              // this function gets a reference to Clause at CRef cr. No modification allowed.
    Clause& get_clause_copied(CRef cr);                   // this function gets a reference to a copied Clause at CRef cr. Modification allowed.
    CRef get_alloc(const vec<Lit>& ps, bool learnt = false);

    std::unordered_map<int, CRef> learnts_map;
    vec<CRef> learnts_copy;                               // if reduceDB gets called, it is easier to copy the whole learnts vector
    bool learnts_copy_is_uninitialized;                   // this is true unless reduceDB gets called. 
    int learnts_size; 
    CRef get_learnts(int x) const;
    void append_learnts(CRef y);
    int  get_learnts_size() const;
    void get_copy_for_learnts();


    // main functions:
    bool     step             (Lit action, float* array);               // make a simulation step toward action, write state to array. return true of array is not empty
    CRef     propagate        ();                                       // Perform unit propagation. Returns possibly conflicting clause.
    void     analyze          (CRef confl, vec<Lit>& learnt, int& bt);  // (bt = backtrack)
    void     cancelUntil      (int level);                              // Backtrack until a certain level.
    void     reduceDB         ();                                       // Reduce the set of learnt clauses.
   
    // other helper functions
    int      get_level        (Var x)   const;
    CRef     get_reason       (Var x)   const;
    lbool    value            (Var x)   const;         // The current value of a variable.
    lbool    value            (Lit p)   const;         // The current value of a literal.
    int      decisionLevel    ()        const;         // Gives the current decisionlevel.
    void     newDecisionLevel ();                                                      // Begins a new decision level.
    void     uncheckedEnqueue (Lit p, CRef from = CRef_Undef);                         // Enqueue a literal. Assumes value of literal is undefined.
    bool     litRedundant     (Lit p);                                                 // (helper method for 'analyze()'
    int      nAssigns         ()        const;         // The current number of assigned literals. 

    // clause_related helper function
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.
    void     attachClause     (CRef cr);                      // Attach a clause to watcher lists.
    void     detachClause     (CRef cr, bool strict = false); // Detach a clause to watcher lists.
    void     removeClause     (CRef cr);                      // Detach and free a clause.                                      
    bool     locked           (const Clause& c) const;        // Returns TRUE if a clause is a reason for some implication in the current state.

    // memory helper functions (current implementation disabled garbage collection)
    virtual void garbageCollect();
    void     checkGarbage(double gf);
    void     checkGarbage();
    uint32_t get_ca_size();
    void     relocAll (ClauseAllocator& to);
};

// inline helper functions
inline int   shadow::nAssigns        ()      const { return trail_size; }
inline void  shadow::newDecisionLevel()            { append_trail_lim(trail_size); }
inline int   shadow::decisionLevel   ()      const { return trail_lim_size; }
inline lbool shadow::value           (Var x) const { return get_assigns(x); } 
inline lbool shadow::value           (Lit p) const { return get_assigns(var(p)) ^ sign(p); } 
inline int   shadow::get_level       (Var x) const { return get_vardata(x).level; }
inline CRef  shadow::get_reason      (Var x) const { return get_vardata(x).reason; }
inline void  shadow::claBumpActivity (Clause& c) {
    if ( (c.activity() += cla_inc) > 1e20 ) {
        // Rescale:
        for (int i = 0; i < get_learnts_size(); i++) {
            get_clause_copied(get_learnts(i)).activity() *= 1e-20;
        }
        cla_inc *= 1e-20; 
    } 
}
inline void  shadow::claDecayActivity()                      { cla_inc *= (1 / clause_decay); }
inline bool  shadow::locked          (const Clause& c) const { 
    return value(c[0]) == l_True && get_reason(var(c[0])) != CRef_Undef && &get_clause(get_reason(var(c[0]))) == &c; 
    // ca.lea(get_reason(var(c[0]))) == &c;  NOTE: not sure if this is equivalent change
}    


// inline functions for difference in state
// IMPORTANT: if parent == NULL, origin must not be NULL, and differences in maps are ignored!! 
inline Lit shadow::get_trail(int x) const { 
    assert ( x < trail_size && "get_trail out of range");
    const shadow* temp = this;
    while(temp->trail_map.count(x) == 0 && temp-> parent != NULL) temp = temp->parent;
    if (temp -> parent == NULL) return temp->origin->trail[x];
    return temp->trail_map.at(x);
} 
inline void shadow::append_trail(Lit y) {
    trail_map[trail_size++] = y;
} 
inline void shadow::trail_map_clear_until(int x){
    trail_size = x;
    for(auto it = trail_map.begin(); it != trail_map.end(); )
        if(it->first >= x)
            it = trail_map.erase(it);
        else
            ++it;
}

inline int shadow::get_trail_lim(int x) const {
    assert (x < trail_lim_size && "get_trail_lim out of range");
    const shadow* temp = this;
    while(temp->trail_lim_map.count(x) == 0 && temp -> parent != NULL) temp = temp -> parent;
    if (temp -> parent == NULL) return temp-> origin -> trail_lim[x];
    return temp -> trail_lim_map.at(x);
}
inline void shadow::append_trail_lim(int y) {
    trail_lim_map[trail_lim_size++] = y;
}
inline void shadow::trail_lim_map_clear_until(int level) { // corresponding to trail_lim.shrink to level
    trail_lim_size = level;
    for (auto it = trail_lim_map.begin(); it != trail_lim_map.end(); )
        if (it -> first >= level)
            it = trail_lim_map.erase(it);
        else
            ++it;
} 

inline lbool shadow::get_assigns(Var x) const { 
    const shadow* temp = this;
    while(temp->assigns_map.count(x) == 0 && temp->parent != NULL) temp = temp->parent;
    if (temp -> parent == NULL) return temp->origin->assigns[x];
    return temp->assigns_map.at(x);
}
inline void shadow::set_assigns(Var x, lbool y) {
    assigns_map[x] = y;
}

inline Solver::VarData shadow::get_vardata(Var x) const {
    const shadow* temp = this;
    while (temp->vardata_map.count(x) == 0 && temp -> parent != NULL) temp = temp -> parent;
    if (temp -> parent == NULL) return temp-> origin -> vardata[x];
    return temp -> vardata_map.at(x);
}
inline void shadow::set_vardata(Var x, Solver::VarData data) {
    vardata_map[x] = data;
}

inline char shadow::get_polarity(Var x) const {
    const shadow* temp = this;
    while (temp->polarity_map.count(x) == 0 && temp -> parent != NULL) temp = temp -> parent;
    if (temp -> parent == NULL) return temp->origin->polarity[x];
    return temp-> polarity_map.at(x);
}
inline void shadow::set_polarity(Var x, char y) {
    polarity_map[x] = y;
}

inline vec<Solver::Watcher>& shadow::get_watches_copied(Lit p_input) {
    int p = p_input.x;
    if (!watches_map.count(p)) {
        watches_map[p] = new vec<Solver::Watcher>();
        shadow* temp = this -> parent; 
        while(temp->watches_map.count(p) == 0 && temp-> parent != NULL) temp = temp -> parent;
        if (temp -> parent == NULL) {
            temp->origin->watches.lookup(p_input).copyVstructTo(*watches_map.at(p)); // Comments by Fei: this is cleaned!
            dirty_map[p] = 0;
        } else {
            temp->watches_map.at(p)->copyVstructTo(*watches_map.at(p)); // Comments by Fei: this may not be clean!!
            dirty_map[p] = temp->dirty_map.at(p);
        }
    }
    return *watches_map.at(p);
}
inline char shadow::get_dirty(Lit p_input) const {
    int p = p_input.x;
    const shadow* temp = this;
    while (temp -> dirty_map.count(p) == 0 && temp -> parent != NULL) temp = temp->parent;
    if (temp -> parent == NULL) return 0; // Comments by Fei: copy vec<Watcher> from origin is bound to be clean.
    return temp ->dirty_map.at(p); 
}
inline void shadow::set_dirty(Lit p, char c) {
    dirty_map[p.x] = c;
}
inline void shadow::clean_watches(Lit p) {
    vec<Solver::Watcher>& ws = *watches_map.at(p.x);
    int  i, j;
    for (i = j = 0; i < ws.size(); i++)
        if (!get_clause(ws[i].cref).mark())
            ws[j++] = ws[i];
    ws.shrink(i - j);
    set_dirty(p, 0);
}

inline const Clause& shadow::get_clause(CRef cr) const {
    const shadow* temp = this;
    while (temp->cref_map.count(cr) == 0 && temp->parent != NULL) temp = temp->parent;
    if (temp -> parent == NULL) return temp->origin->ca[cr];
    return temp->ca_shadow[temp->cref_map.at(cr)];
}
inline Clause& shadow::get_clause_copied(CRef cr) {
    if (cref_map.count(cr)) return ca_shadow[cref_map.at(cr)];
    CRef mapped_cr = ca_shadow.alloc(get_clause(cr));
    cref_map[cr] = mapped_cr;
    return ca_shadow[mapped_cr];
}
inline CRef shadow::get_alloc(const vec<Lit>& ps, bool learnt) {
    CRef outside_value = (CRef)ca_size;
    CRef inside_value = ca_shadow.alloc(ps, learnt);
    cref_map[outside_value] = inside_value;
    ca_size += (ca_shadow.size() - (uint32_t)inside_value);
    return outside_value;
}

inline CRef shadow::get_learnts(int x) const {
    const shadow* temp = this;
    while (temp -> learnts_copy_is_uninitialized && temp -> learnts_map.count(x) == 0 && temp -> parent != NULL) temp = temp -> parent;
    if (temp -> parent == NULL) return temp->origin->learnts[x];
    if (!temp -> learnts_copy_is_uninitialized) return temp->learnts_copy[x];
    return temp-> learnts_map.at(x);
}
inline void shadow::append_learnts(CRef y) {
    assert (learnts_copy_is_uninitialized && "should never been in this state (append to copy)");
    learnts_map[learnts_size++] = y;
}
inline int shadow::get_learnts_size() const {
    if (parent == NULL) return origin -> learnts.size();
    return learnts_copy_is_uninitialized ? learnts_size : learnts_copy.size();
}
inline void shadow::get_copy_for_learnts() {
    assert (learnts_copy_is_uninitialized && "error state (call copy learnts while learnts is already copied");
    learnts_copy.clear();
    learnts_copy.growTo(learnts_size);
    for (vec<CRef>::Size i = 0; i < learnts_size; i++) learnts_copy[i] = get_learnts(i);
    learnts_copy_is_uninitialized = false; 
}


// Memory management functions (currently disabled)
inline uint32_t shadow::get_ca_size() {
    shadow* temp = this;
    while(temp->parent != NULL) temp = temp-> parent;
    return temp->origin->ca.size();
}
inline void shadow::checkGarbage(void){ return checkGarbage(garbage_frac); }
inline void shadow::checkGarbage(double gf){ 
    if (ca_shadow.wasted() > (ca_shadow.size() + get_ca_size()) * gf)
        garbageCollect(); 
    }
}

#endif