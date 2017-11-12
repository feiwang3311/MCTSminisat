#ifndef HYPER_CONST
#define HYPER_CONST

class Hyper_Const
{
public :
	static const int dim0 = 100;           // max_clause
    static const int dim1 = 20;           // max_var
    static const int dim2 = 2;           // nc
    static const int nact = 40;           // nact 
    static const float c_act;     	 // c_act is a hyperparameter for MCTS (decide the level of exploration) 
    static const int MCTS_size_lim = 100; // the size of MCT we want to achieve.
};

#endif