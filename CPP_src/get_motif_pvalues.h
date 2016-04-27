#ifndef gmp_H
#define gmp_H
#include "load.h"
vector<PSSM *> DP_pvalues(vector<PSSM *>,int, vector<double>,bool);

vector<PSSM *> construct_position_specific_pvalues(vector<PSSM *>, int , 
	vector<vector<double>> &  , vector<vector<double>> & , int );


#endif
