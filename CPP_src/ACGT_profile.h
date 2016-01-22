#ifndef ACGT_H
#define ACGT_H
#include <map>
#include <vector>
#include "load.h"

using namespace std;
map<string, double [2000][4]> get_average_ACGT_profile(map<string, vector<segment> > S, 
	vector<PSSM *> PSSMS, int pad ,map<string,  double  >   & NN, map<string, double [2000][4]> & G );


#endif
