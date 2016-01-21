#ifndef ACGT_H
#define ACGT_H
#include <map>
#include <vector>
#include "load.h"

using namespace std;
void get_average_ACGT_profile(map<string, vector<segment> > S, 
	vector<PSSM *> PSSMS, int pad ,map<string,  double  >   & NN, map<string, double [3000][4]> & G );


#endif
