#ifndef ACGT_H
#define ACGT_H
#include <map>
#include <vector>
#include "load.h"

using namespace std;
map<string, vector<vector<double> > > get_average_ACGT_profile(map<string, vector<segment> > S, 
	vector<PSSM *> PSSMS, int pad );


#endif
