#ifndef out_H
#define out_H
#include <map>
#include <vector>
#include "load.h"
using namespace std;
void write_observation_stats(map<string, vector<segment>>,
	string, string, vector<PSSM *>,map<string, vector<vector<double> > >, int);
#endif