#ifndef CSS_H
#define CSS_H
#include "load.h"
#include <map>
#include <vector>
#include <string>
using namespace std;

void collect_sample_stats(map<string, vector<segment>>,
 vector<PSSM *>,
	map<int, vector<double> > &,
	 map<int, vector<double>> &, map<int, map<int, double> > &,int );

#endif
