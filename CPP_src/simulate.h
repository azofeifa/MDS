#ifndef SIM_H
#define SIM_H
#include <map>
#include <string>
#include <vector>
#include "load.h"
using namespace std;
map<string, vector<vector<segment> >  > run_sims(map<string, vector<vector<double> > > ,
 map<string, double>, int );
#endif