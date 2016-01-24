#ifndef SIM_H
#define SIM_H
#include <map>
#include <string>
#include <vector>
#include "load.h"
using namespace std;
map<int, map<int, vector<segment> >> run_sims(map<int, double [2000][4]> ,
 map<int, double>,  vector<PSSM *> , int,int );
#endif