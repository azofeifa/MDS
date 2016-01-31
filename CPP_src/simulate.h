#ifndef SIM_H
#define SIM_H
#include <map>
#include <string>
#include <vector>
#include "load.h"
using namespace std;
void run_sims(map<int, double [2000][4]> ,
 map<int, double>,  vector<PSSM *> , int,int,
 vector<double>, double,
 map<int, vector<vector<double> >> & );
void run_sims2(map<string, vector<segment>> , vector<PSSM *> ,int  , int  , 
	vector<double>  , double, 
	map<int, vector<vector<double> >> & );
#endif