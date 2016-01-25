#ifndef scanner_H
#define scanner_H
#include <vector>
#include <map>
#include "load.h"
using namespace std;
void run_scanner();
vector<double> get_GC_content(map<string, vector<segment>>);
map<string, vector<segment> > run_accross(map<string, vector<segment>>, 
	vector<PSSM *>,vector<double>, double,int);
map<int, map<int, vector<segment> >> scan_simulations(map<int, map<int, vector<segment> >> ,
	vector<PSSM *> , vector<double> ,double );

vector<double> get_sig_positions(int forward[2000], 
	int reverse[2000], int N, PSSM * p, vector<double> background, double pv);


#endif