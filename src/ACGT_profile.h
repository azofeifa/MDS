#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef ACGT_H
#define ACGT_H
#include <map>
#include <vector>
#include "load.h"

using namespace std;
map<int, double [2000][4]> get_average_ACGT_profile(map<string, vector<segment> > S, 
	vector<PSSM *> PSSMS, int pad ,map<int,  double  >   & NN,
	 map<int, double [2000][4]> & G,int );

void get_ACGT_profile_all(map<string, vector<segment> > , 
			  vector<vector<double>> &, int, int);
void get_1st_order_markov(map<string, vector<segment> >,vector<double ** > & , int );
#endif
