#define _GLIBCXX_USE_CXX11_ABI 0
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
void build_cdfs_PSSMs(PSSM *  , int, int , int , int, int, int,double, int);
double get_MD_score(vector<int>, int,bool,int);
vector<double> get_many_MD_scores(vector<int> , int, int);
#endif
