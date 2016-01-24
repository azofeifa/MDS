#ifndef out_H
#define out_H
#include <map>
#include <vector>
#include "load.h"
using namespace std;
void write_out(string out_dir,
	map<int, vector< vector <double> >>, 
	vector<PSSM *>PSSMS);
#endif