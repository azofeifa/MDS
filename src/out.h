#define _GLIBCXX_USE_CXX11_ABI 0

#ifndef out_H
#define out_H
#include <map>
#include <vector>
#include "load.h"
using namespace std;
void write_out(string out_dir,
	map<int, vector< vector <double> >>, 
	vector<PSSM *>PSSMS,string,
	map<int, segment> );
void write_out_2(string out_dir,string ID, vector<PSSM *>,  map<int, vector<double> > , 
	map<int, vector<double>>   , map<int, map<int, double> > ,
	map<int, vector<vector<double> >> , map<int, map<int, vector<double> > >  );
void write_out_3(string ,string,map<int, string>, map<int, vector<double> >,vector<vector<double>>,vector<vector<double>>  );

#endif