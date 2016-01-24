#ifndef MPI_H
#define MPI_H
#include "load.h"
#include <vector>
using namespace std;
void test_send_PSSMS(int, int, vector<PSSM *>,
vector<vector<vector<double>>> & , vector<int>  &, vector<int> &  );

map<int, vector< vector <double> >> collect_PSSM_hits(int  , int  , 
	map<string, vector<segment>>  ,
	map<int, vector<double> >  ,
	map<int, vector<vector<double> >>  ,
	map<int, vector<double>>  	   );





#endif