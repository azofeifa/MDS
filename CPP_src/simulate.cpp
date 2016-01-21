#include "simulate.h"

using namespace std;

void fill_in_sequence(int * X, vector<vector<double> > probs){


}


map<string, vector<vector<segment> >  > run_sims(map<string, vector<vector<double> > > GC, map<string, double> NN, int pad){ //need the ATGC profiles for each motif, return a list of sequences to run on 
	typedef map<string, vector<vector<double> > >::iterator it_type;
	for (it_type p = GC.begin(); p!= GC.end(); p++){
		int sequence[p->second.size()];
		int start 		= 0, stop 	= pad*2;
		string chrom  	= "chrS";
		fill_in_sequence(sequence, p->second);
	}


}

