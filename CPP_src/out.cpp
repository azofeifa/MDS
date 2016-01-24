#include "out.h"
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
vector<double> sort(vector<double> X){
	if (X.size()<2){
		return X;
	}
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < X.size(); i++){
			if (X[i-1]  > X[i] ){
				changed 		= true;
				double copy 	= X[i-1];
				X[i-1] 	= X[i];
				X[i] 	= copy;
			}
		}
	}
	return X;
}


double get_pvalue(double obs, vector<double> null){
	double N 	= float(null.size());
	double S 	= 0.0;
	for (int i = 0 ; i < null.size();i++){
		S++;
		if (null[i] > obs){

			return S/N;
		}
	}
	return 1.0;
}




void write_out(string out_dir,map<int, vector< vector <double> >> collections, 
	vector<PSSM *>PSSMS ){
	ofstream FHW;
	FHW.open(out_dir + "motif_enrichment_statistics.tsv");
	typedef map<int, vector< vector <double> >>::iterator it_type;
	string name;
	for (it_type m = collections.begin(); m!=collections.end(); m++){
		name 	= PSSMS[m->first]->name;
		vector<double> observed_info 	= m->second[0];
		vector<double> null_stats 		= m->second[1];
		int step 						= null_stats.size()/4;
		vector<double> displacements 	= m->second[2];

		
		
		FHW<<name+"\t" ;
		string line 	= "";
		string pvals 	= "";
		int start, stop;
		if (observed_info[observed_info.size()-1]>0){
			for (int i =0; i < observed_info.size();i++){
				start 	= i*step, stop = (i+1)*step;
				line+=to_string(observed_info[i])+"\t";
				vector<double> current(null_stats.begin() + start, null_stats.begin()+stop);
				current 	= sort(current);
				pvals+=to_string(get_pvalue(observed_info[i], current) ) + "\t";
			}
			pvals=pvals.substr(0,pvals.size()-1);
			FHW<<line<<pvals<<endl;
		}else{
			FHW<<"0\t0\t0\t0\t0\t0\t0\t0"<<endl;	
		}
	}
	FHW<<"#null statistics"<<endl;
	for (it_type m = collections.begin(); m!=collections.end(); m++){
		name 	= PSSMS[m->first]->name;
		vector<double> observed_info 	= m->second[0];
		vector<double> null_stats 		= m->second[1];
		int step 						= null_stats.size()/4;
		vector<double> displacements 	= m->second[2];

		
		
		FHW<<name+"\t" ;
		string line 	= "";
		string pvals 	= "";
		int start, stop;
		if (observed_info[observed_info.size()-1]>0){
			for (int i =0; i < observed_info.size();i++){
				start 	= i*step, stop = (i+1)*step;
				vector<double> current(null_stats.begin() + start, null_stats.begin()+stop);
				current 	= sort(current);
				for (int j = 0; j < current.size();j++){
					line+=to_string(current[j])+",";
				}
				line=line.substr(0,line.size()-1) + "\t";
			}
			line=line.substr(0,line.size()-1) ;
			FHW<<line<<endl;
		}else{
			FHW<<"0\t0\t0\t0\t0\t0\t0\t0"<<endl;	
		}
	}
	FHW<<"#displacement data"<<endl;
	for (it_type m = collections.begin(); m!=collections.end(); m++){
		name 	= PSSMS[m->first]->name;
		FHW<<name + "\t";
		vector<double> displacements 	= m->second[2];
		string line = "";
		for (int i = 0; i < displacements.size();i++){
			line+=to_string(displacements[i]) + ",";

		}
		line=line.substr(0,line.size()-1);
		FHW<<line<<endl;

	}
	





}