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

void sort2(vector<double> & X, vector<string> & Y, vector<double> & Z){
	int t 	= 0;
	if (X.size()!=2){

	
		bool changed 	= true;
		while (changed){
			changed = false;
			t+=1;
			for (int i = 1 ; i < X.size(); i++){
				if (X[i-1]  > X[i] ){
					changed 		= true;
					double copy 	= X[i-1];
					string cp_str 	= Y[i-1];
					double a_cp 	= Z[i-1];
					Y[i-1] 	= Y[i];
					X[i-1] 	= X[i];
					Z[i-1] 	= Z[i];
					X[i] 	= copy;
					Y[i] 	= cp_str;
					Z[i] 	= a_cp;
				}
			}
		}
	}
}


double get_pvalue(double obs, vector<double> null){
	double N 	= float(null.size());
	double S 	= 0.0;
	for (int i = 0 ; i < null.size();i++){
		if (null[i] > obs || null[i] == obs ){

			return S/N;
		}
		S++;
	}
	return 1.0;
}




void write_out(string out_dir,map<int, vector< vector <double> >> collections, 
	vector<PSSM *>PSSMS,string ID, vector<segment> intervals ){
	ofstream FHW;
	FHW.open(out_dir + ID+"-motif_enrichment_statistics.tsv");
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
	ofstream FHW2;
	FHW2.open(out_dir + ID+"-motif_hits.bed");
	for (int c = 0 ; c < intervals.size();c++ ){
		FHW2<<intervals[c].chrom<<"\t"<<to_string(intervals[c].start)<<"\t"<<to_string(intervals[c].stop)<<"\t";
		string line 	= "";
		vector<string> motif_names;
		vector<double> motif_positions_abs;
		vector<double> motif_positions_actual;

		typedef map<int, vector<double>>::iterator it_type;
		for (it_type cc=intervals[c].motif_positions.begin(); cc!=intervals[c].motif_positions.end(); cc++){
			for (int d = 0; d < cc->second.size(); d++){
				motif_names.push_back(PSSMS[cc->first]->name);
				motif_positions_abs.push_back(abs(cc->second[d]-1000));
				motif_positions_actual.push_back(cc->second[d]-1000);
				
			}
		}
		sort2(motif_positions_abs, motif_names,motif_positions_actual);
		for (int i = 0 ; i < motif_names.size();i++){
			line+=motif_names[i]+":"+to_string(int(motif_positions_actual[i])) + ",";
		}
		line=line.substr(0,line.size()-1);
		FHW2<<line<<endl;

	}



	





}