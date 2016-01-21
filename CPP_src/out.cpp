#include "out.h"
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
void collect_stats(map<string, vector<segment> > intervals, int p, double & mean, 
	double & var, double & se, double & NN, double & NNN, int pad){
	typedef map<string, vector<segment> >::iterator it_type;
	vector<double> distances;
	for (it_type c = intervals.begin(); c!=intervals.end(); c++){
		for (int i =0; i < c->second.size(); i++){
			if (c->second[i].motif_positions[p].size()){
				NNN+=1.0;
			}
			for (int d = 0; d < c->second[i].motif_positions[p].size(); d++ ){
				distances.push_back(c->second[i].motif_positions[p][d]);
			}
		}
	}
	double S 	= 0;
	double SE 	= 0;
	double S2 	= 0;
	double N 	= 0;
	double d 	= 0.0;
	for (int j = 0; j< distances.size(); j++){
		d 		= distances[j] - pad;

		S+=d;
		S2+=pow(d,2);
		if (abs(d) < 300){
			SE++;
		}
		N++;
	}
	if (N==0){
		se = 0.0, mean 	= 0.0, var = 0.0, NN =0.0;

	}else{
		se = SE / N, mean 	= S/N, var = sqrt(S2 / N), NN = N ;

	}



}



void write_observation_stats(map<string, vector<segment> > intervals, 
	string out_dir, string job_ID , vector<PSSM *>PSSMS,map<string, vector<vector<double> > > G, int pad){
	ofstream FHW;
	FHW.open(out_dir+ job_ID+ "_" + "observed_statistics.tsv" );

	FHW<<"#TF name\tSE ratio\tdisplacement mean\tdisplacement dispersion(0)\tnumber of motif locations\tNumber of Bidirections with motif\n";
	string motif 	= "";



	for (int p = 0; p < PSSMS.size(); p++){
		motif 		= PSSMS[p]->name;
		double mean =0, var=0, se=0, NN=0, NNN=0;
		collect_stats(intervals, p, mean, var, se,NN, NNN, pad);
		FHW<<(motif+"\t" + to_string(se) + "\t" + to_string(mean) + "\t" + to_string(var)+ "\t" + to_string(NN)+ "\t" + to_string(NNN)+ "\n");

	}

	FHW<<"#ACGT Frequency Distributions per motif\n";
	typedef map<string, vector<vector<double> > >::iterator it_type;
	for (it_type m = G.begin(); m!=G.end(); m++){
		FHW<<m->first+"\t";
		string line="";
		string temp="";
		for (int i = 0; i < m->second.size(); i++){
			temp = "";
			for (int j = 0; j < m->second[i].size(); j++){
				temp+=to_string(m->second[i][j])+",";
			}
			line+=temp.substr(0,temp.size()-1) + "\t";
		}
		line=line.substr(0,line.size()-1);
		FHW<<line<<endl;

	}


}