#include "out.h"
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
void collect_stats(map<string, vector<segment> > intervals, int p, double & mean, 
	double & var, double & se, double & NN){
	typedef map<string, vector<segment> >::iterator it_type;
	vector<double> distances;
	for (it_type c = intervals.begin(); c!=intervals.end(); c++){
		for (int i =0; i < c->second.size(); i++){
			for (int d = 0; d < c->second[i].motif_positions[p].size(); d++ ){
				distances.push_back(c->second[i].motif_positions[p][d]);
			}
		}
	}
	double S 	= 0;
	double SE 	= 0;
	double S2 	= 0;
	double N 	= 0;
	for (int d = 0; d < distances.size(); d++){
		S+=distances[d];
		S2+=pow(distances[d],2);
		if (abs(distances[d]) < 300){
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
	string out_dir, string job_ID , vector<PSSM *>PSSMS ){
	ofstream FHW;
	FHW.open(out_dir+ job_ID+ "_" + "observed_statistics.tsv" );
	FHW<<"#TF name\tSE ratio\tdisplacement mean\tdisplacement dispersion(0)\tnumber of motif locations\n";
	string motif 	= "";
	for (int p = 0; p < PSSMS.size(); p++){
		motif 		= PSSMS[p]->name;
		double mean =0, var=0, se=0, NN=0;
		collect_stats(intervals, p, mean, var, se,NN);
		FHW<<(motif+"\t" + to_string(se) + "\t" + to_string(mean) + "\t" + to_string(var)+ "\t" + to_string(NN)+ "\n");

	}



}