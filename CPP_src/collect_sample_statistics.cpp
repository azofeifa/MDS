#include "collect_sample_statistics.h"
#include <cmath>
using namespace std;
void fill_displacements(map<int, vector<double>> & displacements, vector<segment> segments ){
	typedef map<int, vector<double>>::iterator it_type;
	for (int s = 0; s < segments.size(); s++){
		for (it_type p = segments[s].motif_positions.begin();p!=segments[s].motif_positions.end(); p++ ){
			if (p->second.size()){
				for (int i = 0; i < p->second.size(); i++){
					displacements[p->first].push_back(p->second[i]);
				}
			}else if (displacements.find(p->first)==displacements.end() ) {
				displacements[p->first]={};
			}
		}
	}
}

void get_stats(vector<double> displacements, vector<double> & stats    ){
	double mean=0, var=0, se=0, N=0;
	for (int i = 0 ; i < displacements.size(); i++ ){
		double d 	= displacements[i]-1000;
		mean+=abs(d);
		var+=pow(d,2);
		if (abs(d) <100){
			se+=1;
		}
		N+=1;
	}
	if (N>0){
		stats.push_back(se/N);
		stats.push_back(mean/N);
		stats.push_back(sqrt(var/N));
		stats.push_back(N);
	}else{
		stats.push_back(0);
		stats.push_back(0);
		stats.push_back(0);
		stats.push_back(N);	
	}
}




void collect_sample_stats(map<string, vector<segment>> observed,
	map<int, map<int, vector<segment> >> simulated, vector<PSSM *> P,  
	map<int, vector<double> > & observed_statistics,
	map<int, vector<vector<double> >> & observed_null_statistics, map<int, vector<double>> & observed_displacements,
	int rank ){

	//through back three things
	//(3) -> map[motif_ID] 	= binned histogram of displacements?
	vector<segment> combinded;
	typedef map<string, vector<segment>>::iterator it_type; //just go ahead and collapse down the segments from chromosomes
	typedef map<int, map<int, vector<segment> >>::iterator it_type_2; //just go ahead and collapse down the segments from chromosomes
	typedef  map<int, vector<segment> >::iterator it_type_3; //just go ahead and collapse down the segments from chromosomes
	typedef  map<int, vector<double> >::iterator it_type_4; //just go ahead and collapse down the segments from chromosomes
	typedef  map<int ,map<int, vector<double> >>::iterator it_type_5; //just go ahead and collapse down the segments from chromosomes
	for (it_type i = observed.begin(); i!=observed.end(); i++){
		for (int j = 0 ; j < i->second.size(); j++){
			combinded.push_back(i->second[j]);
		}
	}
	

	//this is PSSM_ID->to chuck->to displacements (to_chunk number gives you N_bidir)
	map<int, map<int, vector<double>>> simulated_displacements; //
	fill_displacements(observed_displacements, combinded);
	for (it_type_2 i = simulated.begin(); i!=simulated.end(); i++){ //iterating over PSSM_ID
		for (it_type_3 j = i->second.begin(); j!=i->second.end(); j++){ //iterating over chunk, i->second.size()==sim_N
			bool FOUND=false;
			for (int s 	= 0; s<j->second.size(); s++){//
				if (!j->second[s].motif_positions[i->first].empty()){
					for (int z=0; z< j->second[s].motif_positions[i->first].size();z++){
						FOUND=true;
						simulated_displacements[i->first][j->first].push_back(j->second[s].motif_positions[i->first][z]);
					}
				}
			}
			if (!FOUND){
				simulated_displacements[i->first][j->first]={};			
			}
		}
	}
	for (it_type_4 i = observed_displacements.begin(); i!=observed_displacements.end();i++){
		vector<double> current_stats;
		get_stats(i->second, current_stats);
		observed_statistics[i->first]=current_stats;
	}
	for (it_type_5 i=simulated_displacements.begin(); i!=simulated_displacements.end();i++){
		for (it_type_4 j =i->second.begin(); j!=i->second.end(); j++){
			vector<double> current_stats;
			get_stats(j->second, current_stats);
			observed_null_statistics[i->first].push_back(current_stats)	;
		}
	}

}
