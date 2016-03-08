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
		double d 	= abs(displacements[i])-1000;
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

double get_min_2(vector<double> X){
	if (X.empty()){
		return 3000;
	}
	double MIN 	= abs(X[0]-1000);
	for (int xx = 1 ; xx < X.size(); xx++){
		if (abs(X[xx]-1000) < MIN) {
			MIN 	= abs(X[xx]-1000);
		}
	}
	return MIN;
}




void collect_sample_stats(map<string, vector<segment>> observed,
	 vector<PSSM *> P,  
	map<int, vector<double> > & observed_statistics,
	 map<int, vector<double>> & observed_displacements,map<int, map<int, double> > & observed_co_occurences,
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
	for (int p1 = 0 ; p1 < P.size(); p1++)	{
		for (int p2= 0; p2 < P.size(); p2++){
			observed_co_occurences[p1][p2] 	= 0;
		}
	}
	map<int,int> norms;
	for (int c = 0 ; c < combinded.size(); c++){
		for (it_type_4 p1 = combinded[c].motif_positions.begin(); p1!=combinded[c].motif_positions.end(); p1++){
			if ( get_min_2(p1->second) < 100 ){
				for (it_type_4 p2 = combinded[c].motif_positions.begin(); p2!=combinded[c].motif_positions.end(); p2++){
					if (get_min_2(p2->second) < 100){
						observed_co_occurences[p1->first][p2->first]+=1;
					}
				}
				if (norms.find(p1->first)==norms.end()){
					norms[p1->first]=0;
				}
				norms[p1->first]+=1;
			}

		}
	}
	typedef map<int, map<int, double> >::iterator it_type_7;
	typedef  map<int, double> ::iterator it_type_8;

	for (it_type_7 p1 = observed_co_occurences.begin(); p1!=observed_co_occurences.end(); p1++){
		for (it_type_8 p2 = p1->second.begin(); p2!=p1->second.end(); p2++){
			if (norms.find(p1->first)!=norms.end()){
				observed_co_occurences[p1->first][p2->first]/=	float(norms[p1->first]);
			}
		}		
	}	


	//this is PSSM_ID->to chuck->to displacements (to_chunk number gives you N_bidir)
	fill_displacements(observed_displacements, combinded);
	
	for (it_type_4 i = observed_displacements.begin(); i!=observed_displacements.end();i++){
		vector<double> current_stats;
		get_stats(i->second, current_stats);
		observed_statistics[i->first]=current_stats;
	}
	

}
