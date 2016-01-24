#include "simulate.h"
#include <omp.h>
#include <random>
using namespace std;



map<int, map<int, vector<segment> >> run_sims(map<int, double [2000][4]> GC, 
	map<int, double> NN, vector<PSSM *> P,int sim_N, int rank ){ 
	default_random_engine generator;
	map<int, int>flip;
	flip[0] 	=3, flip[3] = 0, flip[1]=2, flip[2]=1;

	typedef map<int, double [2000][4]>::iterator it_type;
	map<int, map<int, vector<segment> >> S;//PSSM key -> chunk number -> spec simulations
	for (it_type g = GC.begin(); g!=GC.end(); g++){
		//make rand generator
		vector<discrete_distribution<int> > dists(2000);
		for (int f = 0; f < 2000; f++){
			discrete_distribution<int> distribution{GC[g->first][f][0],GC[g->first][f][1],
				GC[g->first][f][2],GC[g->first][f][3]  };
			dists[f] 	= distribution ;
		}
		for (int s = 0; s <  sim_N; s++){//simulate sim_N groups
			vector<segment> currents(int(NN[g->first]));
			int NNN 		= int(NN[g->first]);
			#pragma omp parallel for
			for (int i = 0; i <NNN; i++){ //within each group number of times that motif was found in a bidirectional
				segment current_S 	= segment();
				for (int j=0; j < 2000; j++){
					current_S.forward[j] 	= dists[j](generator);
					current_S.reverse[j] 	= flip[current_S.forward[j]];
				}
				currents[i] 			= current_S;
			}
			if (!NNN){
				S[g->first][s] 	= {};
			
			}else{
				S[g->first][s] 	= currents;
			}
		}
	}
	return S;


}

