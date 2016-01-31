#include "simulate.h"
#include "scanner.h"
#include <omp.h>
#include <random>
using namespace std;

int PSSM_index_2(int i, vector<PSSM *> P){
	for (int j = 0 ; j < P.size(); j++){
		if (P[j]->ID==i){
			return j;
		}
	}
	printf("what???\n");
	return 0;
}

vector<double> get_stats(vector<double> displacements  ){
	double mean=0, var=0, se=0, N=0;
	vector<double > stats; 
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
	return stats;

}

void run_sims(map<int, double [2000][4]> GC, 
	map<int, double> NN, vector<PSSM *> P,int sim_N, int rank, 
	vector<double> background, double pv, 
	map<int, vector<vector<double> >> & observed_null_statistics){ 
	
	default_random_engine generator;
	map<int, int>flip;
	flip[0] 	=3, flip[3] = 0, flip[1]=2, flip[2]=1;

	typedef map<int, double [2000][4]>::iterator it_type;
	map<int, map<int, vector<segment> >> S;//PSSM key -> chunk number -> spec simulations
	for (int b = 0 ; b < background.size(); b++){
		background[b] 	= log(background[b]);
	}
		

	for (it_type g = GC.begin(); g!=GC.end(); g++){
		//make rand generator
		vector<discrete_distribution<int> > dists(2000);
		#pragma omp parallel for
		for (int f = 0; f < 2000; f++){
			discrete_distribution<int> distribution{GC[g->first][f][0],GC[g->first][f][1],
				GC[g->first][f][2],GC[g->first][f][3]  };
			dists[f] 	= distribution ;
		}

		for (int s = 0; s <  sim_N; s++){//simulate sim_N groups
			vector<int *> forwards(int(NN[g->first]));
			vector<int *> reverses(int(NN[g->first]));
			
			int NNN 		= int(NN[g->first]);
			#pragma omp parallel for
			for (int i = 0; i <NNN; i++){ //within each group number of times that motif was found in a bidirectional
				int * forward 	= new int[2000];
				int * reverse 	= new int[2000];
				for (int j=0; j < 2000; j++){
					forward[j] 	= dists[j](generator);
					reverse[j] 	= flip[forward[j]];
				}
				forwards[i] 			= forward;
				reverses[i] 			= reverse;
				
			}
			vector<vector<double>> positions(NNN);

			#pragma omp parallel for
			for (int i = 0 ; i < NNN; i++){
				PSSM * p 	= P[PSSM_index_2(g->first, P) ];
				positions[i] 	= get_sig_positions(forwards[i], 
					reverses[i], 2000, p, background, pv);				
			}
			vector<double> final_positions;
			for (int i = 0 ; i < NNN;i++){
				for (int c = 0; c < positions[i].size(); c++){
					final_positions.push_back(positions[i][c]);
				}
				delete forwards[i];
				delete reverses[i];
			}
			vector<double> stats 		= get_stats(final_positions);
			observed_null_statistics[g->first].push_back(stats);
		}
	}

}

void run_sims2(map<string, vector<segment>> intervals, vector<PSSM *> P,int sim_N, int rank, 
	vector<double> background, double pv, 
	map<int, vector<vector<double> >> & observed_null_statistics){
	for (int b = 0 ; b < background.size(); b++){
		background[b] 	= log(background[b]);
	}
	
	vector<segment> data;
	typedef map<string, vector<segment>>::iterator it_type;
	typedef vector<segment>::iterator it_type_2;
	typedef vector<PSSM *>::iterator it_type_3;

	for (it_type i = intervals.begin(); i!=intervals.end(); i++){
		for (it_type_2 j = i->second.begin(); j!=i->second.end(); j++){
			data.push_back((*j)) ;
		}
	}

	int N 						= data.size();
	int ** forward_matrix 	= new int*[N];
	int ** reverse_matrix 	= new int*[N];
	for (int j = 0 ; j < N; j++){
		forward_matrix[j] 		= new int[2000];
		reverse_matrix[j] 		= new int[2000];
	}
	default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,N-1);
	for (int s = 0; s < sim_N; s++){
		//want to make new forward and reverse, bootstrapping? random shuffling
		for (int j = 0 ; j < N; j++){ //iterate over sequences
			for (int k = 0; k < 2000; k++){ //iterate over positions
				//pick a random integer between 0,N-1
				int l 					= distribution(generator);
				forward_matrix[j][k] 	= data[l].forward[k];
				reverse_matrix[j][k] 	= data[l].reverse[k];
			}
		}
		//okay now scan over forward and reverse _matrix over all PSSMS
		for (it_type_3 p = P.begin(); p!=P.end(); p++){

			vector<vector<double>> positions(N);

			#pragma omp parallel for
			for (int i = 0 ; i < N; i++){
				positions[i] 	= get_sig_positions(forward_matrix[i], 
					reverse_matrix[i], 2000, *p, background, pv);		

			}
			vector<double> final_positions;
			for (int i = 0 ; i < N;i++){
				for (int c = 0; c < positions[i].size(); c++){
					final_positions.push_back(positions[i][c]);
				}
			}
			vector<double> stats 		= get_stats(final_positions);
			observed_null_statistics[(*p)->ID].push_back(stats);
		}

	}




}




