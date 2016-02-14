#include "simulate.h"
#include "scanner.h"
#include <omp.h>
#include <time.h>
#include <random>
#include <mpi.h>

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

vector<double> get_stats_2(vector<double> displacements  ){
	double mean=0, var=0, se=0, N=0;
	vector<double > STATS; 
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
		STATS.push_back(se/N);
		STATS.push_back(mean/N);
		STATS.push_back(sqrt(var/N));
		STATS.push_back(N);
	}else{
		STATS.push_back(0);
		STATS.push_back(0);
		STATS.push_back(0);
		STATS.push_back(N);	
	}
	return STATS;

}


int ** allocate_2D_array(int I,int J, int K){
	int ** X 	= new int*[I];
	for (int i = 0 ; i < I; i++){
		X[i] 		= new int[J*K];
		
	}
	return X;
	
}


double ** allocate_2D_array_double(int I,int J, int K){
	double ** X 	= new double*[I];
	for (int i = 0 ; i < I; i++){
		X[i] 		= new double[J*K];
	}
	return X;
	
}



double get_min(vector<double> X){
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


void FREE_int(int ** X, int I ){
	for (int i = 0 ; i < I; i++){
		delete [] X[i];
	}
	delete []  X;
	
}

void FREE_double(double ** X, int I ){
	for (int i = 0 ; i < I; i++){
		delete [] X[i];
	}
	delete []  X;
	
}


void run_sims2(map<string, vector<segment>> intervals, vector<PSSM *> P,int sim_N, int rank, int nprocs,
	vector<double> background, double pv, 
	map<int, vector<vector<double> >> & observed_null_statistics,map<int, map<int, vector<int> >> & null_co_occur){



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
	std::uniform_int_distribution<int> distribution(0,N-1);
	typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
	uint32_t seed_val 	= time(NULL)*rank;           // populate somehow

	MyRNG rng;                   // e.g. keep one global instance (per thread)
	rng.seed(seed_val);
	//need to make co occurrence 3D matrix
	
	int PN 	= P.size();
	//want to send out start and stops of PSSMS to scan ....
	int count 	= sim_N / nprocs;
	int ** all_co 	= allocate_2D_array(count, P.size(), P.size());
	double ** stats = allocate_2D_array_double(count,P.size(),  4 );
	for (int s = 0; s < count; s++){
		//want to make new forward and reverse, bootstrapping? random shuffling
		for (int j = 0 ; j < N; j++){ //iterate over sequences
			for (int k = 0; k < 2000; k++){ //iterate over positions
				//pick a random integer between 0,N-1
				int l 					= distribution(rng);
			
				forward_matrix[j][k] 	= data[l].forward[k];
				reverse_matrix[j][k] 	= data[l].reverse[k];

				forward_matrix[l][k] 	= data[j].forward[k];
				reverse_matrix[l][k] 	= data[j].reverse[k];

			}

		}

		//okay now scan over forward and reverse _matrix over all PSSMS
		map<int, vector<double> > null_positions;
		vector<vector<vector<double>>> all_hits(N);
		vector<vector<double>> positions_by_motif(P.size());
		for (int i = 0 ; i < N ; i++){

			vector<vector<double>> hits(PN);
			#pragma omp parallel for
			for (int p = 0; p<PN; p++){
				vector<double> positions 	= get_sig_positions(forward_matrix[i], 
					reverse_matrix[i], 2000, P[p], background, pv);		
				hits[p] 		= positions;

			}			
			all_hits[i] 	= hits;
		}
		for (int i = 0 ; i < N; i++){
			for (int j = 0; j < P.size();j++ ){
	
				if ( get_min(all_hits[i][j])<100  ){
					for (int k = 0; k < P.size();k++ ){
						if (  get_min(all_hits[i][k])<100  ){
							all_co[s][j*P.size() + k]+=1;
						}
					}
				}
				positions_by_motif[j].insert(positions_by_motif[j].end(), all_hits[i][j].begin(),all_hits[i][j].end() );
			}
		}
		for (int i = 0 ; i < positions_by_motif.size(); i++){
			vector<double> cs 	= get_stats_2(positions_by_motif[i]);
			stats[s][i*4 + 0] 	=	 cs[0],stats[s][i*4 + 1] 	= cs[1],stats[s][i*4 + 2] 	= cs[2],stats[s][i*4 + 3] 	= cs[3];
		}
	}

	//want to send the all_co and the positions_by_motifs to rank
	vector<int **> all_co_final;
	vector<double **> all_stats_final;
	
	if (rank==0){
		
		all_co_final.push_back(all_co);
		all_stats_final.push_back(stats);
		for (int j = 1; j < nprocs; j++){
			int ** RR 					= allocate_2D_array(count, P.size(), P.size());
			double ** current_stats 	= allocate_2D_array_double( count, P.size(), 4);
			for (int u = 0 ; u < count; u++){
				MPI_Recv(&RR[u][0], P.size()*P.size(), MPI_INT, j, u, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&current_stats[u][0], 4*P.size(), MPI_DOUBLE, j, u+count, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			all_co_final.push_back(RR);
			all_stats_final.push_back(current_stats);
		}
	}else{
		for (int u = 0 ; u < count; u++){
			MPI_Send(&all_co[u][0],  P.size()*P.size(), MPI_INT, 0,u, MPI_COMM_WORLD);
			MPI_Send(&stats[u][0],4*P.size(), MPI_DOUBLE, 0,u+count, MPI_COMM_WORLD);
		}
	}
	for (int p = 0; p < P.size(); p++){
		vector<double> C1,C2,C3,C4;
		observed_null_statistics[P[p]->ID].push_back(C1);
		observed_null_statistics[P[p]->ID].push_back(C2);
		observed_null_statistics[P[p]->ID].push_back(C3);
		observed_null_statistics[P[p]->ID].push_back(C4);	
	}


	if (rank==0){
		for (int ll = 0 ; ll < all_co_final.size(); ll++){
			for (int p = 0; p < P.size(); p++){
				for (int u = 0; u < count;u++){
					observed_null_statistics[P[p]->ID][0].push_back(all_stats_final[ll][u][p*4+0]);
					observed_null_statistics[P[p]->ID][1].push_back(all_stats_final[ll][u][p*4+1]);
					observed_null_statistics[P[p]->ID][2].push_back(all_stats_final[ll][u][p*4+2]);
					observed_null_statistics[P[p]->ID][3].push_back(all_stats_final[ll][u][p*4+3]);
				}
			}
			for (int p1 = 0; p1 < P.size(); p1++){
				for (int p2 = 0; p2 < P.size(); p2++){
					for (int u = 0; u < count; u++){
						null_co_occur[P[p1]->ID][P[p2]->ID].push_back(all_co_final[ll][u][p1*P.size() + p2]);
					}
				}
			}
			FREE_double(all_stats_final[ll], count);
			FREE_int(all_co_final[ll], count);

		}

	}
	else{
		FREE_int(all_co, count);
		FREE_double(stats, count);
	}


}




