#include "simulate.h"
#include "scanner.h"
#include <omp.h>
#include <time.h>
#include <random>
#include <mpi.h>
#include <sys/time.h>

using namespace std;
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
}



mt19937  get_rand_seq(vector<int> & seq,  discrete_distribution<int> B, mt19937 g){
	for (int i = 0 ; i < seq.size(); i++){
		seq[i] 	= B(g);
	}
	return g;
}


int get_pvalue_llr(vector<int> seq, vector<double> background, double threshold, PSSM * p){
	int length 	= seq.size();
	double ll 	= 0;
	for ( int k=0; k < length; k++){
		ll+= p->frequency_table[k][seq[k]];
	}
	double pvalue 	= 1.0-p->get_pvalue(ll*2);
	if (pvalue < threshold){
		return 1;
	}
	return 0; 
}

vector<vector<int>> make_random_draws(	int sim_N, 
										vector<int> seqf, vector<int> seqr, 
										vector<discrete_distribution<int>> D_forward,
										vector<discrete_distribution<int>> D_reverse, mt19937 gen,
										vector<double> background, PSSM * p, double threshold){
	vector<vector<int>> DD(sim_N);
	int hit;
	#pragma omp parallel for
	for (int s = 0 ; s < sim_N; s++){//number of random samples
		vector<int> 	D;
		for (int d 	= 0 ; d < D_forward.size(); d++ ){
			gen 	= get_rand_seq(seqf,D_forward[d], gen );
			gen 	= get_rand_seq(seqr,D_reverse[d], gen );
			hit 	= get_pvalue_llr(seqf, background,threshold, p );
			if (hit){
				D.push_back(d);
			}
			hit 	= get_pvalue_llr(seqr, background,threshold, p );
			if (hit){
				D.push_back(d);
			}
		}
		DD[s]=D;
	}
	return DD;
}
vector<int> to_vector(int * array, int S){
	vector<int> A(S);
	for (int s = 0 ; s < S; s++){
		A[s] 	= array[s];
	}
	return A;
}

void send_out_null_displacement_data(vector<vector<int>> & D, int rank, 
										int nprocs, int ind_N){
	if (rank==0){
		for (int j = 1 ; j < nprocs; j++){
			for (int k = 0 ; k < ind_N; k++){
				int S;
				MPI_Recv(&S, 1, MPI_INT, j, k, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				int * A = new int[S];
				if (S){
					MPI_Recv(&A[0], S, MPI_INT, j, k+ind_N, MPI_COMM_WORLD,MPI_STATUS_IGNORE);					
				}
				D.push_back(to_vector(A,S));
			}
		}
	}else{
		for (int k = 0 ; k < ind_N; k++){
			int S 	= D[k].size();
			MPI_Ssend(&S, 1, MPI_INT, 0,k, MPI_COMM_WORLD);
			int * A = new int[S];
			if (S){
				copy(D[k].begin(), D[k].end(), A);
				MPI_Ssend(&A[0], S, MPI_INT, 0,k+ind_N, MPI_COMM_WORLD);
			}
						
		}
	}
}

string get_dots(int N){
	string line="";
	for (int i = 0 ; i < N; i++){
		line+=".";
	}
	return line;
}



void run_simulations(map<string, vector<segment>> intervals, 
		vector<PSSM *> P, int sim_N,
		vector<vector<double>> background_forward, vector<vector<double>>  background_reverse,
		vector<double> background, double threshold,int rank, int nprocs,Log_File * LG){


	for (int i =0 ; i < 4; i++){
		background[i] 	= log(background[i]);
	}
	

	vector<discrete_distribution<int>> D_forward;
	vector<discrete_distribution<int>> D_reverse;
	for (int i = 0 ; i < background_forward.size(); i++){
		discrete_distribution<int> Df {background_forward[i][0],background_forward[i][1],background_forward[i][2],background_forward[i][3]};
		discrete_distribution<int> Dr {background_reverse[i][0],background_reverse[i][1],background_reverse[i][2],background_reverse[i][3]};
		D_forward.push_back(Df);
		D_reverse.push_back(Dr);
	}

	random_device rd;
	// Initialize Mersenne Twister pseudo-random number generator
	mt19937 gen(rd()*rank);
	
	int ind_N 	= sim_N / nprocs;
	clock_t t;
	
	for (int p = 0 ; p < P.size(); p++){//iterate over every PSSM model
		int WN 	= max(int(44 - P[p]->name.size()), 1);	
		t = clock();
		LG->write(P[p]->name + get_dots(WN), 1);
		int S 	= P[p]->frequency_table.size(), hit=0;
		vector<int> seqf(S);
		vector<int> seqr(S);
		vector<vector<int>> DD 							= make_random_draws(ind_N, seqf, seqr, D_forward, D_reverse, gen, background, P[p], threshold );
		send_out_null_displacement_data(DD, rank, nprocs, ind_N);
		P[p]->null_displacements 	= DD;
		double wall1 = get_wall_time();
		t = clock() - t;

		LG->write("done: " + to_string(float(t)/CLOCKS_PER_SEC) + " seconds (" + to_string(p+1) + "/" + to_string(P.size())+")\n", 1);
	}
	//transform back to frequency space
	for (int p = 0; p < P.size(); p++){ //
		for (int i = 0; i < P[p]->frequency_table.size(); i++){
			for (int j = 0 ; j < 4;j++){
				P[p]->frequency_table[i][j]+=(background[j]);
				P[p]->frequency_table[i][j]=exp(P[p]->frequency_table[i][j]);
			}
		}
	}

}	
































