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


mt19937  get_rand_seq_2( int forward[2000], int reverse[2000],  
							vector<discrete_distribution<int>> D_forward, 
							vector<discrete_distribution<int>> D_reverse,
							mt19937 g){
	for (int i = 0 ; i < 2000; i++){
		forward[i] 	= D_forward[i](g),reverse[i] 	= D_reverse[i](g);
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
										vector<double> background, PSSM * p, double threshold){
	vector<vector<int>> DD(sim_N);
	random_device rd;
	map<int, int> G;
	G[0] 	=3,G[3] 	=0,G[1] 	=2,G[2] 	=1;
	#pragma omp parallel for
	for (int s = 0 ; s < sim_N; s++){//number of random samples
		mt19937 gen(rd()*s);
		int forward[2000], reverse[2000];
		for (int i = 0 ; i < 2000; i++){
			forward[i] 	= D_forward[i](gen);
			reverse[i] 	= G[forward[i]];
		}
		DD[s] 		= get_sig_positions(forward, reverse, 2000, p, threshold);
	}
	return DD;
}

vector<vector<int>> make_random_draws_from_markov(int sim_N, PSSM * p, 
				double threshold,vector<vector<discrete_distribution<int>>> D_forward_MM){
	vector<vector<int>> DD(sim_N);
	random_device rd;
	discrete_distribution<int> U 	= {0.25,0.25,0.25,0.25};
	map<int, int> G;
	G[0] 	=3,G[3] 	=0,G[1] 	=2,G[2] 	=1;
	#pragma omp parallel for
	for (int s = 0 ; s < sim_N; s++){//number of random samples
		mt19937 gen(rd()*s);
		int forward[2000], reverse[2000];
		forward[0] 	= U(gen);
		reverse[0] 	= G[forward[0]];
		for (int i = 1 ; i < 2000; i++){
			forward[i] 	= D_forward_MM[i-1][forward[i-1]](gen);
			reverse[i] 	= G[forward[i]];
		}
		DD[s] 		= get_sig_positions(forward, reverse, 2000, p, threshold);
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

void send_out_null_displacement_data(vector<vector<int>> & D, int rank, int nprocs, int ind_N){
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
		vector<vector<double>> background_D,
		vector<double> background, double threshold,int rank, int nprocs,Log_File * LG, int order,
		vector<double ** > MM, int TSS){


	vector<discrete_distribution<int>> D_forward;
	vector<discrete_distribution<int>> D_reverse;

	vector<vector<discrete_distribution<int>>> D_forward_MM(2000);

	if (order==0){
		for (int i = 0 ; i < background_D.size(); i++){
			discrete_distribution<int> Df {background_D[i][0],background_D[i][1],background_D[i][2],background_D[i][3]};
			D_forward.push_back(Df);
		}
	}else{
		for (int i = 0 ; i < MM.size(); i++){
			vector<discrete_distribution<int>> current(4);
			for (int u = 0 ; u < 4; u++){
				discrete_distribution<int> Df {MM[i][u][0],MM[i][u][1],MM[i][u][2],MM[i][u][3]};
				current[u] 	= Df;
			}
			D_forward_MM[i] 	= current;

		}
	}
	
	int ind_N 	= sim_N / nprocs;
	clock_t t;
	int threads 	= omp_get_max_threads();
	for (int p = 0 ; p < P.size(); p++){//iterate over every PSSM model
		int WN 	= max(int(44 - P[p]->name.size()), 1);	
		t = clock();
		LG->write(P[p]->name + get_dots(WN), 1);
		int S 	= P[p]->frequency_table.size(), hit=0;
		vector<int> seqf(S);
		vector<int> seqr(S);
		vector<vector<int>> DD ;
		if (order==0){
			DD 							= make_random_draws(ind_N, seqf, seqr, 
				D_forward, background, P[p], threshold );
		}else{
			DD 							= make_random_draws_from_markov(ind_N,P[p], threshold, D_forward_MM);
		}
		send_out_null_displacement_data(DD, rank, nprocs, ind_N);
		if (TSS){
			P[p]->null_displacements 		= DD;
		}else{
			P[p]->null_displacements_non 	= DD;			
		}

		double wall1 = get_wall_time();
		t = clock() - t;

		LG->write("done: " + to_string(float(t)/(CLOCKS_PER_SEC)) + " seconds (" + to_string(p+1) + "/" + to_string(P.size())+")\n", 1);
	}
}	
































