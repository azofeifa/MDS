#include "read_in_parameters.h"
#include "load.h"
#include "get_motif_pvalues.h"
#include "scanner.h"
#include "out.h"
#include "ACGT_profile.h"
#include "simulate.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include "MPI_comm.h"
#include "collect_sample_statistics.h"
using namespace std;

int main(int argc,char* argv[]){
	MPI::Init(argc, argv);
	clock_t t1,t2;
	int nprocs		= MPI::COMM_WORLD.Get_size();
	int rank 		= MPI::COMM_WORLD.Get_rank();
    int threads  	= omp_get_max_threads();

	params * P = new params();
	fill_in_options(argv, P, 0);
	if (P->EXIT){
		printf("exiting...\n");
		MPI::Finalize();
		return 0;
	}


	
	if (P->module=="DB"){
		//============================================================
		//necessary user input parameters for DB module
		string fasta_file 			= P->p["-fasta"];
		string bed_file 			= P->p["-bed"];
		string OUT 				= P->p["-o"];
		string PSSM_DB 				= P->p["-DB"];
		string job_ID 				= P->p["-ID"];
		double window 				= 1000;
		double pv 					= stof(P->p["-pv"]);
		int test 					= 1;
		int sim_N 					= stoi(P->p["-sim_N"]);
		int bins 					= stoi(P->p["-br"]);
		int interval_size 			= 0;

		vector<PSSM *> PSSMS;
		//============================================================


		//============================================================
		//....1.... Load PSSM Database
		printf("loading PSSM DB.............................");
		cout.flush();
		PSSMS 		= load_PSSM_DB_new(PSSM_DB, test );
		printf("done\n");
		//============================================================
		//....2.... Compute P-values, approximation error ~ 1.0 / bins
		printf("computing p-values..........................");
		cout.flush();
		vector<double> background 	= {0.25,0.25,0.25,0.25};
		DP_pvalues(PSSMS,bins, background, true);
		printf("done\n");
		

		//============================================================
		//....3.... load intervals from user 
		printf("loading intervals...........................");
		cout.flush();
		map<string, vector<segment>> intervals 	= load_bed_file(bed_file, window,interval_size); 
		printf("done\n");
		
		//============================================================
		//....4.... insert the fasta sequnece into the provided intervals
		printf("inserting fasta.............................");
		cout.flush();
		intervals 								= insert_fasta_sequence(fasta_file, intervals,test);
		printf("done\n");
		
		//============================================================
		//....5.... get generic GC content across intervals for background model
		vector<vector<double>> background_forward, background_reverse;
		printf("computing ACGT profile......................");
		cout.flush();
		get_ACGT_profile_all(intervals, 
			background_forward,background_reverse, rank);
		printf("done\n");

		//============================================================
		//....6.... perform simulations and get random MD scores
		printf("running simulations for null................");
		cout.flush();
		run_sim3(intervals,PSSMS, sim_N,background_forward, background_reverse,background, pv);
		printf("done\n");
		//============================================================
		//....6....write out to DB file
		printf("writing out simulations.....................");
		cout.flush();
		write_out_null_stats( PSSMS,OUT,  P, background);
		printf("done\n");
	}else if (P->module=="EVAL"){
		//============================================================
		//necessary user input parameters for DB module
		string fasta_file 			= P->p["-fasta"];
		string bed_file 			= P->p["-bed"];
		string DB_file 				= P->p["-DB"];
		string OUT 					= P->p["-o"]; 
		int BSN 					= stoi(P->p["-bsn"]);
		int test 					= 1;
		double window 				= 1000;
		int interval_size 			= 0;
		//the rest of the parameters will be in the DB file
		//============================================================


		//============================================================
		//....1.... Load PSSM Database
		printf("loading PSSM DB.............................");
		cout.flush();
		vector<double> background;
		vector<PSSM *> PSSMS 					= load_personal_DB_file(DB_file, P,background);
		printf("done\n");
		int bins 		= stoi(P->p["-bins"]);
		double pv 		= stod(P->p["-pv"]);

		//============================================================
		//....2.... Load user provided intervals
		printf("loading intervals...........................");
		cout.flush();
		map<string, vector<segment>> intervals 	= load_bed_file(bed_file, window,interval_size); 
		printf("done\n");
		
		//============================================================
		//....3.... Computed LLR distribution for motifs
		printf("computing p-values..........................");
		cout.flush();
		DP_pvalues(PSSMS,bins, background,false);
		printf("done\n");


		
		//============================================================
		//....5.... insert the fasta sequnece into the provided intervals
		printf("inserting fasta.............................");
		cout.flush();
		intervals 								= insert_fasta_sequence(fasta_file, intervals,0);
		printf("done\n");
		
		//============================================================
		//....6.... scan the provided intervals
		printf("scanning intervals..........................");
		cout.flush();
		run_accross(intervals , PSSMS,  background, pv, interval_size, BSN);
		printf("done\n");
		
		//============================================================
		//....7.... assess significance
		write_out_stats(PSSMS, OUT, P);






	}
	



	MPI::Finalize();


	return 0;
}