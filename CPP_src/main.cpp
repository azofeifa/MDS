#include "read_in_parameters.h"
#include "load.h"
#include "get_motif_pvalues.h"
#include "scanner.h"
#include "out.h"
#include "ACGT_profile.h"
#include "simulate.h"
#include "error_stdo_logging.h"
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
		string OUT 					= P->p["-o"];
		string PSSM_DB 				= P->p["-DB"];
		int job_ID 					= 1;
		double window 				= 1000;
		double pv 					= stof(P->p["-pv"]);
		int test 					= 0;
		int sim_N 					= stoi(P->p["-sim_N"]);
		int bins 					= stoi(P->p["-br"]);
		int interval_size 			= 0;
		int verbose 				= 1;
		vector<PSSM *> PSSMS;
		//============================================================

	    Log_File * LG 	= new  Log_File(rank, job_ID, P->p["-ID"], P->p["-log_out"]);

	    LG->write(P->get_header(), verbose);
		//============================================================
		//....1.... Load PSSM Database

		LG->write("loading PSSM DB.............................", verbose);
		PSSMS 		= load_PSSM_DB_new(PSSM_DB, test );
		LG->write("done\n", verbose);
		//============================================================
		//....2.... Compute P-values, approximation error ~ 1.0 / bins
		LG->write("computing p-values..........................",verbose);
		vector<double> background 	= {0.25,0.25,0.25,0.25};
		DP_pvalues(PSSMS,bins, background, true);
		LG->write("done\n", verbose);
		

		//============================================================
		//....3.... load intervals from user 
		LG->write("loading intervals...........................", verbose);
		map<string, vector<segment>> intervals 	= load_bed_file(bed_file, window,interval_size); 
		LG->write("done\n", verbose);
		
		//============================================================
		//....4.... insert the fasta sequnece into the provided intervals
		LG->write("inserting fasta.............................", verbose);
		intervals 								= insert_fasta_sequence(fasta_file, intervals,test);
		LG->write("done\n", verbose);
		
		//============================================================
		//....5.... get generic GC content across intervals for background model
		vector<vector<double>> background_forward, background_reverse;
		LG->write("computing ACGT profile......................", verbose);
		get_ACGT_profile_all(intervals, 
			background_forward,background_reverse, rank);
		LG->write("done\n",verbose);

		//============================================================
		//....6.... perform simulations and get random MD scores
		LG->write("\n     running simulations for null model\n\n",verbose);
		run_simulations(intervals,PSSMS, sim_N,background_forward, background_reverse,background, pv, rank, nprocs, LG);
		//============================================================
		//....6....write out to DB file
		if (rank==0){
			LG->write("writing out simulations.....................",verbose);
			write_out_null_stats( PSSMS,OUT,  P, background);
			LG->write("done\n",verbose);
		}
		if (rank==0){
			collect_all_tmp_files(P->p["-log_out"], P->p["-ID"], nprocs, job_ID);
		}

	}else if (P->module=="EVAL"){
		//============================================================
		//necessary user input parameters for DB module
		string fasta_file 			= P->p["-fasta"];
		string bed_file 			= P->p["-bed"];
		string DB_file 				= P->p["-DB"];
		string OUT 					= P->p["-o"]; 
		int BSN 					= stoi(P->p["-bsn"]);
		int test 					= 1;
		int job_ID 					= 1;
		double window 				= 1000;
		int verbose 				= 1;
		int interval_size 			= 0;
		//the rest of the parameters will be in the DB file
		//============================================================

		Log_File * LG 	= new  Log_File(rank, job_ID, P->p["-ID"], P->p["-log_out"]);
		LG->write(P->get_header(), verbose);
		

		//============================================================
		//....1.... Load PSSM Database
		LG->write("loading PSSM DB.............................", verbose);
		vector<double> background;
		vector<PSSM *> PSSMS 					= load_personal_DB_file(DB_file, P,background);
		if (PSSMS.empty()){
			if (rank==0){
				collect_all_tmp_files(P->p["-log_out"], P->p["-ID"], nprocs, job_ID);
			}
			LG->write("exiting...\n", verbose);
			MPI::Finalize();
			return 0;
		}
		LG->write("done\n",verbose);
		
		int bins 		= stoi(P->p["-bins"]);
		double pv 		= stod(P->p["-pv"]);

		//============================================================
		//....2.... Load user provided intervals
		LG->write("loading intervals...........................", verbose);
		map<string, vector<segment>> intervals 	= load_bed_file(bed_file, window,interval_size); 
		LG->write("done\n", verbose);
		
		//============================================================
		//....3.... Computed LLR distribution for motifs
		LG->write("computing p-values..........................", verbose);
		DP_pvalues(PSSMS,bins, background,false);
		LG->write("done\n", verbose);


		
		//============================================================
		//....5.... insert the fasta sequnece into the provided intervals
		LG->write("inserting fasta.............................", verbose);
		intervals 								= insert_fasta_sequence(fasta_file, intervals,test);
		LG->write("done\n",verbose);
		
		//============================================================
		//....6.... scan the provided intervals
		LG->write("\n         scanning intervals\n\n",verbose);
		
		scan_intervals(intervals, PSSMS, background, 
 								pv,  interval_size,   BSN, 
 								rank,   nprocs,   LG);

		
		//============================================================
		//....7.... assess significance
		if (rank==0){
			write_out_stats(PSSMS, OUT, P);
			collect_all_tmp_files(P->p["-log_out"], P->p["-ID"], nprocs, job_ID);
		}
	}
	


	MPI::Finalize();


	return 0;
}