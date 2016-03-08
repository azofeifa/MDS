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
	if (P->check()){
		printf("exiting...\n");
		return 0;
	}
	vector<PSSM *> PSSMS;
	vector<vector<vector<double>>> streamed_PSSMS;
	vector<int>  PSSM_IDS;
	vector<int>  PSSM_N;


	//============================================================
	//necessary user input parameters
	string fasta_file 			= P->p["-fasta"];
	string bed_file 			= P->p["-bed"];
	string out_dir 				= P->p["-o"];
	string log_out 				= P->p["-log_out"];
	string PSSM_DB 				= P->p["-DB"];
	string job_ID 				= P->p["-ID"];
	string bed_out 				= P->p["-bed_out"];
	int site_br 				= stoi(P->p["-site_br"]);
	int simN 					= stoi(P->p["-simN"]);
	int pad 					= 1000;
	int bins 					= stoi(P->p["-br"]);
	double pv 					= stod(P->p["-pv"]);
	//============================================================

	ofstream FHW; //progress log file

	PSSMS 		= load_PSSM_DB(PSSM_DB, nprocs, rank);
	if (rank==0){
		P->display();
		FHW.open(log_out + job_ID + "-gTFI_log_file.txt");
	
	}

		
	
	//============================================================
	//....1.... load intervals from user, currently each interval is set to 2kb long
	t1=clock();
	map<string, vector<segment>> intervals 	= load_bed_file(bed_file, pad); 
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("loaded bed_file: %f seconds\n",(float(t2)-float(t1))/CLOCKS_PER_SEC  );
		FHW<<"loaded bed_file: "<<to_string(t)<<" seconds"<<endl;
		FHW.flush();	
	}
	//============================================================
	//....2.... insert the fasta sequnece into the provided intervals
	t1=clock();
	intervals 								= insert_fasta_sequence(fasta_file, intervals);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("loaded fasta file into bed file: %f seconds\n",(float(t2)-float(t1))/CLOCKS_PER_SEC );
		FHW<<"loaded fasta file into bed file: "<<to_string(t)<<" seconds"<<endl; 
		FHW.flush();	
	}
	if (intervals.empty()){
		printf("fasta file sequences was not inserted into data intervals....\n");
		return 0;
	}
	//============================================================
	//....3.... get generic GC content across intervals for background model
	t1=clock();
	vector<double> background  				= get_GC_content(intervals);

	vector<vector<double>> background_forward, background_reverse;

	get_ACGT_profile_all(intervals, 
		background_forward,background_reverse, rank);
		
	PSSMS 	=  construct_position_specific_pvalues(PSSMS, bins , background_forward  , background_reverse, site_br  );

	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("loaded background: %f seconds\n",t );
		FHW<<"loaded background: "<<to_string(t)<<" seconds"<<endl;
		FHW.flush();
	}
	t1=clock();

	// //============================================================
	// //....4.... do the dynamic programming for pvalue calculations
	// PSSMS 	= DP_pvalues(PSSMS,bins, background);
	// if (rank==0){
	// 	t2=clock();
	// 	float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
	// 	printf("computed p-values: %f seconds\n",(float(t2)-float(t1))/CLOCKS_PER_SEC );
	// 	FHW<<"computed p_values: "<<to_string(t)<<" seconds"<<endl;
	// 	FHW.flush();
	// }
	
	//============================================================
	//....5.... now scan for PSSM hits on provided intervals, <pv
	t1=clock();
	//intervals 	= run_accross(intervals, PSSMS,background, pv,rank);
	intervals 	= run_accross2(intervals, PSSMS,background_forward, 
		background_reverse, pv,rank, bed_out, job_ID);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("scanned intervals: %f seconds\n",t);
		FHW<<"scanned intervals: "<<to_string(t)<<" seconds"<<endl;
		FHW.flush();
	}

	// map<int, double> NN;
	// map<int, double [2000][4]> G;
	// cout.flush();
	// t1=clock();
	
	// //============================================================
	// //....7.... run simulations, independent multinomial draw, 
	// //collect stats on the fly, save on memory....

	// t1=clock();
	// map<int, vector<vector<double> >> observed_null_statistics;
	// map<int, map<int, vector<double> > > null_co_occur;
	
	// run_sims2(intervals, PSSMS, simN, rank, nprocs, background, pv, observed_null_statistics, null_co_occur);
	// if (rank==0){
	// 	t2=clock();
	// 	float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
	// 	printf("finished simulations: %f seconds\n", t );
	// 	FHW<<"finished simulations: "<<to_string(t)<<" seconds"<<endl;
	// 	FHW.flush();
	// }
	// t1=clock();
	
	// //========================================================================
	// //....8.... collect sample statistics on observations
	map<int, vector<double> > observed_statistics;
	map<int, vector<double>>  observed_displacements;
	map<int, map<int, double> > observed_co_occurrences; 

	collect_sample_stats(intervals, PSSMS,observed_statistics,
		observed_displacements,observed_co_occurrences, rank);

	map<int, vector<double> > GGG =  send_collect_observed_statistics( rank,  nprocs, PSSMS, observed_displacements);

	if (rank==0){
		map<int, string> 	 AA;
		load_PSSM_ID_names_only(PSSM_DB, AA);

		write_out_3(out_dir, job_ID,  AA, GGG);
	}
	// t1=clock();
	// if (rank==0){
	// 	collect_sample_stats(intervals, PSSMS,observed_statistics,
	// 		observed_displacements,observed_co_occurrences, rank);
	// 	if (rank==0){
	// 		t2=clock();
	// 		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
	// 		printf("finished collecting sample stats: %f seconds\n", t );
	// 		FHW<<"finished collecting sample stats: "<<to_string(t)<<" seconds"<<endl;
	// 		FHW.flush();
	// 	}
	// 	//now write out results
	// 	write_out_2( out_dir, job_ID,  PSSMS ,	  observed_statistics, 
	// 	  observed_displacements,  observed_co_occurrences,
	// 		 observed_null_statistics,   null_co_occur);
	// }

	// t1=clock();
	

	MPI::Finalize();


	return 0;
}