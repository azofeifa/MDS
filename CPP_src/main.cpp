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
	int simN 					= stoi(P->p["-simN"]);
	int pad 					= 1000;
	int bins 					= stoi(P->p["-br"]);
	double pv 					= stod(P->p["-pv"]);
	//============================================================

	ofstream FHW; //progress log file


	if (rank==0){ //load the PSSM table and send it out to everybody
		P->display();
		FHW.open(log_out + job_ID + "-gTFI_log_file.txt");
		
		PSSMS 		= load_PSSM_DB(PSSM_DB);
		if (PSSMS.empty()){
			printf("PSSMs were not loaded, exiting...\n");
			return 0;
		}
	
		FHW<<"loaded PSSMS"<<endl;
		FHW.flush();
		test_send_PSSMS(rank, nprocs,PSSMS,streamed_PSSMS,PSSM_IDS,PSSM_N);
	}else{
		test_send_PSSMS(rank, nprocs, PSSMS,streamed_PSSMS,PSSM_IDS,PSSM_N);
	}

	//============================================================
	//....1.... load intervals from user, currently each interval is set to 2kb long
	t1=clock();
	map<string, vector<segment>> intervals 	= load_bed_file(bed_file, pad); 
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("loaded bed_file: %f seconds\n",(float(t2)-float(t1))/CLOCKS_PER_SEC  );
		FHW<<"loaded bed_file: "<<to_string(t)<<"seconds"<<endl;
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
		FHW<<"loaded fasta file into bed file: "<<to_string(t)<<"seconds"<<endl; 
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
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("loaded background: %f seconds\n",t );
		FHW<<"loaded background: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}
	vector<PSSM *> 	new_PSSMS 		= convert_streatmed_to_vector(streamed_PSSMS, PSSM_IDS, PSSM_N);
	t1=clock();
	//============================================================
	//....4.... do the dynamic programming for pvalue calculations
	new_PSSMS 	= DP_pvalues(new_PSSMS,bins, background);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("computed p-values: %f seconds\n",(float(t2)-float(t1))/CLOCKS_PER_SEC );
		FHW<<"computed p_values: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}
	
	//============================================================
	//....5.... now scan for PSSM hits on provided intervals, <pv
	t1=clock();
	intervals 	= run_accross(intervals, new_PSSMS,background, pv,rank);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("scanned intervals: %f seconds\n",t);
		FHW<<"scanned intervals: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}

	map<int, double> NN;
	map<int, double [2000][4]> G;
	cout.flush();
	t1=clock();
	//============================================================
	//....6.... get the average ACGT profile according to each motif
	G 	= get_average_ACGT_profile(intervals, new_PSSMS,pad, NN, G, rank);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("gathering ACGT frequency profiles: %f seconds\n",t);
		FHW<<"gathered ACGT profiles for each motif: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}
	//============================================================
	//....7.... run simulations, independent multinomial draw

	t1=clock();
	map<int, map<int, vector<segment> >> S = run_sims(G , NN, new_PSSMS, simN,rank );
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("finished simulations: %f seconds\n", t );
		FHW<<"finished simulations: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}
	t1=clock();
	//============================================================
	//....8.... scan simulations for PSSM hits
	S 	= scan_simulations(S,new_PSSMS, background, pv );
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("finished scanning simulations: %f seconds\n", t );
		FHW<<"finished scanning simulations: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}
	//========================================================================
	//....9.... collect sample statistics on simulations and observations
	map<int, vector<double> > observed_statistics;
	map<int, vector<vector<double> >> observed_null_statistics;
	map<int, vector<double>>  observed_displacements;
	t1=clock();
	collect_sample_stats(intervals, S, new_PSSMS,observed_statistics,observed_null_statistics,observed_displacements,rank);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("finished collecting sample stats: %f seconds\n", t );
		FHW<<"finished collecting sample stats: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}

	t1=clock();
	map<int, vector< vector <double> >> collections 	= collect_PSSM_hits(rank, nprocs, intervals, observed_statistics, observed_null_statistics,observed_displacements);
	if (rank==0){
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("finished gathering results from MPI jobs: %f seconds\n", t );
		FHW<<"finished gathering results from MPI jobs: "<<to_string(t)<<"seconds"<<endl;
		FHW.flush();
	}

	if (rank==0){
		t1=clock();
		write_out(out_dir, collections, PSSMS);
		t2=clock();
		float t 	= (float(t2)-float(t1))/CLOCKS_PER_SEC ;
		printf("finished writing results: %f seconds\n", t );
		FHW<<"finished writing results: "<<to_string(t)<<"seconds"<<endl;
		printf("DONE :)\n");
		FHW<<"DONE :)"<<endl;
		FHW.flush();
	}

	MPI::Finalize();


	return 0;
}