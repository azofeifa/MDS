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
using namespace std;

int main(int argc,char* argv[]){
	params * P = new params();
	fill_in_options(argv, P, 0);
	if (P->check()){
		printf("exiting...\n");
		return 0;
	}
	P->display();
	string fasta_file 			= P->p["-fasta"];
	string bed_file 			= P->p["-bed"];
	string out_dir 				= P->p["-o"];
	string log_out 				= P->p["-log_out"];
	string PSSM_DB 				= P->p["-DB"];
	string job_ID 				= P->p["-ID"];
	ofstream FHW;
	FHW.open(log_out + job_ID + "-gTFI_log_file.txt");
	int pad 					= 1500;
	int bins 					= stoi(P->p["-br"]);
	map<string, vector<segment>> intervals 	= load_bed_file(bed_file, pad);
	FHW<<"loaded bed_file"<<endl;
	FHW.flush();
	if (intervals.empty()){
		printf("user provided intervals were not loaded, exiting....\n");
	}
	vector<PSSM *> PSSMS 		= load_PSSM_DB(PSSM_DB);
	FHW<<"loaded PSSMS"<<endl;
	FHW.flush();
	if (PSSMS.empty()){
		printf("PSSMs were not loaded, exiting...\n");
		return 0;
	}
	intervals 								= insert_fasta_sequence(fasta_file, intervals);
	FHW<<"loaded intervals"<<endl;
	FHW.flush();
	if (intervals.empty()){
		printf("fasta file sequences was not inserted into data intervals....\n");
	}
	vector<double> background  				= get_GC_content(intervals);
	FHW<<"loaded background"<<endl;
	FHW.flush();


	DP_pvalues(PSSMS,bins, background);
	FHW<<"computed p_values"<<endl;
	FHW.flush();
	intervals 	= run_accross(intervals, PSSMS,background);
	FHW<<"ran across intervals"<<endl;
	FHW.flush();
	map<string, double> NN;
	map<string, double [3000][4]> G;
	get_average_ACGT_profile(intervals, PSSMS,pad, NN, G);
	FHW<<"gathered ACGT profiles for each motif"<<endl;
	FHW.flush();


	// //map<string, vector<vector<segment> >  > SIMS 	=  run_sims(G , NN, pad );


	// write_observation_stats(intervals, out_dir, job_ID, PSSMS, G, pad);
	// FHW<<"DONE :)"<<endl;
	// FHW.flush();


	return 0;
}