#include "read_in_parameters.h"
#include "load.h"
#include "get_motif_pvalues.h"
#include "scanner.h"
#include "out.h"
#include <vector>
#include <string>
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
	int pad 					= stoi(P->p["-pad"]);
	int bins 					= stoi(P->p["-br"]);
	map<string, vector<segment>> intervals 	= load_bed_file(bed_file, pad);
	if (intervals.empty()){
		printf("user provided intervals were not loaded, exiting....\n");
	}
	vector<PSSM *> PSSMS 		= load_PSSM_DB(PSSM_DB);
	if (PSSMS.empty()){
		printf("PSSMs were not loaded, exiting...\n");
		return 0;
	}
	intervals 								= insert_fasta_sequence(fasta_file, intervals);
	if (intervals.empty()){
		printf("fasta file sequences was not inserted into data intervals....\n");
	}
	vector<double> background  				= get_GC_content(intervals);


	DP_pvalues(PSSMS,bins, background);
	intervals 	= run_accross(intervals, PSSMS,background);

	write_observation_stats(intervals, out_dir, job_ID, PSSMS);


	return 0;
}