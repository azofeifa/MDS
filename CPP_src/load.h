#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>
#include "read_in_parameters.h"

using namespace std;
class PSSM{
public:
	string name;
	int N;
	int ID;
	int SN;
	vector<vector<double>> frequency_table;
	vector<vector<double>> pvalues;
	vector<vector<int>> null_displacements;
	vector<vector<int>> null_displacements_non;
	
	vector<int> null_displacements_2;
	vector<int> null_displacements_2_non;

	
	vector<int> observed_displacements; 
	vector<int> observed_displacements_TSS; 
	vector<int> observed_displacements_non; 


	vector<int> binned_observed_displacements; 
	vector<int> binned_observed_displacements_TSS; 
	vector<int> binned_observed_displacements_non; 


	
	vector<int> binned_null_displacements; 
	vector<int> binned_null_displacements_non; 


	double zeros; 
	double zeros_non;

	double TSS_association;
	double stationary_pvalue;

	map<int, vector<vector<double>>> position_specific_pvalues_forward;
	map<int, vector<vector<double>>> position_specific_pvalues_reverse;
	
	vector<vector<double> > MD_CDF;
	vector<vector<double> > MD_CDF_NON;
	vector<vector<double> > MD_CDF_TSS;


	vector<vector<double> > ENRICH_CDF;

	vector<double> non_pvalues, stationary_p_values;
	double ENRICH_score;
	
	double MD_score;
	double MD_score_TSS;
	double MD_score_NON;

	double pv_enrich_score_rt, pv_MD_score_rt, pv_enrich_score_lt, pv_MD_score_lt;

	double total;
	double total_TSS;
	double total_NON;


	PSSM();
	PSSM(string);
	PSSM(int);
	void bin_observations();
	void bin_null_displacements();
	void get_pvalue_stats(double);
	void transform_ft();
	double get_pvalue(double);
	double get_pvalue2_f(double,int, int);
	double get_pvalue2_r(double,int, int);
	double get_position_specific_pvalue(int , double);
	string get_consensus();

};
class segment{

public:
	string chrom; 
	int start, stop; 
	int rstart, rstop;
	int forward[2000];
	int reverse[2000];
	int position;
	string seq; 
	bool TSS;
	int N;
	map<int, vector<double>> motif_positions;
	map<string, int> motif_hits;
	segment();
	segment(string, int, int,int, int, int);
	bool transform();


};

vector<PSSM *> load_PSSM_DB(string, int, int, int);

map<string, vector<segment>> load_bed_file(string, int, int &, double &) ;

map<string, vector<segment> > insert_fasta_sequence(string , map<string, vector<segment> > , int);

vector<PSSM *> convert_streatmed_to_vector(vector<vector<vector<double>>>,
	vector<int>, vector<int>);
void load_PSSM_ID_names_only(string , map<int, string> &  );
vector<PSSM *> load_PSSM_DB_new(string, int);

void write_out_null_stats(vector<PSSM *> , string, params * , vector<double> ,vector<vector<double>>,
					vector<vector<double>>,vector<double ** >,vector<double ** > );




vector<PSSM *> load_personal_DB_file(string , params *, vector<double> & );
void write_out_stats(vector<PSSM *>, string, params *, double, double, double);

void collect_all_tmp_files(string , string , int  , int );

map<string, vector<segment>> label_TSS(map<string, vector<segment>> , map<string, vector<segment>>,double & );


void write_out_bed_file(vector<segment>, string,int );


#endif
