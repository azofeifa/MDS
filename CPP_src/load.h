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
	vector<int> null_displacements_2;
	vector<int> observed_displacements; 
	vector<int> binned_observed_displacements; 
	vector<int> binned_null_displacements; 

	double zeros; 
	map<int, vector<vector<double>>> position_specific_pvalues_forward;
	map<int, vector<vector<double>>> position_specific_pvalues_reverse;
	vector<vector<double> > MD_CDF;
	vector<vector<double> > ENRICH_CDF;
	double ENRICH_score;
	double MD_score;
	double pv_enrich_score_rt, pv_MD_score_rt, pv_enrich_score_lt, pv_MD_score_lt,total;
	PSSM();
	PSSM(string);
	PSSM(int);
	void bin_observations();
	void bin_null_displacements();
	void get_pvalue_stats();
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
	int N;
	map<int, vector<double>> motif_positions;
	segment();
	segment(string, int, int,int, int, int);
	bool transform();


};

vector<PSSM *> load_PSSM_DB(string, int, int, int);

map<string, vector<segment>> load_bed_file(string, int, int &) ;

map<string, vector<segment> > insert_fasta_sequence(string , map<string, vector<segment> > , int);

vector<PSSM *> convert_streatmed_to_vector(vector<vector<vector<double>>>,
	vector<int>, vector<int>);
void load_PSSM_ID_names_only(string , map<int, string> &  );
vector<PSSM *> load_PSSM_DB_new(string, int);
void write_out_null_stats(vector<PSSM *> , string, params * , vector<double> );
vector<PSSM *> load_personal_DB_file(string , params *, vector<double> & );
void write_out_stats(vector<PSSM *>, string, params *);

void collect_all_tmp_files(string , string , int  , int );



#endif
