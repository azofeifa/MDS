#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>
using namespace std;
class PSSM{
public:
	string name;
	int N;
	int SN;
	vector<vector<double>> frequency_table;
	double pvalues[200][2];
	PSSM();
	PSSM(string);
	void transform_ft();
	double get_pvalue(double);
	string get_consensus();

};
class segment{

public:
	string chrom; 
	int start, stop; 
	int forward[3000];
	int reverse[3000];
	string seq; 
	int N;
	vector<vector< double>> motif_positions;
	segment();
	segment(string, int, int);
	bool transform();


};

vector<PSSM *> load_PSSM_DB(string);

map<string, vector<segment>> load_bed_file(string, int) ;

map<string, vector<segment> > insert_fasta_sequence(string , map<string, vector<segment> > );
#endif
