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
	int ID;
	int SN;
	vector<vector<double>> frequency_table;
	vector<vector<double>> pvalues;
	PSSM();
	PSSM(string);
	PSSM(int);
	void transform_ft();
	double get_pvalue(double);
	string get_consensus();

};
class segment{

public:
	string chrom; 
	int start, stop; 
	int forward[2000];
	int reverse[2000];
	int position;
	string seq; 
	int N;
	map<int, vector<double>> motif_positions;
	segment();
	segment(string, int, int,int);
	bool transform();


};

vector<PSSM *> load_PSSM_DB(string);

map<string, vector<segment>> load_bed_file(string, int) ;

map<string, vector<segment> > insert_fasta_sequence(string , map<string, vector<segment> > );

vector<PSSM *> convert_streatmed_to_vector(vector<vector<vector<double>>>,
	vector<int>, vector<int>);

#endif
