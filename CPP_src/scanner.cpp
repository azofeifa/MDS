#include "scanner.h"
#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;

vector<double> get_GC_content(map<string, vector<segment>> S){
	vector<double> ATGC = {1,1,1,1};
	map<char, int> table;
	table['A'] 	= 0, table['C']=1, table['G']=2, table['T']=3;
	table['a'] 	= 0, table['c']=1, table['g']=2, table['t']=3;
	
	string white_space 	= "";
	for (int i = 0 ; i < S.size(); i++){
		white_space+= " ";
	}
	int counter=0;
	string stars 	= "";

	typedef map<string, vector<segment>>::iterator it_type; 
	for (it_type c = S.begin(); c!=S.end(); c++){
		
		printf("\rcomputing ACGT Frequency|%s%s|", stars.c_str(), white_space.substr(0,white_space.size()-counter).c_str());
		cout.flush();
		for (int i = 0 ; i < c->second.size(); i++ ){
			string seq 	= c->second[i].seq;
			for (int j = 0 ; j < seq.size(); j++){
				if (table.find(seq[j]) != table.end()){
					ATGC[table[seq[j]]]+=1;
				}
			}
		}
		stars+="*";
		counter+=1;
	}
	double SUM 	= 0;
	for (int i = 0 ; i < ATGC.size(); i++){
		SUM+=ATGC[i];
	}
	printf("\rcomputing ACGT Frequency|%s%s|done, ", stars.c_str(), white_space.substr(0,white_space.size()-counter).c_str());
	for (int i = 0 ; i < ATGC.size(); i++){
		ATGC[i]/=SUM;
		printf("%f, ",ATGC[i] );
	}
	printf("\n");
	return ATGC;
}


vector<double> get_sig_positions(int * forward, 
	int * reverse, int N, PSSM * p, vector<double> background){

	vector<double> locs_pvs;
	int length 	= p->frequency_table.size();
	double pvaluef, pvaluer,llr,llf;
	int k,j;
	bool collect;
	int BOUND 	= N-length;
	int i 		= 0;
	while (i < BOUND){
		llf 	= 0;
		llr 	= 0;
		k 		= 0;
		collect = true;
		j 		= i;
		while (k < length){
			llf+= (p->frequency_table[k][forward[j]]) - (background[forward[j]]);
			llr+= (p->frequency_table[k][reverse[j]]) - (background[reverse[j]]);
			j++;
			k++;
		}
		if (collect){
			pvaluef 	= 1.0-p->get_pvalue(llf*2);
			pvaluer 	= 1.0-p->get_pvalue(llr*2);
			if (pvaluef < pow(10,-6)){
				locs_pvs.push_back(i);
			}
			if (pvaluer < pow(10,-6)){
				locs_pvs.push_back(i);
			}
		}
		i++;
	}
	return locs_pvs;
}

vector<vector<double>> wrapper(segment & S, vector<PSSM *> PSSMS, vector<double> background){
	vector<vector<double>> sig_positions(PSSMS.size());
	#pragma omp parallel for
	for (int p = 0 ; p < PSSMS.size(); p++){
		sig_positions[p]= get_sig_positions(S.forward, S.reverse, S.N, PSSMS[p], background);
	}
	return sig_positions;
}




map<string, vector<segment> > run_accross(map<string, vector<segment>> S , vector<PSSM *> PSSMS, vector<double> background){
	typedef map<string, vector<segment>>::iterator it_type; 
	for (int p = 0 ; p < PSSMS.size(); p++){
		for (int i = 0 ; i<PSSMS[p]->frequency_table.size(); i++){
			for (int j = 0 ; j<PSSMS[p]->frequency_table[i].size(); j++){
				PSSMS[p]->frequency_table[i][j] 	= log(PSSMS[p]->frequency_table[i][j]);
			}
		}
	}	
	for (int b = 0 ; b < background.size(); b++){
		background[b] 	= log(background[b]);
	}
	map<string, vector<segment> > newS;
	for (it_type c = S.begin(); c!=S.end(); c++){
		for (int i = 0 ; i < c->second.size(); i++){
			c->second[i].motif_positions = wrapper(c->second[i], PSSMS,background);
			newS[c->first].push_back(c->second[i]);
		}
	}
	return newS;
}
