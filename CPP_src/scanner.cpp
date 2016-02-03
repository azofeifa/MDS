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
	for (int i = 0 ; i < ATGC.size(); i++){
		ATGC[i]/=SUM;
	}
	return ATGC;
}


vector<double> get_sig_positions(int forward[2000], 
	int reverse[2000], int N, PSSM * p, vector<double> background, double pv){

	vector<double> locs_pvs;
	int length 	= p->frequency_table.size();
	double pvaluef, pvaluer,llr,llf;
	int k,j;
	bool collect;
	int BOUND 	= N-length;
	int i;
	int f, r, l;
	for (i =0 ; i  < BOUND; i++){
		llf 	= 0;
		llr 	= 0;
		k 		= 0;
		j 		= i;
		l 		= i + length-1;
		for (k=0; k < length; k++){
			f 	= forward[j], r 	= reverse[l];
			llf+= (p->frequency_table[k][f]) - (background[f]);
			llr+= (p->frequency_table[k][r]) - (background[r]);
			j++;
			l--;
		}
		pvaluef 	= 1.0-p->get_pvalue(llf*2);
		pvaluer 	= 1.0-p->get_pvalue(llr*2);
		if (pvaluef < pv){
			locs_pvs.push_back(i);
		}
		if (pvaluer < pv){
			locs_pvs.push_back(i);
		}
	}
	return locs_pvs;
}

map<int, vector<double>> wrapper(segment & S, vector<PSSM *> PSSMS, 
	vector<double> background, double pv){
	vector<vector<double>> sig_positions(PSSMS.size());
	//#pragma omp parallel for
	for (int p = 0 ; p < PSSMS.size(); p++){
		sig_positions[p]= get_sig_positions(S.forward, S.reverse, 2000, 
			PSSMS[p], background, pv);
	}
	map<int, vector<double>> sig_positions_map;
	for (int i = 0 ; i < PSSMS.size();i++){
		sig_positions_map[PSSMS[i]->ID]=sig_positions[i];
	}
	return sig_positions_map;
}








map<string, vector<segment> > run_accross(map<string, vector<segment>> S ,
 vector<PSSM *> PSSMS, vector<double> background, double pv, int rank){
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
			c->second[i].motif_positions = wrapper(c->second[i], PSSMS,background, pv);
			newS[c->first].push_back(c->second[i]);
		}
	}
	return newS;
}



int PSSM_index(int i, vector<PSSM *> P){
	for (int j = 0 ; j < P.size(); j++){
		if (P[j]->ID==i){
			return j;
		}
	}
	printf("what???\n");
	return 0;
}

map<int, map<int, vector<segment> >> scan_simulations(map<int, map<int, vector<segment> >> S,
	vector<PSSM *> P, vector<double> background,double pv ){
	for (int b = 0 ; b < background.size(); b++){
		background[b] 	= log(background[b]);
	}
	typedef map<int, map<int, vector<segment> >>::iterator it_type; 
	typedef map<int, vector<segment> >::iterator it_type_2; 
	
	map<int, map<int, vector<segment> >> newS;
	for (it_type s = S.begin(); s!=S.end(); s++ ){
		int p_index 		= PSSM_index(s->first, P);

		for (it_type_2 cn 	= s->second.begin(); cn!=s->second.end(); cn++){
			vector<segment> currents(cn->second.size());
			#pragma omp parallel for
			for (int i = 0 ; i < cn->second.size(); i++){
				cn->second[i].motif_positions[s->first]=get_sig_positions(cn->second[i].forward,
					cn->second[i].reverse, 2000,P[p_index], background, pv );
				currents[i] 	= cn->second[i];
			}

			newS[s->first][cn->first]= currents;

		}
	}
	return newS;
	

}
























