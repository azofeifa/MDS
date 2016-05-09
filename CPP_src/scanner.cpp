#include "scanner.h"
#include <iostream>
#include <fstream>
#include "collect_sample_statistics.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>

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


vector<int> get_sig_positions(int forward[2000], 
	int reverse[2000], int N, PSSM * p, double pv){

	vector<int> locs_pvs;
	int length 	= p->frequency_table.size();
	double pvaluef, pvaluer,llr,llf;

	int BOUND 	= N-length;
	int k,j,i,f, r, l;
	for (i =0 ; i  < BOUND; i++){
		llf 	= 0;
		llr 	= 0;
		k 		= 0;
		j 		= i;
		l 		= i + length-1;
		for (k=0; k < length; k++){
			f 	= forward[j], r 	= reverse[l];
			llf+= (p->frequency_table[k][f]) ;
			llr+= (p->frequency_table[k][r]) ;
			j++;
			l--;
		}


		pvaluef 	= 1.0-p->get_pvalue(llf);
		pvaluer 	= 1.0-p->get_pvalue(llr);
		if (pvaluef < pv){
			locs_pvs.push_back(i);	
		}
		if (pvaluer < pv){
			locs_pvs.push_back(i);
		}
	}
	return locs_pvs;
}



vector<segment> collapse_map(map<string, vector<segment>> S){
	typedef map<string, vector<segment>>::iterator it_type; 
	vector<segment> DD ;	
	for (it_type c = S.begin(); c!=S.end(); c++){
		for (int i = 0 ; i < c->second.size(); i++){
			DD.push_back(c->second[i]);
		}
	}
	return DD;
}
vector<int> to_vector_2(int * array, int S){
	vector<int> A(S);
	for (int s = 0 ; s < S; s++){
		A[s] 	= array[s];
	}
	return A;
}
void fill_array(vector<int> DD, int * A){
	for (int i = 0; i < DD.size(); i++){
		A[i] 	= DD[i];
	}
}


double send_out_displacement_data(vector<int> & DD, int rank, int nprocs, double TSS_spec_association){
	double final_association 	= TSS_spec_association;
	if (rank==0){
		for (int j = 1 ; j < nprocs; j++){
			int S;
			double S2;
			MPI_Recv(&S, 1, MPI_INT, j, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&S2, 1, MPI_DOUBLE, j, 2, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			final_association+=S2;
			if (S>0){
				int * A = new int[S];
				MPI_Recv(&A[0], S, MPI_INT, j, 3, MPI_COMM_WORLD,MPI_STATUS_IGNORE);					
				vector<int> temp 	= to_vector_2(A, S);
				DD.insert(DD.end(), temp.begin(), temp.end());
			}
		}
	}else{
		int S 	= DD.size();

		MPI_Ssend(&S, 1, MPI_INT, 0,1, MPI_COMM_WORLD);
		MPI_Ssend(&TSS_spec_association, 1, MPI_DOUBLE, 0,2, MPI_COMM_WORLD);
		
		if (S>0){
			int * A = new int[S];
			fill_array(DD, A);
			MPI_Ssend(&A[0], S, MPI_INT, 0,3, MPI_COMM_WORLD);
		}
	}
	return final_association/nprocs;
}

string get_dots_2(int N){
	string line="";
	for (int i = 0 ; i < N; i++){
		line+=".";
	}
	return line;
}


void scan_intervals(map<string, vector<segment>> S ,
 												vector<PSSM *> PSSMS, vector<double> background, 
 												double pv, int interval_size, int bsn, 
 												int rank, int nprocs, Log_File * LG, int MD_window,double TSS_association){

	vector<segment> D 	= collapse_map(S);
	int threads  		= omp_get_max_threads();
	int diff 			= D.size() / nprocs;
	int start 	= rank*diff, stop 	= min(int( (rank+1)*diff),int(D.size()));
	if (rank+1==nprocs){
		stop 	= D.size();
	}
	clock_t t;
	LG->write("(this MPI call will process " + to_string(stop-start) + " intervals)\n\n", 1);
	vector<vector<int>> array_of_final_displacements(PSSMS.size());
	for (int p = 0 ; p < PSSMS.size(); p++){
		vector<vector<int>> displacements(int(stop-start));
		int WN 	= max(int(44 - PSSMS[p]->name.size()), 1);	
		t = clock();
		LG->write(PSSMS[p]->name + get_dots_2(WN), 1);
		#pragma omp parallel for
		for (int i = 0 ; i < stop-start; i++ ){
			displacements[i] 	= get_sig_positions(D[start+i].forward, D[start+i].reverse, 2000, PSSMS[p], pv);
		}
		vector<int> final_displacements;
		double TSS_spec_association 	= 0;
		for (int i =0 ; i < displacements.size(); i++){
			if(D[i].TSS){
				TSS_spec_association+=int(displacements[i].empty());
			}
			for (int j = 0 ; j < displacements[i].size(); j++ ){
				final_displacements.push_back(displacements[i][j]);
			}
		}
		TSS_spec_association/=displacements.size();
		TSS_spec_association=send_out_displacement_data(final_displacements, rank, nprocs,TSS_spec_association );

		PSSMS[p]->TSS_association 	= TSS_spec_association;
		if (rank==0){
			array_of_final_displacements[p] 	= final_displacements;
		}
		t = clock() - t;
		LG->write("done: " + to_string(float(t)/(CLOCKS_PER_SEC)) + " seconds (" + to_string(p+1) + "/" + to_string(PSSMS.size())+"), " + to_string(TSS_spec_association)+"\n", 1);
	}
	if (rank==0){
		LG->write("computing boostrapped distribution..........", 1);

		t = clock();
		#pragma omp parallel for
		for (int p = 0 ; p < PSSMS.size(); p++){
			vector<int> final_displacements 	= array_of_final_displacements[p];
			double MD_score 		= get_MD_score(final_displacements,MD_window,true);
			double ENRICH_score 	= get_MD_score(final_displacements,MD_window,false);
			double NN 				= final_displacements.size();
			PSSMS[p]->MD_score 		= MD_score;
			PSSMS[p]->ENRICH_score 	= ENRICH_score;		
			PSSMS[p]->total 		= NN;
			
			build_cdfs_PSSMs(PSSMS[p], bsn, interval_size, NN, MD_window,PSSMS[p]->TSS_association);
			PSSMS[p]->get_pvalue_stats();
			PSSMS[p]->observed_displacements 	= final_displacements;
		}
		t = clock() - t;
		LG->write("done: " + to_string(float(t)/(CLOCKS_PER_SEC*threads)) + " seconds\n", 1);
	}



}



























