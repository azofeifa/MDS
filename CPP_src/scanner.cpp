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


vector<int> get_sig_positions(vector<int> forward, 
			      vector<int> reverse, int N, PSSM * p, double pv){

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


		if (llf > p->ll_thresh){
			locs_pvs.push_back(i);	
		}
		else if (llr > p->ll_thresh){
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


double send_out_displacement_data(vector<int> & DD, vector<int> & DD_TSS,
	vector<int> & DD_NON , int rank, int nprocs, double TSS_spec_association, vector<int> & special_hits){
	double final_association 	= TSS_spec_association;
	if (rank==0){
		for (int j = 1 ; j < nprocs; j++){
			//===========================================
			//recieve All

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
			int S3;
			//===========================================
			//recieve spec interval hits (for bed_out)
			MPI_Recv(&S3, 1, MPI_INT, j, 4, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if (S3>0){
				int * B = new int[S3];
				MPI_Recv(&B[0], S3, MPI_INT, j, 5, MPI_COMM_WORLD,MPI_STATUS_IGNORE);					
				vector<int> temp 	= to_vector_2(B, S3);
				special_hits.insert(special_hits.end(), temp.begin(), temp.end());
					
			}
			
			int S_TSS, S_NON;
			
			//===================================
			//revieve TSS displacements
			MPI_Recv(&S_TSS, 1, MPI_INT, j, 6, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if (S_TSS>0){
				int * A = new int[S_TSS];
				MPI_Recv(&A[0], S_TSS, MPI_INT, j, 7, MPI_COMM_WORLD,MPI_STATUS_IGNORE);					
				vector<int> temp 	= to_vector_2(A, S_TSS);
				DD_TSS.insert(DD_TSS.end(), temp.begin(), temp.end());
			}

			//===================================
			//revieve nonTSS displacements
			MPI_Recv(&S_NON, 1, MPI_INT, j, 8, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if (S_NON>0){
				int * A = new int[S_NON];
				MPI_Recv(&A[0], S_NON, MPI_INT, j, 9, MPI_COMM_WORLD,MPI_STATUS_IGNORE);					
				vector<int> temp 	= to_vector_2(A, S_NON);
				DD_NON.insert(DD_NON.end(), temp.begin(), temp.end());
			}
			




		}
	}else{
		int S 		= DD.size();
		int S_TSS 	= DD_TSS.size();
		int S_NON 	= DD_NON.size();
		
		//===========================================
		//Send All
		MPI_Ssend(&S, 1, MPI_INT, 0,1, MPI_COMM_WORLD);
		MPI_Ssend(&TSS_spec_association, 1, MPI_DOUBLE, 0,2, MPI_COMM_WORLD);
		
		if (S>0){
			int * A = new int[S];
			fill_array(DD, A);
			MPI_Ssend(&A[0], S, MPI_INT, 0,3, MPI_COMM_WORLD);
		}
		//===========================================
		//send spec interval hits (for bed_out)

		int S3 	=	 special_hits.size();
		MPI_Ssend(&S3, 1, MPI_INT, 0,4, MPI_COMM_WORLD);
		if (S3>0){
			int * B = new int[S3];
			fill_array(special_hits, B);
			MPI_Ssend(&B[0], S3, MPI_INT, 0,5, MPI_COMM_WORLD);
		}
		//===================================
		//send TSS displacements
		MPI_Ssend(&S_TSS, 1, MPI_INT, 0,6, MPI_COMM_WORLD);
		
		if (S_TSS>0){
			int * A = new int[S_TSS];
			fill_array(DD_TSS, A);
			MPI_Ssend(&A[0], S_TSS, MPI_INT, 0,7, MPI_COMM_WORLD);
		}

		

		//===================================
		//send nonTSS displacements
		MPI_Ssend(&S_NON, 1, MPI_INT, 0,8, MPI_COMM_WORLD);
		
		if (S_NON>0){
			int * A = new int[S_NON];
			fill_array(DD_NON, A);
			MPI_Ssend(&A[0],S_NON, MPI_INT, 0,9, MPI_COMM_WORLD);
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
		    int rank, int nprocs, Log_File * LG, int MD_window,double TSS_association, string out2){

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
	vector<vector<int>> array_of_final_displacements_TSS(PSSMS.size());
	vector<vector<int>> array_of_final_displacements_NON(PSSMS.size());



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
		vector<int> final_displacements_NON;
		vector<int> final_displacements_TSS;
		
		vector<int> special_hits;
		double TSS_spec_association 	= 0;
		double N 	= 0.01;
		for (int i =0 ; i < displacements.size(); i++){
			double MIN 	= 2000;
			for (int j = 0 ; j < displacements[i].size(); j++ ){
				if (D[start+i].TSS){
					TSS_spec_association++;
				}
				N++;			
				if (abs(displacements[i][j] - 1000) < MIN){
					MIN 	= abs(displacements[i][j] - 1000);
				}
				final_displacements.push_back(displacements[i][j]);
				if (D[start+i].TSS){
					final_displacements_TSS.push_back(displacements[i][j]);
				}else{
					final_displacements_NON.push_back(displacements[i][j]);	
				}
			}
			if (MIN < MD_window){
				special_hits.push_back(start + i);
			}
		}
		TSS_spec_association/=N;
		TSS_spec_association=send_out_displacement_data(final_displacements,
			final_displacements_TSS,final_displacements_NON, rank, nprocs,TSS_spec_association, special_hits );

		PSSMS[p]->TSS_association 	= TSS_spec_association;
		if (rank==0){
			array_of_final_displacements[p] 		= final_displacements;
			array_of_final_displacements_TSS[p] 	= final_displacements_TSS;
			array_of_final_displacements_NON[p] 	= final_displacements_NON;

			for (int i = 0 ; i < special_hits.size(); i++){
				D[special_hits[i]].motif_hits[PSSMS[p]->name] 	= 1;
			}
		}
		t = clock() - t;
		LG->write("done: " + to_string(float(t)/(CLOCKS_PER_SEC)) + " seconds (" + to_string(p+1) + "/" + to_string(PSSMS.size())+"), " + to_string(TSS_spec_association)+"\n", 1);
	}
	if (rank==0){
		if (not out2.empty()){
			LG->write("writing out bed hits........................", 1);
			write_out_bed_file(D, out2 , MD_window);
			LG->write("done\n",1);
		}
		LG->write("computing boostrapped distribution..........", 1);

		t = clock();
		#pragma omp parallel for
		for (int p = 0 ; p < PSSMS.size(); p++){
			vector<int> final_displacements 	= array_of_final_displacements[p];
			vector<int> final_displacements_TSS = array_of_final_displacements_TSS[p];
			vector<int> final_displacements_NON = array_of_final_displacements_NON[p];

			double MD_score 			= get_MD_score(final_displacements,MD_window,true);
			double MD_score_TSS 		= get_MD_score(final_displacements_TSS,MD_window,true);
			double MD_score_NON 		= get_MD_score(final_displacements_NON,MD_window,true);

			vector<double> md_scores_many 		= get_many_MD_scores(final_displacements, 10 );
			vector<double> md_scores_many_TSS 	= get_many_MD_scores(final_displacements_TSS, 10);
			vector<double> md_scores_many_NON 	= get_many_MD_scores(final_displacements_NON, 10);


			double NN 				= final_displacements.size();
			double NN_TSS 			= final_displacements_TSS.size();
			double NN_NON 			= final_displacements_NON.size();
			
			PSSMS[p]->MD_score 		= MD_score;
			PSSMS[p]->MD_score_TSS 	= MD_score_TSS;
			PSSMS[p]->MD_score_NON 	= MD_score_NON;


			PSSMS[p]->md_many 		= md_scores_many;
			PSSMS[p]->md_many_TSS 	= md_scores_many_TSS;
			PSSMS[p]->md_many_non 	= md_scores_many_NON;



			PSSMS[p]->total 		= NN;
			PSSMS[p]->total_TSS 	= NN_TSS;
			PSSMS[p]->total_NON 	= NN_NON;

			build_cdfs_PSSMs(PSSMS[p], bsn, interval_size, NN, NN_TSS, NN_NON, MD_window,0.01);
			


			PSSMS[p]->get_pvalue_stats(2*MD_window/2000.0);
			PSSMS[p]->observed_displacements 		= final_displacements;
			PSSMS[p]->observed_displacements_TSS 	= final_displacements_TSS;
			PSSMS[p]->observed_displacements_non 	= final_displacements_NON;
			
		}
		t = clock() - t;
		LG->write("done: " + to_string(float(t)/(CLOCKS_PER_SEC*threads)) + " seconds\n", 1);
	}
}

void scan_accross_3(map<string, vector<segment> > S , PSSM * p,string out_dir, double pv,vector<double> background ){
	map<char, int> table;
	map<int, int>flip;
	map<char, char> UC;
	map<char, char> SS;



	flip[0] 	=3, flip[3] = 0, flip[1]=2, flip[2]=1;
		

	UC['A'] 	= 'A';
	UC['T'] 	= 'T';
	UC['G'] 	= 'G';
	UC['C'] 	= 'C';

	UC['a'] 	= 'A';
	UC['t'] 	= 'T';
	UC['g'] 	= 'G';
	UC['c'] 	= 'C';

	SS['A'] 	= 'T';
	SS['T'] 	= 'A';
	SS['G'] 	= 'C';
	SS['C'] 	= 'G';



	table['A'] 	= 0, table['C']=1, table['G']=2, table['T']=3;
	table['a'] 	= 0, table['c']=1, table['g']=2, table['t']=3;





	ofstream FHW(out_dir+p->name + ".bed");
	double threshold  	= p->get_threshold(pv);
	FHW<<"#p-value\t"<<to_string(log10(pv))<<endl;
	FHW<<"#llr\t"<<to_string(threshold)<<endl;
	
	FHW<<"#motif\t" <<p->name<<endl;
	for (int i = 0 ; i < p->frequency_table.size();i++){
		string 	linef 	= "#";
		double x;
		for (int s =0; s < 3;s++){
			x = (p->frequency_table[i][s]/2.) + log(background[s]);
			x 		= exp(x);
			linef+=to_string(x) + ",";
		}
		x = (p->frequency_table[i][3]/2.) + log(background[3]);
		x 		= exp(x);
		linef+=to_string(x) + "\n";
		FHW<<linef;

	}
	vector<int> locs_pvs;
	int length 	= p->frequency_table.size();
	double pvaluef, pvaluer,llr,llf;
	typedef map<string, vector<segment> >::iterator it_type_1;
	for (it_type_1 c=S.begin(); c!=S.end();c++){
		for (int s = 0 ; s < c->second.size(); s++){
			string seq 	= c->second[s].seq;
			int N 		= seq.size();
			int BOUND 	= N-length;
			int k,j,i,f, r, l;
			string current 	= "",line="",rcurrent="";
			for (i =0 ; i  < BOUND; i++){
				llf 	= 0;
				llr 	= 0;
				k 		= 0;
				j 		= i;
				l 		= i + length-1;
				current = "";
				rcurrent= "";
				for (k=0; k < length; k++){
					if (table.find(seq[j])!=table.end() ){
						f 	= table[seq[j]], r 	= flip[table[seq[l]]];
						llf+= (p->frequency_table[k][f]) ;
						llr+= (p->frequency_table[k][r]) ;
						j++;
						l--;
					}else{
						llf=-10000,llr=-10000;
						break;
					}
				}
				// pvaluef 	= 1.0-p->get_pvalue(llf);
				// pvaluer 	= 1.0-p->get_pvalue(llr);
				if (threshold < llf){
					j 	= i;
					for (k=0; k < length; k++){
						if (table.find(seq[j])!=table.end() ){
							current+=UC[seq[j]];
						}
						j++;
					}					


					line 	 = c->second[s].chrom+"\t" +to_string(i+c->second[s].start) + "\t";
					line 	+= to_string(i+length+c->second[s].start) + "\t";
					line	+= current+"\t" + "+" +"\t" + to_string(llf)+ "\n";
					FHW<<line;
				}
				if (threshold < llr){
					l 		= i + length-1;
					for (k=0; k < length; k++){
						if (table.find(seq[j])!=table.end() ){
							rcurrent+=SS[UC[seq[l]]];
						}
						l--;
					}					

					line 	 = c->second[s].chrom+"\t" +to_string(i+c->second[s].start) + "\t";
					line 	+= to_string(i+length+c->second[s].start) + "\t";
					line	+= rcurrent+"\t" + "-" + "\t" + to_string(llr)+"\n";
					FHW<<line;
				}


			}					


		}
	}

	FHW.close();
}




void scan_intervals_genome_wide(map<string, vector<segment>> S, vector<PSSM *> PSSMS, vector<double> background, 
									double pv, 
									int rank, int nprocs, Log_File * LG, string out_dir ){
	typedef map<string, vector<segment>>::iterator it_type;
	#pragma omp parallel for	
	for (int p=0; p < PSSMS.size(); p++){
		scan_accross_3(S, PSSMS[p], out_dir, pv,background );
	}
}



























