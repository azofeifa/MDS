#include "load.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>

#include "split.h"
using namespace std;

segment::segment(){};
segment::segment(string chr, int st, int sp, int ID, int rst, int rsp){
	chrom=chr, start=st, stop=sp;
	seq 	= "";
	position=ID;
	rstart = rst, rstop = rsp;
	TSS 	= false;
}




vector<segment> sort(vector<segment> segments){
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < segments.size(); i++){
			if (segments[i-1].stop > segments[i].stop){
				changed 		= true;
				segment copy 	= segments[i-1];
				segments[i-1] 	= segments[i];
				segments[i] 	= copy;
			}
		}
	}
	return segments;


}

bool segment::transform(){
	map<char, int> table;
	map<int, int>flip;
	table['A'] 	= 0, table['C']=1, table['G']=2, table['T']=3;
	table['a'] 	= 0, table['c']=1, table['g']=2, table['t']=3;
	flip[0] 	=3, flip[3] = 0, flip[1]=2, flip[2]=1;
	N 			= seq.size();
	bool keep 	= true;
	for (int j = 0 ; j < seq.size(); j++){

		if (table.find(seq[j]) != table.end()){
			forward[j] 	= table[seq[j]];
			reverse[j] 	= flip[table[seq[j]]];
		}else{
			forward[j] 	= 5;
			reverse[j] 	= 5;
			keep 		= false;
		}
	}
	return keep;
}



map<string, vector<segment>> load_bed_file(string FILE, int pad, int & N){
	map<string, vector<segment>> S ;
	ifstream FH(FILE);

	if (FH){
		string line, chrom;
		vector<string> line_array;
		int start, stop;
		double x;
		int i 	= 0;
		while (getline(FH, line)){
			if (line.substr(0,1)!="#"){
				line_array=splitter(line, "\t");
				chrom 	= line_array[0], start = stoi(line_array[1]), stop = stoi(line_array[2]);
				N++;
				x 	= (start + stop)/2.;
				S[chrom].push_back(segment(chrom, x-pad, x+pad,i, start, stop));
				i++;
			}
		}
		typedef map<string, vector<segment> >::iterator it_type;
		for (it_type i = S.begin(); i  !=S.end();i++ ){
			S[i->first] 	= sort(i->second);
		}

	}else{
		printf("couldn't open %s \n", FILE.c_str() );
	}

	return S;	
}
map<string, vector<segment> > insert_fasta_sequence(string fasta_file, map<string, vector<segment> > S, int test){
	ifstream FH(fasta_file);
	map<string, vector<segment> > newS;
	if (FH){
		string line;
		int start = 0, stop = 0, N=0;
		string chrom 	= "";
		int i=0,n=0, b=0, u =0, l=0;
		bool collect 	= false;
		vector<segment> current;
		string white_space 	= "                         ";
		string stars 	= "";
		int counter 	= 0;
		while(getline(FH,line)){
			n 	= line.size();
			if (line.substr(0,1)==">"){
				if (!current.empty()){
					S[chrom] 	= current;
					if (test){
						break;
					}
				}
				chrom 	= line.substr(1,line.size());
				start 	= 0, i = 0;
				N 			= S[chrom].size();
				current 	= S[chrom];
				
				if (!chrom.empty()){
					collect 	= true;
				}else{
					collect 	= false;
				}

			}else if (collect){
				while (i < N and current[i].stop < start ){
					i+=1;
				}
				b 	= start + n;
				if (i< N and b > current[i].start){
					l 	= 0;
					for (int k = start; k < b; k++ ){
						u 	= i;
						while (u < N and k > current[u].start){
							if (k <= current[u].stop){
								current[u].seq+=line[l];
							}
							u++;
						}
						l++;
					}

				}
				start+=n;
			}
		}
		if (!current.empty()){
			S[chrom] 	= current;
		}		
		typedef map<string, vector<segment> >::iterator it_type;
		int LOSS 	= 0;
		for (it_type i = S.begin(); i!=S.end(); i++){
			for (int j = 0; j < i->second.size();j++){
				int d 	= i->second[j].stop-i->second[j].start;
				if (i->second[j].seq.size() == d ){
					bool keep 	= i->second[j].transform();
					if (keep){
						newS[i->second[j].chrom].push_back(i->second[j]);
					}else{
						LOSS+=1;
					}
				}else{
					printf("WARNING: ignoring %s:%d-%d, not found in fasta file, %d,%d\n",i->second[j].chrom.c_str(), i->second[j].start, i->second[j].stop,d,i->second[j].seq.size());
				}
			}
		}

	}else{
		printf("couldn't open %s \n", fasta_file.c_str() );
	}
	return newS;
}
PSSM::PSSM(){
	ENRICH_score=0,MD_score=0;
};
PSSM::PSSM(string ID){
	name 	= ID;
	ENRICH_score=0,MD_score=0;
};
PSSM::PSSM(int  id){
	ID 	= id;
	ENRICH_score=0,MD_score=0;
};


string PSSM::get_consensus(){
	string consens 	= "";
	for (int i =0 ; i < frequency_table.size(); i++){
		double max 	= -10000000000;
		int argmax 	= 0;
		for (int j = 0; j < frequency_table[j].size(); j++){
			if (frequency_table[i][j]>max){
				max 	= frequency_table[i][j];
				argmax 	= j;
			}
		}
		consens+=to_string(argmax);
	}
	return consens;
}

double find_closest(double obs, vector<vector<double>> x  ){
	int k;
	int a 	= 0;
	int b 	= x.size();
	int t 	= 0;
	while ((b-a)>2){
		k 	= (b+a)/2;
		if ( obs < x[k][0]   ){
			b 	= k;
		}else{
			a 	= k;
		}
		t+=1;
	}
	return x[k][1];
}


double PSSM::get_pvalue(double obs){
	int k;
	int a 	= 0;
	int b 	= SN;
	while ((b-a)>2){
		k 	= (b+a)/2;
		if ( obs < pvalues[k][0]   ){
			b 	= k;
		}else{
			a 	= k;
		}
	}
	return pvalues[k][1];

}


double PSSM::get_pvalue2_f(double obs, int i, int s){
	int k;
	int a 	= 0;
	int b 	= SN;
	while ((b-a) > 2 ){
		k 	= (b+a)/2;
		if ( obs < position_specific_pvalues_forward[i][k][0]   ){
			b 	= k;
		}else{
			a 	= k;
		}
	}
	return position_specific_pvalues_forward[i][k][1];
}

double PSSM::get_pvalue2_r(double obs, int i, int s){
	int k;
	int a 	= 0;
	int b 	= SN;
	while ((b-a) > 2 ){
		k 	= (b+a)/2;
		if ( obs < position_specific_pvalues_reverse[i][k][0]   ){
			b 	= k;
		}else{
			a 	= k;
		}
	}
	return position_specific_pvalues_reverse[i][k][1];
}

void PSSM::get_pvalue_stats(){
	int i 	= 0;
	int J=0, K = 0;

	while (i < MD_CDF.size()-1){
		if (MD_CDF[i][0] >= MD_score  ){
			J++;
			
		}
		if (MD_CDF[i][0] <=  MD_score ){
			K++;
		}

		
		i++;
	}
	pv_MD_score_rt 	= float(J) / MD_CDF.size();
	pv_MD_score_lt 	= float(K) / MD_CDF.size();
	i 				= 0;
	J=0, K = 0;
	while (i < ENRICH_CDF.size()-1){
		if (ENRICH_CDF[i][0] >= ENRICH_score  ){
			J++;
		}
		if (ENRICH_CDF[i][0] <= ENRICH_score  ){
			K++;
		}
		i++;
	}
	pv_enrich_score_rt 	= float(J) / ENRICH_CDF.size();
	pv_enrich_score_lt 	= float(K) / ENRICH_CDF.size();

}
void PSSM::bin_observations(){
	for (int i = 0 ; i < 2000; i++){
		binned_observed_displacements.push_back(0);
	}
	for (int i = 0 ; i < observed_displacements.size(); i++){
		int x 	= observed_displacements[i];
		binned_observed_displacements[x]++;
	}
}

void PSSM::bin_null_displacements(){
	for (int i = 0 ; i < 2000; i++){
		binned_null_displacements.push_back(0);
	}
	for (int i = 0 ; i < null_displacements.size(); i++){
		for (int j = 0 ; j < null_displacements[i].size(); j++){
			int x 	= null_displacements[i][j];
			binned_null_displacements[x]++;
		}
	}
}







vector<PSSM *> load_PSSM_DB_new(string FILE, int threshold){
	ifstream FH(FILE);
	vector<PSSM *> 	PS;
	PSSM * P 	= NULL;
	bool test  	= threshold > 0;
	if (FH){
		vector<string>line_array;
		string MOTIF 	= "", line 	= "";
		int ID 			= 0, t= 0;

		while(getline(FH,line)){
			if (line.substr(0,1)==">"){
				if (P!=NULL){
					PS.push_back(P);
				}
				line_array 	= split_by_comma(line.substr(1,line.size()-1), " ");				
				if (0<line_array.size() and line_array.size()<3){
					if (test and t > threshold){
						return PS;	
					}
					
					P 	= new PSSM(line_array[0]);
					if (line_array.size()>1){
						P->N 		= stod(line_array[1]);
					}else{
						P->N 		= 0.0;
					}
					t+=1;
				}else{
					P 				= NULL;
				}
			}else if (P!=NULL){
				line_array 	= split_by_comma(line, " ");
				if (line_array.size()==4){
					vector<double> x;
				 	for (int i = 0 ; i < 4; i++){
				 		x.push_back(stod(line_array[i]));
				 	}
				 	P->frequency_table.push_back(x);	
				}else{
					P 	= NULL;
				}

			}
		}
	}else{
		printf("Could not open %s\n",FILE.c_str() );
	}
	if (P!=NULL){
		PS.push_back(P);
	}
	return PS; 
}
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
}


void write_out_null_stats(vector<PSSM *> PSSMS, string OUT, params * PP, 
	vector<double> background, vector<vector<double>> MAP_background,vector<vector<double>> MAP_background_NON,
	 vector<double ** > M,vector<double ** > MN){
	
	string fasta_file 			= PP->p["-fasta"];
	string bed_file 			= PP->p["-bed"];
	string out_dir 				= PP->p["-o"];
	string PSSM_DB 				= PP->p["-DB"];
	string job_ID 				= PP->p["-ID"];
	double window 				= 1000;
	double pv 					= stof(PP->p["-pv"]);
	int test 					= 1;
	int sim_N 					= stoi(PP->p["-sim_N"]);
	int bins 					= stoi(PP->p["-br"]);
	

	ofstream FHW(OUT);
	FHW<<"#=====================================================================\n";
	FHW<<"#Database File\n";
	FHW<<"#Data/Time    "<<currentDateTime()<<endl;
	FHW<<"#-window      "<<to_string(window)<<endl;
	FHW<<"#-bed         "<<(bed_file)<<endl;
	FHW<<"#-fasta       "<<(fasta_file)<<endl;
	FHW<<"#-sim_N       "<<to_string(sim_N)<<endl;
	FHW<<"#-bins        "<<to_string(bins)<<endl;
	FHW<<"#-pv          "<<to_string(pv)<<endl;
	FHW<<"#-background  "<<to_string(background[0])<<","<<to_string(background[1])<<","<<to_string(background[2])<<","<<to_string(background[3])<<endl;
	FHW<<"#=====================================================================\n";
	
	for (int p = 0 ; p < PSSMS.size(); p++){
		FHW<<">" + PSSMS[p]->name <<endl;
		for (int i = 0 ; i < PSSMS[p]->frequency_table.size(); i++){
			string line;
			for (int j = 0 ; j < 3; j++){
				line+=to_string(PSSMS[p]->frequency_table[i][j])+",";
			}
			line+=to_string(PSSMS[p]->frequency_table[i][3])+"\n";
			FHW<<line;
		}
		FHW<<"~";
		//count up the number of zeros
		int zero 	= 1;
		for (int s = 0 ; s < PSSMS[p]->null_displacements.size(); s++){
			if (PSSMS[p]->null_displacements[s].empty()){
				zero++;
			}
		}
		PSSMS[p]->bin_null_displacements();
		FHW<<to_string(zero) + "|" ;

		string line 	= "";

		for (int i = 0 ; i < PSSMS[p]->binned_null_displacements.size() ; i ++  ){
			if (i+1 < PSSMS[p]->binned_null_displacements.size() ){
				line+=to_string(PSSMS[p]->binned_null_displacements[i])+",";
			}else{
				line+=to_string(PSSMS[p]->binned_null_displacements[i]);	
			}
		}

		FHW<<line+"\n";	

	}
	FHW<<"#Estimated Background Distribution (TSS)\n";
	for (int i = 0 ; i < MAP_background.size(); i++){
		FHW<<"#\t"+ to_string(MAP_background[i][0])+","+to_string(MAP_background[i][1])+","+to_string(MAP_background[i][2])+","+to_string(MAP_background[i][3])+ "\n";
	}
	FHW<<"#Estimated Background Distribution (~TSS)\n";
	for (int i = 0 ; i < MAP_background_NON.size(); i++){
		FHW<<"#\t"+ to_string(MAP_background_NON[i][0])+","+to_string(MAP_background_NON[i][1])+","+to_string(MAP_background_NON[i][2])+","+to_string(MAP_background_NON[i][3])+ "\n";
	}
	FHW<<"#Estimate Markov Chain (TSS)\n";
	for (int i = 0 ; i < M.size(); i++){
		FHW<<"#$" + to_string(i)+"\n";
		for (int u = 0 ; u < 4; u++){
			FHW<<"#\t" + to_string(M[i][u][0])+","+to_string(M[i][u][1]) + "," + to_string(M[i][u][2]) + "," + to_string(M[i][u][3])+ "\n";
		}
	}
	FHW<<"#Estimate Markov Chain (~TSS)\n";
	for (int i = 0 ; i < MN.size(); i++){
		FHW<<"#$" + to_string(i)+"\n";
		for (int u = 0 ; u < 4; u++){
			FHW<<"#\t" + to_string(MN[i][u][0])+","+to_string(MN[i][u][1]) + "," + to_string(MN[i][u][2]) + "," + to_string(MN[i][u][3])+ "\n";
		}
	}

}

vector<PSSM *> sort_PSSMS(vector<PSSM *> PSSMS){
	bool switched	= true;
	while (switched){
		switched=false;
		for (int i = 0 ; i < PSSMS.size()-1; i++){
			if (PSSMS[i]->pv_MD_score_rt > PSSMS[i+1]->pv_MD_score_rt){
				PSSM * copy 	= PSSMS[i];
				PSSMS[i] 		= PSSMS[i+1];
				PSSMS[i+1] 		= copy;
				switched 		= true;
			}
		}
	}
	return PSSMS;
}

void write_out_stats(vector<PSSM *> PSSMS, string OUT, params *P){

	ofstream FHW(OUT);
	PSSMS 	= sort_PSSMS(PSSMS);
	FHW<<P->get_header();
	FHW<<"#alpha level 0.01\tadjusted alpha level " + to_string(0.01/PSSMS.size())+"\n";
	FHW<<"#motif identifier\tTotal (2KB)\tMD score\tenrichment score\tMD score p-value (left, right tail)\tenrichment score p-value (left, right tail)\tMD score significance (left, right adj)\tenrichment score significance (left, right adj)\n";
	for (int p = 0 ; p < PSSMS.size(); p++){
		FHW<<PSSMS[p]->name<<"\t"<<to_string(int(PSSMS[p]->total))<<"\t"<<to_string(PSSMS[p]->MD_score)+"\t"+to_string(PSSMS[p]->ENRICH_score)+"\t";
		FHW<<to_string(PSSMS[p]->pv_MD_score_lt)+ "," +to_string(PSSMS[p]->pv_MD_score_rt) +"\t"+to_string(PSSMS[p]->pv_enrich_score_lt)+ "," +to_string(PSSMS[p]->pv_enrich_score_rt)+"\t";
		int MDL	= 0, MDR=0, ENL=0, ENR=0;



		if (PSSMS[p]->pv_MD_score_lt < (0.01/PSSMS.size())){
			MDL=-1;
		}
		if (PSSMS[p]->pv_MD_score_rt <  (0.01/PSSMS.size() )){
			MDR=1;
		}
	
		if (PSSMS[p]->pv_enrich_score_lt < (0.01/PSSMS.size())){
			ENL=-1;
		}
		if (PSSMS[p]->pv_enrich_score_rt <  (0.01/PSSMS.size() )){
			ENR=1;
		}
		FHW<<to_string(MDL) + "," + to_string(MDR) + "\t" + to_string(ENL) + "," + to_string(ENR)+"\n";
	}
	FHW<<"#Binned Observation statistics range={-1000,..,1000}\n";
	for (int p =0 ; p < PSSMS.size(); p++){
		FHW<<PSSMS[p]->name<<"\t";
		PSSMS[p]->bin_observations();
		string line="";
		for (int i = 0; i < PSSMS[p]->binned_observed_displacements.size(); i++ ){
			if (i+1 <  PSSMS[p]->binned_observed_displacements.size()){
				line+=to_string(PSSMS[p]->binned_observed_displacements[i]) + ",";
			}else{
				line+=to_string(PSSMS[p]->binned_observed_displacements[i]);				
			}
		}
		FHW<<line<<endl;

	}


	FHW<<"#Empiracle Bootstrapped Distribution\n";
	FHW<<"#Motif Identifier\tMD Score\tEnrichment Number\n";
	for (int p = 0 ; p < PSSMS.size(); p++){
		string line = "", line2="";

		FHW<<PSSMS[p]->name<<"\t";
		for (int i = 0 ; i < PSSMS[p]->MD_CDF.size();i++){
			if (i+1 < PSSMS[p]->MD_CDF.size()){
				line+=to_string(PSSMS[p]->MD_CDF[i][0])+",";
				line2+=to_string(PSSMS[p]->ENRICH_CDF[i][0])+",";
			}else{
				line+=to_string(PSSMS[p]->MD_CDF[i][0]);	
				line2+=to_string(PSSMS[p]->ENRICH_CDF[i][0]);
			}
		}
		FHW<<line+ "\t" + line2<<endl;
	}	
}

void load_PSSM_ID_names_only(string FILE, map<int, string> & G){
	ifstream FH(FILE);
	if (FH){
		vector<string>line_array;
		int N 	= 0;
		string MOTIF, line;

		while(getline(FH,line)){
			if (line.substr(0,5)=="MOTIF"){
				line_array 	= split_by_ws(line, " ");
				if (line_array.size()>1){
					MOTIF  	= line_array[1];			
					G[N] 	= MOTIF;
					N++;
				}
			}
		}
	}else{

	}
}

vector<PSSM *> convert_streatmed_to_vector(vector<vector<vector<double>>> streamed,vector<int> IDS,
	vector<int> NS){
	int N 	= IDS.size();
	vector<PSSM *> PSSMS;
	for (int i = 0 ; i < N; i++){
		PSSM * current 	= new PSSM(IDS[i]);
		current->frequency_table 	= streamed[i];
		current->N 	= NS[i];
		PSSMS.push_back(current);
	}
	return PSSMS;
}
void collect_all_tmp_files(string dir, string job_name, int nprocs, int job_ID){
	int c 	= 0;
	time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%m_%d_%H_%M", &tstruct);
    string DT 	= buf;
	string OUT 		= dir+ job_name + "-" + to_string(job_ID) +"_" + DT+ ".log";
	ofstream FHW(OUT);
	FHW<<"=======Application Error/Standard Out Log File=======\n\n\n";
	for (int rank = 0; rank < nprocs; rank++){
		string FILE 	= dir+"tmp_" + job_name+ "-" +to_string(job_ID) + "_" + to_string(rank) + ".log";
		string line;
		ifstream FH(FILE);
		if (FH){
			if (rank!=0){
				FHW<<"=======MPI Call: " + to_string(rank) +"=======\n";
			}
			while (getline(FH, line)){
				if ("#" != line.substr(0,1) or rank==0){
					FHW<<line<<endl;
				}
			}

			FH.close();
			remove( FILE.c_str()) ;
			c++;
		}
	}
}


vector<PSSM *> load_personal_DB_file(string FILE, params * P, vector<double> & background){
	vector<PSSM *> PSSMS;
	ifstream FH(FILE);
	PSSM * pssm 	= NULL;
	if (FH){
		string line;
		vector<string> line_array;
		vector<string> line_array_2;
		
		bool header 	= true;
		while (getline(FH,line)){
			if (line.substr(0,1)=="#" and not header){
				//fill in P
				line_array 	= split_by_ws(line.substr(1,line.size()-1), "");
				string F 	= line_array[0];
				string val 	= line_array[line_array.size()-1];
				if (F.substr(0,F.size())=="-pv" or F.substr(0,F.size())=="-bins"){
					P->p[F] = val;
				}else if (F.substr(0,F.size())=="-background"){
					line_array 	= split_by_comma(val, "");
					for (int i = 0 ; i < 4;i++){
						background.push_back(stod(line_array[i]));
					}
				}

			}else if (line.substr(0,1)==">" and line.substr(0,1)!="#"){
				if (pssm != NULL){
					PSSMS.push_back(pssm);
				}
				pssm 	= new PSSM(line.substr(1,line.size()-1));
			}else if(pssm != NULL and line.substr(0,1)!="~" and line.substr(0,1)!="#"){
				line_array 	= split_by_comma(line, "");
				vector<double> row;
				for (int i = 0 ; i < 4;i++){
					row.push_back(stod(line_array[i]));
				}
				pssm->frequency_table.push_back(row);
				
			}else if(line.substr(0,1)=="~" and pssm!=NULL and line.substr(0,1)!="#"){
				line_array 		= split_by_bar(line.substr(1,line.size()-1), " ");
				double zeros 	= stod(line_array[0]);
				line_array 		= split_by_comma(line_array[1], " ");
				for (int i = 0 ; i < line_array.size(); i++){
					int x 		= stoi(line_array[i]);
					pssm->null_displacements_2.push_back(x);						
				}
				pssm->zeros 	= zeros;
			}
			header 		= false;

		}

	}else{
		printf("\n\nCould not open %s\n",FILE.c_str() );
	}
	if (pssm!=NULL){
		PSSMS.push_back(pssm);
	}
	return PSSMS;
}

vector<segment> merge(vector<segment> S, string chr){
	vector<segment> newS;
	int i = 0, o_st = 0, o_sp = 0;
	int ID 	= 0;
	while (i < S.size()){
		o_st 	= S[i].start , o_sp = S[i].stop;
		while (i < S.size() and S[i].start < o_sp and S[i].stop > o_st ){
			o_st 	= min(S[i].start, o_st) , o_sp 	= max(S[i].stop, o_sp);
			i++;
		}
		segment nS(chr, o_st, o_sp, ID, o_st, o_sp);
		ID++;
		newS.push_back(nS);
	}
	return newS;
}




map<string, vector<segment>>  label_TSS(map<string, vector<segment>> S, map<string, vector<segment>> TSS){
	//need to sort both
	typedef map<string, vector<segment>>::iterator it_type;
	//merge and sort TSS
	for (it_type c = TSS.begin(); c!=TSS.end(); c++){
		TSS[c->first] = merge(c->second, c->first);
	}	
	for (it_type c = S.begin(); c!=S.end(); c++){
		if (TSS.find(c->first)!=TSS.end()){
			int j = 0, N = TSS[c->first].size();
			for (int i = 0 ; i < c->second.size(); i++){
				while (j < N and TSS[c->first][j].stop < c->second[i].start){
					j++;
				}
				if (j < N and TSS[c->first][j].start < c->second[i].stop ){
					c->second[i].TSS 	= true;
				}
			}
		}
		S[c->first] 	= c->second;
	}
	int tss = 0 , NOT = 0;
	for (it_type c = S.begin(); c!=S.end(); c++){
		for (int i = 0 ; i < c->second.size(); i++){
			if (c->second[i].TSS){
				tss++;
			}else{
				NOT++;
			}
		}
	}
	return S;

}








