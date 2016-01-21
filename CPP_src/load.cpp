#include "load.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "split.h"
using namespace std;

segment::segment(){};
segment::segment(string chr, int st, int sp){
	chrom=chr, start=st, stop=sp;
	seq 	= "";
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
	forward 	= new int[seq.size()];
	reverse 	= new int[seq.size()];
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



map<string, vector<segment>> load_bed_file(string FILE, int pad){
	map<string, vector<segment>> S ;
	ifstream FH(FILE);

	if (FH){
		string line, chrom;
		vector<string> line_array;
		int start, stop;
		double x;
		while (getline(FH, line)){
			if (line.substr(0,1)!="#"){
				line_array=splitter(line, "\t");
				chrom 	= line_array[0], start = stoi(line_array[1]), stop = stoi(line_array[2]);
				if (pad != 0){
					x 	= (start + stop)/2.;
					S[chrom].push_back(segment(chrom, x-pad, x+pad));
				}else{
					S[chrom].push_back(segment(chrom, start, stop));					
				}
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
map<string, vector<segment> > insert_fasta_sequence(string fasta_file, map<string, vector<segment> > S){
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
					
				//	break;
				}
				counter+=1;
				stars+="*";
				printf("\rinserting fasta|%s%s|", stars.c_str(), white_space.substr(0,white_space.size()-counter).c_str());
				chrom 	= line.substr(1,line.size());
				cout.flush();
				start 	= 0, i = 0;
				N 	= S[chrom].size();
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
		printf("\rinserting fasta|%s%s|done\n", stars.c_str(),  "");
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
					printf("WARNING: ignoring %s:%d-%d, not found in fasta file, %d\n",i->second[j].chrom.c_str(), i->second[j].start, i->second[j].stop,d);
				}
			}
		}

	}else{
		printf("couldn't open %s \n", fasta_file.c_str() );
	}
	return newS;
}
PSSM::PSSM(){};
PSSM::PSSM(string ID){
	name 	= ID;
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

vector<PSSM *> load_PSSM_DB(string FILE){
	ifstream FH(FILE);
	vector<PSSM *> all_motifs;
	if (FH){
		string line;
		vector<string>line_array;
		string MOTIF 	= "";
		PSSM * P 	= NULL;
		while(getline(FH,line)){

			if (line.substr(0,5)=="MOTIF"){
				MOTIF 		= "";
				P 			= NULL;
				line_array 	= split_by_ws(line, " ");
				if (line_array.size()>1){
					line_array 	= split_by_us(line_array[1], " ");
					MOTIF  	= line_array[0];			
					P 		= new PSSM(MOTIF);
				}

			}else if (line.substr(0,6)=="letter" and P!=NULL){
				line_array 	= split_by_ws(line, " ");
				if (line_array.size()>7){
					P->N 	= stoi(line_array[7]);
				}else{
					P 		= NULL;
				}
			}else if(P!=NULL and line.substr(0,3)!="URL"){
				line_array 	= split_by_tab(line, " ");

				if (line_array.size() == 4){
				 	vector<double> x;
				 	for (int i = 0 ; i < 4; i++){
				 		x.push_back(stod(line_array[i]));
				 	}
				 	P->frequency_table.push_back(x);
				}

			}else if(line.substr(0,3)=="URL" and P!= NULL){
				all_motifs.push_back(P);
				P 	= NULL;
			}
		}	
	}else{
		printf("couldn't open %s\n",FILE.c_str() );
	}
	return all_motifs;
}





