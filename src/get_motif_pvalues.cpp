#include "get_motif_pvalues.h"
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <omp.h>
double get_order_stat(vector<vector<double>> X, vector<double> background, bool MIN){
	double current;
	double val 	= 0;
	for (int i = 0 ; i < X.size(); i++){
		if (MIN){
			current 	= log(pow(10,1000));
		}else{
			current 	= log(pow(10,-1000));
		}
		for (int j = 0; j < X[i].size(); j++){
			if (MIN){
				current = min(current, (log(X[i][j]) - log(background[j])));
			}else{
				current = max(current, (log(X[i][j]) - log(background[j])));	
			}
		}
		val+=current;
	}
	return 2*val;
}
void smooth_frequence_table(PSSM * p){
	vector<vector<double>> X 	= p->frequency_table;
	for (int i = 0 ; i < X.size(); i++){
		double SUM 	= 0.0;
		for (int j = 0; j < X[i].size(); j++){
			X[i][j]=(X[i][j]*double(p->N) + 1);
			SUM+=X[i][j];
		}
		for (int j = 0; j < X[i].size(); j++){
			X[i][j]=X[i][j]/SUM;			
		}
	}
	p->frequency_table 	= X;
}

void sort_xy(vector<double>& x, vector<double>& y){
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < x.size(); i++){
			if (x[i-1]  > x[i]){
				changed 		= true;
				double X 	= x[i-1];
				double Y 	= y[i-1];
				x[i-1] 	= x[i];
				x[i] 	= X;
				y[i-1] 	= y[i];
				y[i] 	= Y;
				
			}
		}
	}

}

void histogram(vector<double> x, vector<double> y, int bins, 
	double min_x, double max_x, vector<double>& edges,vector<double>& counts )
{
	double delta 	= (max_x - min_x) / double(bins);
	sort_xy(x,y);

	for (int i = 0; i < bins; i++){
		edges.push_back(min_x + (i*delta));
		counts.push_back(0);
	}
	int j = 0;
	int N = edges.size();
	for (int i = 0 ; i < x.size(); i++){
		while ((j+1) < N and edges[j] <= x[i]){
			j++;
		}
		if (j > 0){
			counts[j-1]+=y[i];
		}
	}
}

vector<vector<double>> compute_pvalues(PSSM * p, vector<double > background,int bins){
	vector<vector<double>>p_values;
	int N 	= p->frequency_table.size(); //number of positions in the PSSM
	double min_score 	= get_order_stat(p->frequency_table, background, true);
	double max_score 	= get_order_stat(p->frequency_table, background, false);

	vector<vector<double>> A 	= p->frequency_table;
	vector<double> counts,edges;
	for (int i = 0 ; i < N; i++){
		vector<double> nvls, weights;

		for (int j = 0; j < A[i].size(); j++){
			if (i==0){
				nvls.push_back(2*(log(A[i][j]) - log(background[j])  ));
				weights.push_back(1);

			}else{
				for (int k = 0; k < edges.size(); k++){
					if (counts[k]>0){
						nvls.push_back(edges[k]+2*(log(A[i][j]) - log(background[j])  ));
						weights.push_back(counts[k]);
					}
				}
			}
		}
		edges.clear();
		background.clear();
		histogram(nvls, weights, bins, min_score, max_score, edges, counts);
	}
	sort_xy(edges, counts);
	double S 	= 0.0;
	double NN 	= 0.0;
	for (int i = 0 ; i < edges.size(); i++){
		NN+=counts[i];
	}
	for (int i = 0 ; i < edges.size(); i++){
		S+=counts[i];
		vector<double> current(2);
		current[0] 	= edges[i], current[1] = double(S)/double(NN);
		p_values.push_back(current);
	}

	return p_values;
}


vector<PSSM *> DP_pvalues(vector<PSSM *> P, int bins,
		vector<double> background, bool SMOOTH, double pval){
	double current 				= 0.0;
	string stars 				= "";
	string white_space 			= "";
	double val 					= 0.1;
	int  ws_counter 			= 1.0/val;
	#pragma omp parallel for
	for (int i = 0 ; i < P.size(); i++){
		if (SMOOTH){
			smooth_frequence_table(P[i]);
		}
	
		int BINS 							= P[i]->frequency_table.size()*bins;
		vector<vector<double>> p_values 	= compute_pvalues(P[i], background,BINS);

		for (int k = 0 ; k < BINS; k++){
			P[i]->pvalues.push_back(vector<double>(2));
			for (int j = 0 ; j < 2;j++){
				P[i]->pvalues[k][j]  	= p_values[k][j];
			}
		}
		P[i]->SN 							= p_values.size();
	}

	#pragma omp parallel for
	for (int p = 0 ; p < P.size();p++){
		for (int s = 0 ; s < P[p]->frequency_table.size();s++){
			for (int i = 0 ; i < 4; i++){
				P[p]->frequency_table[s][i] 	= 2*(log(P[p]->frequency_table[s][i]) - log(background[i]));
			}
		}
		P[p]->get_ll_threshold(pval);
	}


	return P;
}

int get_random_multinomial(vector<double> frequencies, double U ){
	int I 			= 0;
	while (I < frequencies.size() and U > frequencies[I]){
		I+=1;
	}
	if (I==4){
		I=3;
	}
	return I;
}

double get_score(vector<vector<double>> frequency_table, vector<int> obs ){
	double SCORE 	= 0;
	for (int o = 0 ; o < obs.size(); o++){
		SCORE+=log(frequency_table[o][obs[o]]);
	}
	return SCORE;
}

vector<double> bs_sort(vector<double> X){
	bool changed 	= true;
	while (changed){
		changed 	= false;
		for (int i = 1 ; i < X.size(); i++){
			if (X[i-1] > X[i] ){
				double cp 	= X[i];
				X[i] 		= X[i-1];
				X[i-1] 		= cp;
				changed 	= true;
			}
		}
	}
	return X;

}




vector<vector<double>> BIN(vector<double> X, int bins){
	vector<vector<double>> BINNED ;
//	X 				= bs_sort(X);
	double min_x 	= X[0];
	double max_x 	= X[X.size()-1];
	double step 	= (max_x-min_x) / bins;
	int j 			= 0;
	for (int b = 0 ; b < bins; b++){
		vector<double> current(2);
		current[0] 	= min_x + step*b;
		current[1] 	= 0.0;
		BINNED.push_back(current);
	}
	for (int i = 0 ; i< X.size(); i++){
		while ( j < BINNED.size() and BINNED[j][0] < X[i] ){
			j++;
		}
		if (j < BINNED.size()){
			BINNED[j][1]+=1;
		}else{
			BINNED[j-1][1]+=1;	
		}
	}
	for (int b = 0 ; b < BINNED.size();b++){
		BINNED[b][1]=BINNED[b][1]/X.size();
	}
	double S 	= 0.0;
	for (int b = 0 ; b < BINNED.size(); b++){
		S 				+= 	BINNED[b][1];
		BINNED[b][1] 	=  	S;
	}
	return BINNED;
}

void insert(double x, vector<double> & array){
	int a 	= 0, b = array.size();
	int i 	= b+a/2, previ 	= -1000;
	if (i){
		while (true){ 
			i 	= (b+a)/2;
			if (i==previ){
				break;
			}
			if (x >= array[i]){
				a 	= i;
			}else if (x <= array[i]){
				b 	= i;
			}
			previ 	= i;
		}
		i 	= (b+a)/2;
		if (x > array[i]){
			i++;
		}
		array.insert(array.begin()+i, x);
	}else{
		array.push_back(x);
	}
	
}


void compute_pvalues_simulation(PSSM * p, vector<vector<double>> & background, int bins, int site_br, int s){
	smooth_frequence_table(p);

	map<int,vector<vector<double>> > binned_positions_specific;
	int delta 		= background.size()/site_br;


	for (int b = 0 ; b < site_br; b++){
		int st 	= b*delta, sp = (b+1)*delta;
		vector<double> background_current(4);
		double NN 				= 0.0;
		background_current[0] 	= 0,background_current[1] 	= 0,background_current[2] 	= 0,background_current[3] 	= 0;
		for (int i = st ; i < sp; i++ ){

			for (int s = 0 ; s < 4; s++){
				background_current[s]+=background[i][s];
			}
			NN++;
		}
		for (int s = 0 ; s < 4; s++){
			background_current[s]/=NN;
		}
		for (int i = st ; i < sp ; i++)	{
			background[i] 						= background_current;
		}
		
		binned_positions_specific[b] 			= compute_pvalues(p, background_current, bins);

	}
	int current_stop 	= delta;
	int B 				= 0;
	for (int b = 0 ; b < background.size(); b++){
		int BINS 							= p->frequency_table.size()*bins;
		if (b > delta){
			delta+=delta;
			B++;
		}


		if (s==1){
			p->position_specific_pvalues_forward[b] = binned_positions_specific[B];
		}else{
			p->position_specific_pvalues_reverse[b] = binned_positions_specific[B];
		
		}
		p->SN 							= p->position_specific_pvalues_forward[b].size();		
	}
}


vector<PSSM *> construct_position_specific_pvalues(vector<PSSM * > P, int bins, 
	vector<vector<double>> & background_forward, vector<vector<double>> & background_reverse, int site_br){
	//smooth background
	for (int b = 0 ; b < background_forward.size(); b++){
		vector<double> background_current(4);
		background_current[0] 	= 0,background_current[1] 	= 0,background_current[2] 	= 0,background_current[3] 	= 0;
		vector<double> background_current_r(4);
		background_current_r[0] 	= 0,background_current_r[1] 	= 0,background_current_r[2] 	= 0,background_current_r[3] 	= 0;
		double N 	= 0;
		int SP 		= b+20;
		if ((b + 20) > background_forward.size()){
			SP  	= background_forward.size();
		}
		for (int j = b; j < SP; j++){
			for (int s = 0 ; s < 4; s++){
				background_current[s]+=background_forward[j][s];
				background_current_r[s]+=background_reverse[j][s];
			}
			N++;
		}
		for (int s = 0 ; s < 4; s++){
			background_current[s]/=N;
			background_current_r[s]/=N;
		}	
		background_forward[b] 	= background_current;
		background_reverse[b] 	= background_current_r;

	}


	#pragma omp parallel for
	for (int i = 0 ; i < P.size(); i++){
		compute_pvalues_simulation(P[i], background_forward,  P[i]->frequency_table.size()*bins, site_br,1);
		compute_pvalues_simulation(P[i], background_reverse, P[i]->frequency_table.size()*bins, site_br,-1);
	}

	return P;

}
