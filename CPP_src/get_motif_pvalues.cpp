#include "get_motif_pvalues.h"
#include <iostream>
#include <cmath>
#include <algorithm>
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

void histogram(vector<double> x, vector<double> y, int bins, double min_x, double max_x, vector<double>& edges,vector<double>& counts ){
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
	smooth_frequence_table(p);
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


vector<PSSM *> DP_pvalues(vector<PSSM *> P, int bins,vector<double> background){
	double current 				= 0.0;
	string stars 				= "";
	string white_space 			= "";
	double val 					= 0.1;
	int  ws_counter 			= 1.0/val;
	for (int i = 0 ; i < P.size(); i++){
		vector<vector<double>> p_values 	= compute_pvalues(P[i], background,200);
		for (int k = 0 ; k < 200; k++){
			for (int j = 0 ; j < 2;j++){
				P[i]->pvalues[k][j]  	= p_values[k][j];
			}
		}

		P[i]->SN 							= p_values.size();
		if (double(i)/P.size() >= current){
			white_space 	= "";
			for (int d = 0; d < ws_counter; d++){
				white_space+=" ";
			}
			ws_counter-=1;

			current+=val;
			stars 		+= "*";
			printf("\rcalculating LLRs|%s%s|", stars.c_str(), white_space.c_str());
			cout.flush();
		}
	}
	stars 		+= "*";
	printf("\rcalculating LLRs|%s%s| done\n", stars.c_str(), "");
	cout.flush();
	return P;

}
