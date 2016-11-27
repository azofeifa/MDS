#include "collect_sample_statistics.h"
#include <cmath>
#include <random>
#include <iostream>
using namespace std;
void fill_displacements(map<int, vector<double>> & displacements, vector<segment> segments ){
  typedef map<int, vector<double>>::iterator it_type;
  for (int s = 0; s < segments.size(); s++){
    for (it_type p = segments[s].motif_positions.begin();p!=segments[s].motif_positions.end(); p++ ){
      if (p->second.size()){
	for (int i = 0; i < p->second.size(); i++){
	  displacements[p->first].push_back(p->second[i]);
	}
      }else if (displacements.find(p->first)==displacements.end() ) {
	displacements[p->first]={};
      }
    }
  }
}

void get_stats(vector<double> displacements, vector<double> & stats    ){
  double mean=0, var=0, se=0, N=0;
  for (int i = 0 ; i < displacements.size(); i++ ){
    double d 	= abs(displacements[i])-1000;
    mean+=abs(d);
    var+=pow(d,2);
    if (abs(d) <100){
      se+=1;
    }
    N+=1;
  }
  if (N>0){
    stats.push_back(se/N);
    stats.push_back(mean/N);
    stats.push_back(sqrt(var/N));
    stats.push_back(N);
  }else{
    stats.push_back(0);
    stats.push_back(0);
    stats.push_back(0);
    stats.push_back(N);	
  }
}

double get_min_2(vector<double> X){
  if (X.empty()){
    return 3000;
  }
  double MIN 	= abs(X[0]-1000);
  for (int xx = 1 ; xx < X.size(); xx++){
    if (abs(X[xx]-1000) < MIN) {
      MIN 	= abs(X[xx]-1000);
    }
  }
  return MIN;
}



/*
void collect_sample_stats(map<string, vector<segment>> observed,
			  vector<PSSM *> P,  
			  map<int, vector<double> > & observed_statistics,
			  map<int, vector<double>> & observed_displacements,map<int, map<int, double> > & observed_co_occurences,
			  int rank ){

  //through back three things
  //(3) -> map[motif_ID] 	= binned histogram of displacements?
  vector<segment> combinded;
  typedef map<string, vector<segment>>::iterator it_type; //just go ahead and collapse down the segments from chromosomes
  typedef map<int, map<int, vector<segment> >>::iterator it_type_2; //just go ahead and collapse down the segments from chromosomes
  typedef  map<int, vector<segment> >::iterator it_type_3; //just go ahead and collapse down the segments from chromosomes
  typedef  map<int, vector<double> >::iterator it_type_4; //just go ahead and collapse down the segments from chromosomes
  typedef  map<int ,map<int, vector<double> >>::iterator it_type_5; //just go ahead and collapse down the segments from chromosomes
  for (it_type i = observed.begin(); i!=observed.end(); i++){
    for (int j = 0 ; j < i->second.size(); j++){
      combinded.push_back(i->second[j]);
    }
  }
  for (int p1 = 0 ; p1 < P.size(); p1++)	{
    for (int p2= 0; p2 < P.size(); p2++){
      observed_co_occurences[p1][p2] 	= 0;
    }
  }
  map<int,int> norms;
  for (int c = 0 ; c < combinded.size(); c++){
    for (it_type_4 p1 = combinded[c].motif_positions.begin(); p1!=combinded[c].motif_positions.end(); p1++){
      if ( get_min_2(p1->second) < 100 ){
	for (it_type_4 p2 = combinded[c].motif_positions.begin(); p2!=combinded[c].motif_positions.end(); p2++){
	  if (get_min_2(p2->second) < 100){
	    observed_co_occurences[p1->first][p2->first]+=1;
	  }
	}
	if (norms.find(p1->first)==norms.end()){
	  norms[p1->first]=0;
	}
	norms[p1->first]+=1;
      }

    }
  }
  typedef map<int, map<int, double> >::iterator it_type_7;
  typedef  map<int, double> ::iterator it_type_8;

  for (it_type_7 p1 = observed_co_occurences.begin(); p1!=observed_co_occurences.end(); p1++){
    for (it_type_8 p2 = p1->second.begin(); p2!=p1->second.end(); p2++){
      if (norms.find(p1->first)!=norms.end()){
	observed_co_occurences[p1->first][p2->first]/=	float(norms[p1->first]);
      }
    }		
  }	


  //this is PSSM_ID->to chuck->to displacements (to_chunk number gives you N_bidir)
  fill_displacements(observed_displacements, combinded);
	
  for (it_type_4 i = observed_displacements.begin(); i!=observed_displacements.end();i++){
    vector<double> current_stats;
    get_stats(i->second, current_stats);
    observed_statistics[i->first]=current_stats;
  }
}
*/
double get_MD_score(vector<int> D, int window, bool NORM, int large_window){
  double N 	= 0, S 	= 0, center 	= large_window; //this is assumed;NO ITS NOT!!!
  double a 	= center-window, b = center + window;
  for (int i = 0 ; i < D.size();i++){
    if ( D[i] >=  a and D[i] <= b){
      S++;
    }
    N++;
	
  }
  if (N>0 and NORM){
    return S/ N;
  }
  return S;
}

vector<double> get_many_MD_scores(vector<int> D, int step_size, int large_window){
  vector<double> md_scores;
  for (int window = 1 ; window < large_window; window+=step_size){
    md_scores.push_back(get_MD_score(D, window, 1, large_window));
  }
  return md_scores;

}




vector<double> bubble_sort(vector<double> X){
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



vector<vector<double> > make_CDF(vector<double> X){
  X 	= bubble_sort(X);
  double N 	= X.size();
  vector<vector<double>> Y;
  for (int i = 0 ; i < X.size(); i++) {
    vector<double> x 	= {X[i], i/N};
    Y.push_back(x);
  }
  return Y;
}

void turn_to_CDF(vector<double> & CDF, vector<int> X, double & N){
  for (int i = 0 ; i < X.size(); i++){
    N+=X[i];
  }
  double new_n=0;
  for (int i = 0 ; i < X.size(); i++){
    new_n+=X[i];
    CDF.push_back(X[i]/N);
  }
}
int get_sample(vector<double> CDF,double U ){
  int i=0;
  while (i+ 1 < CDF.size() and U > CDF[i]){
    i++;
  }
  return i; 
}





void build_cdfs_PSSMs(PSSM *  P, int bsn, int interval_size, int hit_size, int hit_size_TSS, 
		      int hit_size_NON, int MD_window, double TSS_association, int large_window){
  random_device rd;
  // Initialize Mersenne Twister pseudo-random number generator
  mt19937 gen(rd());
	
	
  double total_N 		= 0;
  double total_N_non 	= 0;

  vector<double> CDF;
  vector<double> CDF_non;

  turn_to_CDF(CDF,P->null_displacements_2, total_N);
  turn_to_CDF(CDF_non,P->null_displacements_2_non, total_N_non);


  discrete_distribution<int> distribution(CDF.begin(),CDF.end());
  discrete_distribution<int> distribution_non(CDF_non.begin(),CDF_non.end());


  uniform_real_distribution<double> distribution_2(0,1);


  vector<double> MD_scores;
  vector<double> MD_scores_NON;
  vector<double> MD_scores_TSS;



  for (int b = 0 ; b < bsn; b++){//make this number of random collections
    double enriched 	= 0;
    vector<int> current_collection_spec;
    vector<int> current_collection_spec_NON;
    vector<int> current_collection_spec_TSS;
    int spec_N 	= 0;
    int k 		= 0;
    int t 		= 0;
    double NN 	= 0.0;
    double U;
    for (int i = 0 ; i < hit_size; i++ ){


      U 	= distribution_2(gen);
      if (U < TSS_association){
	k 	= distribution(gen);
	if ( i < hit_size_TSS){
	  current_collection_spec_TSS.push_back(k);
	}
      }else{
	k 	= distribution_non(gen);
	if (i < hit_size_NON){
	  current_collection_spec_NON.push_back(k);
	}
      }
      current_collection_spec.push_back(k);				
    }
    double MD_score 		= get_MD_score(current_collection_spec, MD_window,true,large_window);
    double MD_score_TSS 	= get_MD_score(current_collection_spec_TSS, MD_window,true,large_window);
    double MD_score_NON 	= get_MD_score(current_collection_spec_NON, MD_window,true, large_window);
		
    MD_scores.push_back(MD_score);
    MD_scores_TSS.push_back(MD_score_TSS);
    MD_scores_NON.push_back(MD_score_NON);

  }
  P->MD_CDF 			= make_CDF(MD_scores);
  P->MD_CDF_TSS 		= make_CDF(MD_scores_TSS);
  P->MD_CDF_NON 		= make_CDF(MD_scores_NON);
}




















