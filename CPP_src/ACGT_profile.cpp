#include "ACGT_profile.h"

using namespace std;
map<int, double [2000][4]> get_average_ACGT_profile(map<string, vector<segment> > S, 
	vector<PSSM *> PSSMS, int pad,map<int,   double>   & NN , 
	map<int, double [2000][4]> & G, int rank )
{
	typedef map<string, vector<segment> >::iterator it_type;
	typedef vector<segment>::iterator it_type_2;
	double current[2000][4];
	for (int p = 0; p < PSSMS.size(); p++){
		for (int i = 0 ; i < 2000; i++ ){
			current[i][0]=1.0,current[i][1]=1.0,current[i][2]=1.0,current[i][3]=1.0;
		}
		NN[PSSMS[p]->ID] 			= 0.0;
		bool FOUND 					= false;
		for (it_type c = S.begin(); c!=S.end(); c++){
			for (it_type_2 s=c->second.begin(); s!=c->second.end(); s++){
				if ( !s->motif_positions[PSSMS[p]->ID].empty() and s->N == 2000 ){
					FOUND 	= true;
					for (int i = 0 ; i <2000 ; i++ ){
						current[i][s->forward[i]]=1+current[i][s->forward[i]];
					}
					NN[PSSMS[p]->ID]+=1;
				}
				
			}
						

		}
		for (int i = 0 ; i <2000; i++ ){
			double S = 0.0;
			for (int j = 0 ; j <4; j++ ){
				S+=current[i][j];
			}
			for (int j = 0 ; j <4; j++ ){
				current[i][j]=current[i][j]/S;
			}
		}
		for (int i = 0 ; i < 2000; i++){
			for (int j = 0; j < 4; j++){
				G[PSSMS[p]->ID][i][j] 	= current[i][j];
			}
		}
	}
	return G;
}

void get_ACGT_profile_all(map<string, vector<segment> > S, 
	vector<vector<double>> & background_forward,vector<vector<double>> & background_reverse, int rank){
	for (int i = 0 ; i < 2000; i++){
		vector<double> current 	= {10,10,10,10};
		background_forward.push_back(current);
		background_reverse.push_back(current);
		
	}

	typedef map<string, vector<segment> >::iterator it_type;
	typedef vector<segment>::iterator it_type_2;
	double N 	= 0;
	for (it_type c = S.begin(); c!=S.end(); c++){
		for (it_type_2 s=c->second.begin(); s!=c->second.end(); s++){
			for (int i = 0 ; i <2000 ; i++ ){
				background_forward[i][s->forward[i]]++;
				background_reverse[i][s->reverse[i]]++;
			}
		}
	}
	for (int i = 0 ; i < 2000; i ++){
		double SUM_f= 0, SUM_r=0;
		for (int j = 0; j < 4; j++ ){
			SUM_f+=background_forward[i][j];
			SUM_r+=background_reverse[i][j];
		}
		for (int j = 0; j < 4; j++ ){
			background_forward[i][j]/=SUM_f;
			background_reverse[i][j]/=SUM_f;
		}
	}



}