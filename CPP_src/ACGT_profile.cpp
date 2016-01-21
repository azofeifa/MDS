#include "ACGT_profile.h"

using namespace std;
map<string, vector<vector<double> > > get_average_ACGT_profile(map<string, vector<segment> > S, 
	vector<PSSM *> PSSMS, int pad ){
	
	map<string, vector<vector<double> > > G;
	typedef map<string, vector<segment> >::iterator it_type;
	typedef vector<segment>::iterator it_type_2;
	for (int p = 0; p < PSSMS.size(); p++){
		for (int i = 0; i < pad*2; i++){
			vector<double> current 	= {1.0,1.0,1.0,1.0};
			G[PSSMS[p]->name].push_back(current);
		}
		for (it_type c = S.begin(); c!=S.end(); c++){
			for (it_type_2 s=c->second.begin(); s!=c->second.end(); s++){
				if (  s->motif_positions[p].empty() ){
					for (int i = 0 ; i <s->N ; i++ ){
						G[PSSMS[p]->name][i][s->forward[i]]+=1;
					}
				}
			}
		}
		for (int i = 0 ; i <G[PSSMS[p]->name].size(); i++ ){
			double S = 0.0;
			for (int j = 0 ; j <G[PSSMS[p]->name][i].size(); j++ ){
				S+=G[PSSMS[p]->name][i][j];
			}
			for (int j = 0 ; j <G[PSSMS[p]->name][i].size(); j++ ){
				G[PSSMS[p]->name][i][j]/=S;
			}
		}
	}



	return G;
}