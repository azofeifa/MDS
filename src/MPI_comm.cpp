#include "MPI_comm.h"
#include <mpi.h>
using namespace std;

double * allocated_PSSM_array(int rows, int cols, PSSM * p){
	double * data 	= new double[rows*4];
	int k 			= 0;
	for (int r = 0; r < rows; r++){
		for (int c = 0; c < 4; c++){
			data[k] 	= p->frequency_table[r][c];
			k++;
		}
	}
	return data;
}


void test_send_PSSMS(int rank, int nprocs, vector<PSSM *> PSSMS, 
	vector<vector<vector<double>>> & streamed_PSSMS, 
	vector<int>  & PSSM_IDS, vector<int> & PSSM_N){
	int S;
	if (rank==0){
		int start, stop;
		int N 		= PSSMS.size();
		int count 	= N / nprocs;
		if (count == 0){
			count 	= 1;
		}
		
		for (int j =0; j < nprocs; j++){
			start 	= j*count;
			stop 	= (j+1)*count;
			stop 		= min(stop, N);
			if (j == nprocs-1){
				stop 	= N;
			}
			if (start > stop){
				stop 	= start; //in this case there are more nodes than segments...so
			}
			S 		= stop-start;
			if (j > 0){
				MPI_Send(&S, 1, MPI_INT, j,1, MPI_COMM_WORLD);
				int k 	= 0;
				for (int i = start ; i < stop; i++){
					int PN 			= PSSMS[i]->frequency_table.size();
					int sN 			= PSSMS[i]->N;
					double * P 	= allocated_PSSM_array(PN, 4, PSSMS[i]);
					MPI_Send(&sN, 1, MPI_INT, j,k+3*S, MPI_COMM_WORLD);
					MPI_Send(&i, 1, MPI_INT, j,k+2*S, MPI_COMM_WORLD);
					MPI_Send(&PN, 1, MPI_INT, j,k+S, MPI_COMM_WORLD);
					MPI_Send(&P[0], 
						PSSMS[i]->frequency_table.size()*4, MPI_DOUBLE, j, k, MPI_COMM_WORLD);
					k++;
				}
			}else{
				for (int i = start; i < stop;i++){
					streamed_PSSMS.push_back(PSSMS[i]->frequency_table);
					PSSM_IDS.push_back(i);
					PSSM_N.push_back(PSSMS[i]->N);
				}
			}
		}		

	}else{
		MPI_Recv(&S, 1, MPI_INT, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		for (int i = 0 ; i < S; i++){
			int sN;
			MPI_Recv(&sN, 1, MPI_INT, 0, i+3*S, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			
			int ID;
			MPI_Recv(&ID, 1, MPI_INT, 0, i+2*S, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			int PN;
			MPI_Recv(&PN, 1, MPI_INT, 0, i+S, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			double * P 	= new double[PN*4];
			MPI_Recv(&P[0], PN*4, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//convert back to 2d
			int k 	= 0;
			vector<vector<double>> new_PSSM;
			for (int r=0; r < PN; r++){
				vector<double> current(4);
				for (int c = 0; c < 4; c++){
					current[c] 	= P[k];
					k++;
				}
				new_PSSM.push_back(current);
			}
			streamed_PSSMS.push_back(new_PSSM);
			PSSM_IDS.push_back(ID);
			PSSM_N.push_back(sN);

		}
	
	}
}


void make_array_1D(vector<double> observed_statistics, double * info ){
	for (int i = 0; i < observed_statistics.size(); i++){
		info[i] 	= observed_statistics[i];
	}
}

void make_array_2D(vector<vector<double>> null_statistics, double * info ){
	int k = 0;
	for (int j = 0 ; j < 4; j++){
	
		for (int i = 0; i < null_statistics.size(); i++){
			info[k] 	= null_statistics[i][j];
			k++;		
		}
	}
}

vector<double> transform_back(double * X, int N){
	vector<double> x;
	for (int i = 0 ; i < N; i++){
		x.push_back(X[i]);
	}
	return x;
}



map<int, vector< vector <double> >> collect_PSSM_hits(int rank, int nprocs, 
	map<string, vector<segment>> intervals,
	map<int, vector<double> > observed_statistics,
	map<int, vector<vector<double> >> observed_null_statistics,
	map<int, vector<double>>  observed_displacements	   ){


	map<int, vector< vector <double> >> collections;

	if (rank==0){//will be doing all the receiving
		for (int j = 1; j < nprocs; j++ ){
			int S;
			int tag 	= 0;
			MPI_Recv(&S, 1, MPI_INT, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
			tag++;
			for (int p =0; p < S; p++){ //S here is the number of PSSMS
				int * info 	= new int[4];//info[0]=PSSM ID, info[1] number of observed statistics, info[2] number of flattened null statitics

				MPI_Recv(&info[0], 4, MPI_INT, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				tag++;
				double * obs_stats 	= new double[info[1]];
				double * null_stats = new double[info[2]];
				double * displacements = new double[info[3]];
				
				MPI_Recv(&obs_stats[0], info[1], MPI_DOUBLE, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
				tag++;
				MPI_Recv(&null_stats[0], info[2], MPI_DOUBLE, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
				tag++;
				MPI_Recv(&displacements[0], info[3], MPI_DOUBLE, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
				tag++;
				collections[info[0]].push_back(transform_back(obs_stats, info[1]));
				collections[info[0]].push_back(transform_back(null_stats,info[2]));
				collections[info[0]].push_back(transform_back(displacements,info[3]));
				
			}
	
		}
		typedef map<int, vector<double> >::iterator it_type; 
		typedef map<int, vector<vector<double> >>::iterator it_type_2; 
		
		for (it_type m = observed_statistics.begin(); m!=observed_statistics.end(); m++){
			int * info 	= new int[4];
			int km 				= m->second.size()*observed_null_statistics[m->first].size();

			int od 				= observed_displacements[m->first].size();
			

			//info[0]=PSSM_ID, info[1]= # observed stats, info[2] = # number flattened, info[3]= # of motif hists
			info[0] = m->first, info[1] = m->second.size(), info[2]=km, info[3]=od;

			double * observed_info 	= new double[info[1]];
			double * null_info 	= new double[info[2]];
			double * displacements 	= new double[info[3]];
			
			make_array_1D(m->second, observed_info);
			make_array_1D(observed_displacements[m->first], displacements);
			make_array_2D(observed_null_statistics[m->first], null_info);
			collections[info[0]].push_back(transform_back(observed_info, info[1]));
			collections[info[0]].push_back(transform_back(null_info,info[2]));
			collections[info[0]].push_back(transform_back(displacements,info[3]));

		}
	}else{
		typedef map<int, vector<double> >::iterator it_type; 
		typedef map<int, vector<vector<double> >>::iterator it_type_2; 
		int S;
		S 	= observed_statistics.size();
		int tag = 0;
		MPI_Ssend(&S, 1, MPI_INT, 0,tag, MPI_COMM_WORLD);

		int p =0;
		
		tag++;

		for (it_type m = observed_statistics.begin(); m!=observed_statistics.end(); m++){
			int * info 	= new int[4];
			int km 				= m->second.size()*observed_null_statistics[m->first].size();

			int od 				= observed_displacements[m->first].size();
			

			//info[0]=PSSM_ID, info[1]= # observed stats, info[2] = # number flattened, info[3]= # of motif hists
			info[0] = m->first, info[1] = m->second.size(), info[2]=km, info[3]=od;
			MPI_Ssend(&(info[0]), 4, MPI_INT, 0,tag, MPI_COMM_WORLD);
			tag++;
			double * observed_info 	= new double[info[1]];
			double * null_info 	= new double[info[2]];
			double * displacements 	= new double[info[3]];
			
			make_array_1D(m->second, observed_info);
			make_array_1D(observed_displacements[m->first], displacements);
			make_array_2D(observed_null_statistics[m->first], null_info);

			MPI_Ssend(&(observed_info[0]), m->second.size(), MPI_DOUBLE, 0,tag, MPI_COMM_WORLD);
			tag++;
			MPI_Ssend(&(null_info[0]), km, MPI_DOUBLE, 0,tag, MPI_COMM_WORLD);
			tag++;
			MPI_Ssend(&(displacements[0]), od, MPI_DOUBLE, 0,tag, MPI_COMM_WORLD);
			tag++;

			p++;
			
		}


	}
	return collections;
}

map<int, segment> gather_PSSM_hits_by_bidirectional(int rank, int nprocs, map<string, vector<segment>> intervals){
	typedef map<string, vector<segment>>::iterator it_type;
	
	map<int, segment> m_collapsed;
	vector<segment> collapsed;
	for (it_type i = intervals.begin(); i!=intervals.end();i++){
		for (int j = 0 ; j < i->second.size(); j++){
			collapsed.push_back(i->second[j]);
			m_collapsed[i->second[j].position] 	= i->second[j];
		}
	}
	if (rank==0){
		for (int j = 1 ; j < nprocs; j++ ){
			int tag = 0;
			for (int s = 0; s < collapsed.size(); s++){
				//want to recieve the number of PSSMS in motif positions 
				int P;
				MPI_Recv(&P, 1, MPI_INT, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
				tag++;
				for (int p = 0; p < P; p++){
					int x[3]; //want to recieve the positions to expcet
					MPI_Recv(&(x), 3, MPI_INT, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
					tag++;

					for (int d = 0; d < x[1]; d++){
						double DD;
						MPI_Recv(&DD, 1, MPI_DOUBLE, j, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
						m_collapsed[x[2]].motif_positions[x[0]].push_back(DD);
						tag++;
					}
						
				}
			}
		}




	}else{
		typedef map<int, vector<double>> ::iterator it_type_2;
		int tag 	= 0;
		for (int s = 0; s < collapsed.size();s++){
			int P 	= collapsed[s].motif_positions.size();
			MPI_Ssend(&P, 1, MPI_INT, 0,tag, MPI_COMM_WORLD);
			int p = 0;
			tag++;
			for (it_type_2 m = collapsed[s].motif_positions.begin(); m!=collapsed[s].motif_positions.end(); m++){
				int x[3];
				x[0] 	= m->first, x[1] 	= m->second.size(), x[2]=collapsed[s].position;
				MPI_Ssend(&x[0], 3, MPI_INT, 0,tag, MPI_COMM_WORLD);
				tag++;
				for (int d = 0; d< x[1];d++){
					double DD 	= m->second[d];
					MPI_Ssend(&DD, 1, MPI_DOUBLE, 0,tag, MPI_COMM_WORLD);
					tag++;

				}
				p++;
			}
		}
	}
	return m_collapsed;
}

map<int, vector<double> > send_collect_observed_statistics(int rank, int nprocs, vector<PSSM *> PSSMS, 
	map<int, vector<double>>  observed_displacements){
	map<int, vector<double> > G ;
	if (rank==0){
		for (int j = 1 ; j < nprocs ; j ++ ){
			int S; 
			MPI_Recv(&S, 1, MPI_INT, j, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
			for (int p = 0 ; p < S; p++){
				//(1) want to receive size of data displacement array
				//(1) want to receive PSSMs[p]->ID as well
				//(2) then want to receive the array
				int ID, S2;
				MPI_Recv(&S2, 1, MPI_INT, j, p+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&ID, 1, MPI_INT, j, p+1+S,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				double * displacement 	= new double[S2];
				MPI_Recv(&displacement[0], S2, MPI_DOUBLE, j, p+1+S*2, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
				G[ID] 		= transform_back(displacement, S2);
			}
		}
		for (int p = 0 ; p < PSSMS.size(); p++){
			G[PSSMS[p]->ID] 	= observed_displacements[PSSMS[p]->ID];
		}
	}else{
		//want to send size of PSSMs
		int S 	= PSSMS.size();
		MPI_Ssend(&S, 1, MPI_INT, 0,0, MPI_COMM_WORLD);
		for (int p = 0; p < PSSMS.size(); p++){
			//(1) want to send the PSSM[p]->ID
			//(2) want to send the displacement data observed
			int S2 	= observed_displacements[PSSMS[p]->ID].size();
			MPI_Ssend(&S2, 1, MPI_INT, 0,p+1, MPI_COMM_WORLD);
			int S3 	= PSSMS[p]->ID;
			MPI_Ssend(&S3, 1, MPI_INT, 0,p+1+S, MPI_COMM_WORLD);
			
			double * displacement 	= new double[S2];
			make_array_1D(observed_displacements[PSSMS[p]->ID], displacement   );
			MPI_Ssend(&(displacement[0]), S2, MPI_DOUBLE, 0,p+1+S*2, MPI_COMM_WORLD);


		}
	}
	return G;
}	
























