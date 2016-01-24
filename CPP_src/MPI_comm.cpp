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
	vector<vector<vector<double>>> & streamed_PSSMS, vector<int>  & PSSM_IDS, vector<int> & PSSM_N){
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
			MPI_Recv(&S, 1, MPI_INT, j, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
			for (int p =0; p < S; p++){ //S here is the number of PSSMS
				int * info 	= new int[4];//info[0]=PSSM ID, info[1] number of observed statistics, info[2] number of flattened null statitics
				MPI_Recv(&info[0], 4, MPI_INT, j, p, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				double * obs_stats 	= new double[info[1]];
				double * null_stats = new double[info[2]];
				double * displacements = new double[info[3]];
				
				MPI_Recv(&obs_stats[0], info[1], MPI_DOUBLE, j, p+S, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
				MPI_Recv(&null_stats[0], info[2], MPI_DOUBLE, j, p+2*S, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
				MPI_Recv(&displacements[0], info[3], MPI_DOUBLE, j, p+3*S, MPI_COMM_WORLD,MPI_STATUS_IGNORE);//receiving number of PSSMs scanned
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
		MPI_Ssend(&S, 1, MPI_INT, 0,1, MPI_COMM_WORLD);
		int p =0;
		


		for (it_type m = observed_statistics.begin(); m!=observed_statistics.end(); m++){
			int * info 	= new int[4];
			int km 				= m->second.size()*observed_null_statistics[m->first].size();

			int od 				= observed_displacements[m->first].size();
			

			//info[0]=PSSM_ID, info[1]= # observed stats, info[2] = # number flattened, info[3]= # of motif hists
			info[0] = m->first, info[1] = m->second.size(), info[2]=km, info[3]=od;
			MPI_Ssend(&(info[0]), 4, MPI_INT, 0,p, MPI_COMM_WORLD);
			double * observed_info 	= new double[info[1]];
			double * null_info 	= new double[info[2]];
			double * displacements 	= new double[info[3]];
			
			make_array_1D(m->second, observed_info);
			make_array_1D(observed_displacements[m->first], displacements);
			make_array_2D(observed_null_statistics[m->first], null_info);

			MPI_Ssend(&(observed_info[0]), m->second.size(), MPI_DOUBLE, 0,p+S, MPI_COMM_WORLD);
			MPI_Ssend(&(null_info[0]), km, MPI_DOUBLE, 0,p+S*2, MPI_COMM_WORLD);
			MPI_Ssend(&(displacements[0]), od, MPI_DOUBLE, 0,p+S*3, MPI_COMM_WORLD);



			p++;
			
		}


	}

	return collections;

}

