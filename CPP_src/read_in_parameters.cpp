#include "read_in_parameters.h"
#include <iostream>
using namespace std;
params::params(){
	p["-bed"] 			= "";
	p["-fasta"] 		= "";
	p["-DB"] 			= "";
	p["-o"] 			= "";
	p["-ID"] 			= "1";
	p["-log_out"] 		= "";
	p["-br"] 			= "100";
	p["-pv"] 			= "0.0000001";
	p["-simN"] 			= "10";
}
bool params::check(){
	typedef map<string, string >::iterator it_type; 
	bool EXIT = false;
	for (it_type i = p.begin(); i!= p.end(); i++){
		if (i->second.empty()){
			printf("%s was not specified\n",i->first.c_str() );
			EXIT=true;
		}
	}
	return EXIT;
}
void params::help(){

}

void params::display(){
	printf("------------------------------------------------\n" );
	typedef map<string, string >::iterator it_type; 
	for (it_type i = p.begin(); i!= p.end(); i++){
		int dir 	= 20-i->first.size();
		string fillers 	= "";
		for (int i =0 ; i < dir; i++){
			fillers+=" ";
		}
		printf("%s  %s: %s\n", i->first.c_str() , fillers.c_str(), i->second.c_str() );
	}	
	printf("------------------------------------------------\n" );
}

void fill_in_options(char* argv[],params * P, int rank){
	string F 		= "";
	char * COM 		= "-";
	string current 	= "";
	while (*argv){
		F 	= string(*argv); 
		if(!current.empty()){
			P->p[current] 	= F;
			current 		= "";
		}
		if ((*argv)[0] == COM[0]){
			current 	= F;
			if ( P->p.find(F) ==P->p.end() ){
				if (rank == 0){
					printf("Unknown user option: %s\n", F.c_str() );
				}
				current 	= "";
			}
			else if(P->p.find(F) !=P->p.end() ){
				current 	= F;
			}
		}
		if (F.substr(0,2) == "-h" or F.substr(0,7)=="--help" ){
			P->help();
		}
		
		argv++;
	}

	
	if (!current.empty()){
		P->p[current] 	= F;
	}

}




















