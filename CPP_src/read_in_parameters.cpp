#include "read_in_parameters.h"
#include <iostream>
using namespace std;
params::params(){
	p["-o"] 	= "";
	p["-bed"] 	= "";
	p["-fasta"] = "";
	p["-DB"] 	= "";
	p["-null"] 	= "1";
	p["-pv"] 	= "0.0000001";
	p["-ID"] 	= "gTFI";
	p["-sim_N"] = "10";
	p["-br"] 	= "10";
	p["-bsn"] 	= "100";

	module 		= "";
	EXIT 		= false;

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
	argv++;
	if (*argv){
		P->module  = string(*argv);
		if (P->module!= "DB" and P->module!="EVAL"){
			printf("do not under module, either DB or EVAL\n");
			P->EXIT 	= true;
		}
	}else{
		printf("No module specified\n");
		P->EXIT 		= true;
	}
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




















