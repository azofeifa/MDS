#include "read_in_parameters.h"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <cctype>
using namespace std;
params::params(){
	p["-TSS"] 	= "";
	p["-o"] 	= "";
	p["-hits"] 	= "0";
	p["-bed"] 	= "";
	p["-fasta"] 	= "";
	p["-DB"] 	= "";
	p["-pv"] 	= "0.0000001";
	p["-ID"] 	= "gTFI";
	p["-sim_N"] 	= "10";
	p["-br"] 	= "100";
	p["-bsn"] 	= "100";
	p["-log_out"] 	= "";
	p["-t"] 	= "0";
	p["-h"] 	= "150";
	p["-H"]         = "3000";

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
	string header 	= "";
	header+="--------------------------------------------------------------------------------------\n";
	header+="                    Global Transcription Factor Inference (tINF)                      \n";
	printf("%s\n",header.c_str() );
	printf("                   ....description of application modules....             \n");
	printf("                                  (critical)                            \n\n");
	
	printf("DB        : must be provided immediately following the application call \"SE\"\n");
	printf("              this creates a database of random sequence similuations\n");
	printf("              each of length 2KB following the natural ATGC frequencey\n");
	printf("              estimated from the set of provided intervals \n");
	printf("EVAL      : must be provided immediately following the application call \"SE\"\n");
	printf("              this will compute enrichment and motif displace scores\n");
	printf("              relative to the provided intervals and the provided database\n");
	printf("              file generated from the \"DB\" module\n");
	
	printf("\n\n");
	header="";
	header+="                ....description of non-default parameters....          \n";
	header+="                                 (critical)                             \n";
	printf("%s\n",header.c_str() );
	printf("-bed      : /path/to/the/bed/file\n");
	printf("              this file should be bedg formatted\n");
	printf("              chromosome[tab]start[tab]stop[\\n]\n");
	printf("              please adhere strictly to that format\n");
	printf("-fasta    : /path/to/the/genome/file\n");
	printf("              this file should be fasta formatted\n");
	printf("              >[chromosome identifier][\\n]ATGCCC....[\\n]\n");
	printf("              this chromosome identifier should match chromosome\n");
	printf("              identifier in bed file otherwise it will be ignored\n");
	printf("-DB       : /path/to/PSSM/motif/file (in case of DB module)\n");
	printf("              PSSM motif file (in case of DB module) should be\n");
	printf("              formatted as >ID,N[\\n] where ID can be anything and N indicates\n");
	printf("              the number TF ChIP-sites this motif was found (N may be ignored)\n");
	printf("              below the >, the PSSM should be formated as a probability matrix\n");
	printf("              rows indicate position in matrix and columns ACGT frequence (in that order) \n");
	printf("          : /path/to/\"DB\"/module/output/file (in case of EVAL module)\n");
	printf("-o        : /path/to/output/file/name\n");
	printf("              the \"DB\" module will create a DB used for \"EVAL\" module\n");
	printf("              the \"EVAL\" module will create a set of enrichment and \n");
	printf("              motif displacement statistics relative to the set of intervals\n");
	printf("\n");
	printf("                    ....description of default parameters....          \n");	
	printf("                                (non-critical)                         \n");
	printf("\n");


	printf("-pv       : (to be used for the \"DB\" module) pvalue threshold cut off for\n");
	printf("              motif hit, recomeneded < 10^-6 (default)\n");
	printf("-sim_N    : (to be used for the \"DB\" module) number of random sequences to\n");
	printf("              generate, recomeneded >10^7 (default)\n");
	printf("-bsn      : (to be used for the \"EVAL\" module) number of bootstrapped samples\n");
	printf("              recomeneded >10^4 (default)\n");
	printf("\nQuestions/Bugs? joseph[dot]azofeifa[at]colorado[dot]edu\n" );		
	printf("--------------------------------------------------------------------------------------------\n");

}
const std::string currentDateTime2() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
}

string params::get_header(){
	string line 	= "";

	if (module.size()==6 and module.substr(0,6)=="GENOME"){
		line+="#================================================================\n";
		line+="#Genome File for Motif Scanner\n";
		line+="#Date/Time      "+currentDateTime2()+"\n";
		line+="#-fasta         "+p["-fasta"]+"\n";
		line+="#-bed           "+p["-bed"]+"\n";
		line+="#-DB            "+p["-DB"]+"\n";
		line+="#-o             "+p["-o"]+"\n";
		line+="#-log_out       "+p["-log_out"]+"\n";
		line+="#-pv            "+p["-pv"]+"\n";
		line+="#================================================================\n";
		
	}

	if (module.size()==4 and module.substr(0,4) == "EVAL"){
		line+="#================================================================\n";
		line+="#Statistics File for Motif Displacement and Enrichment\n";
		line+="#Date/Time      "+currentDateTime2()+"\n";
		line+="#-ID            "+p["-ID"]+"\n";
		line+="#-bed           "+p["-bed"]+"\n";
		line+="#-fasta         "+p["-fasta"]+"\n";
		line+="#-DB            "+p["-DB"]+"\n";
		line+="#-TSS           "+p["-TSS"]+"\n";
		line+="#-o             "+p["-o"]+"\n";
		line+="#-log_out       "+p["-log_out"]+"\n";
		line+="#-pv            "+p["-pv"]+"\n";
		line+="#-bsn           "+p["-bsn"]+"\n";
		line+="#-hits          "+p["-hits"]+"\n";
		line+="#-h             "+p["-h"]+"\n";
		line+="#-H             "+p["-H"]+"\n";
		line+="#================================================================\n";
	}else if(module.size()==2 and module.substr(0,2) == "DB"){
		line+="#================================================================\n";
		line+="#Database File for Motif Displacement and Enrichment\n";
		line+="#Date/Time      "+currentDateTime2()+"\n";
		line+="#-ID            "+p["-ID"]+"\n";
		line+="#-bed           "+p["-bed"]+"\n";
		line+="#-fasta         "+p["-fasta"]+"\n";
		line+="#-DB            "+p["-DB"]+"\n";
		line+="#-TSS           "+p["-TSS"]+"\n";
		line+="#-o             "+p["-o"]+"\n";
		line+="#-log_out       "+p["-log_out"]+"\n";
		line+="#-pv            "+p["-pv"]+"\n";
		line+="#-sim_N         "+p["-sim_N"]+"\n";
		line+="#-H             "+p["-H"]+"\n";
		line+="#================================================================\n";		
	}
	return line;

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

string does_file_exist(string FILE, bool & EXIT){
  ifstream infile(FILE);
  if (infile.good()){
    EXIT=false;
    return "yes";
  }
  EXIT= true;
  return "no";
}

string does_dir_exist(string pathname, bool & EXIT){
  struct stat info;
  if( stat( pathname.c_str(), &info ) != 0 ){
    EXIT =true;
    return "no";
  }
  else if( info.st_mode & S_IFDIR ){ 
    EXIT = false;
    return "yes";
  }
  EXIT =true;
  return "no";
}

string is_number(string x, bool & EXIT){
  for (int i = 0 ; i < x.size(); i++){
    if (!isdigit( x[ i ] ) and x[i]!='.' ){
      EXIT=true;
      return "no";
    }
  }
  EXIT=false;
  return "yes";
}


void fill_in_options(int argc, char* argv[], params * P, int rank){
  if (argc < 2){
    printf("No module specified\ntry -h or --help for quick reference software usage\n");
    P->EXIT                 = true;
  }else{
    string param=argv[1];
    P->module=param;
    if (not (P->module== "DB" or P->module=="EVAL" or P->module=="GENOME" or P->module=="--help" or P->module=="-h"  )){
      printf("did not specify module, either DB or EVAL\ntry -h or --help for quick reference software usage\n");
      P->EXIT         = true;
    }
    
    for (int i = 2 ; i < argc; i++){
      param=argv[i];
      if (P->p.find(param)==P->p.end()){
	printf("Unknown user option: %s\ntry -h or --help for quick reference software usage\n", param.c_str() );
	P->EXIT=true;
      }else{
	if (param=="-hits"){
	  P->p[param]="1";
	}
	if (param!="-hits" and i+1 < argc ){
	  P->p[param]=argv[i+1];
	  i+=1;
	}
      }
    }
  }
  if (not P->EXIT){
    bool temp    = false;
    string result="";
    result       = does_file_exist(P->p["-fasta"], temp);
    if (rank==0){
      printf("checking if fasta file (-fasta) exists................................%s\n", result.c_str());
    }
    P->EXIT      = P->EXIT or temp;
    
    result       = does_file_exist(P->p["-DB"], temp);
    if (rank==0){
      printf("checking if database file (-db) exists................................%s\n", result.c_str());
    }
    P->EXIT      = P->EXIT or temp;

    result       = does_file_exist(P->p["-bed"], temp);
    if (rank==0){
      printf("checking if bed file (-bed) exists....................................%s\n", result.c_str());
    }
    P->EXIT      = P->EXIT or temp;

    result       = does_file_exist(P->p["-TSS"], temp);
    if (rank==0){
      printf("checking if promoter/TSS location bed file (-TSS) exists..............%s\n", result.c_str());
    }
    P->EXIT      = P->EXIT or temp;

    result       = does_dir_exist(P->p["-o"], temp);
    if (rank==0){
      printf("checking if out directory (-o) exists.................................%s\n", result.c_str());
    }
    P->EXIT      = P->EXIT or temp;
    if (P->module=="DB"){
      result       = is_number(P->p["-pv"],temp);
      P->EXIT      = P->EXIT or temp;
      if (not temp){
	float x    = stod(P->p["-pv"]);
	if (x <= 0 or x >= 1){
	  result   = "no", P->EXIT=true;
	}
      }
      if (rank==0){
	printf("checking if pvalue threshold (-pv) is a number in [0,1]...............%s\n", result.c_str());
      }
      
      result       = is_number(P->p["-H"],temp);
      if (rank==0){
	printf("checking if local background radius (-H) is a number..................%s\n", result.c_str());
      }
      P->EXIT      = P->EXIT or temp;
      
      result       = is_number(P->p["-sim_N"],temp);
      if (rank==0){
	printf("checking if number of sequence simulations (-sim_N) is a number.......%s\n", result.c_str());
      }
      P->EXIT      = P->EXIT or temp;
    }
    else if (P->module == "EVAL"){
      result       = is_number(P->p["-h"],temp);
      if (rank==0){
	printf("checking if MDS radius (-h) is a number...............................%s\n", result.c_str());
      }
      P->EXIT      = P->EXIT or temp;

      result       = is_number(P->p["-BSN"],temp);
      if (rank==0){
	printf("checking if MD simulation generation (-bsn) is a number...............%s\n", result.c_str());
      }
      P->EXIT      = P->EXIT or temp;
    }

  }

}




















