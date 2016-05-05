#ifndef SIM_H
#define SIM_H
#include <map>
#include <string>
#include <vector>
#include "load.h"
#include "error_stdo_logging.h"
using namespace std;
void run_simulations(map<string, vector<segment>>, vector<PSSM *>, int,vector<vector<double>> ,vector<vector<double>>,
					vector<double> , double , int, int ,Log_File * , int,vector<double ** > );
#endif