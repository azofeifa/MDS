#ifndef read_in_parameters_H
#define read_in_parameters_H
#include <map>
using namespace std;

class params{
public:
	map<string, string> p;
	params();
	bool check();
	void help();
	void display();
};
void fill_in_options(char** ,params *  , int);


#endif