#ifndef read_in_parameters_H
#define read_in_parameters_H
#include <map>
using namespace std;

class params{
public:
	map<string, string> p;
	string module;
	bool EXIT;
	params();
	bool check();
	void help();
	void display();
	string get_header();
};
void fill_in_options(int, char** ,params *  , int);


#endif
