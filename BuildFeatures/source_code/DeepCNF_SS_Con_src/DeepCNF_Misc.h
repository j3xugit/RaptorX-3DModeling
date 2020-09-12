#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

// ----- definitions ------//
#ifdef _FEAT_FLOAT
	typedef float FeatScore;
#else
	typedef double FeatScore;
#endif


// ---- parse string -----//
extern int Parse_Str(string &in,vector <int> &out, char separator);
extern int Parse_Str_Double(string &in,vector <double> &out, char separator);
extern void Parse_Double(vector <double> &in, string &out,char separator);

//------- parse feature range --------//
extern int str2int_(const char *str,int &num);
extern void Parse_Single_Range(string &in_str,int &start,int &end);
extern int Parse_Feature_Range(string &in_str_,vector <int> &range_out);

//------- load training data ---------//
extern int Parse_FeatNum(string &in_str);
extern void Parse_Selection(string &in_str, string &out_str, vector <int> &range_out);
extern int LoadData(string &input_file, int num_procs,int proc_id, 
	int local_num_ori, vector <int> &range_out,
	vector <vector <string> > & feat_in, vector <vector <int> > & label_in);
extern int LoadData(string &input_file, int num_procs,int proc_id, 
	int local_num_ori, vector <int> &range_out,
	vector <vector <vector <FeatScore> > > & feat_in, vector <vector <int> > & label_in);
