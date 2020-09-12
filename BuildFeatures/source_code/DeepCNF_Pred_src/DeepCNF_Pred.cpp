#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <string.h>
#include <omp.h>
#include "getopt.h"
#include "DeepCNF.h"
#include "DeepCNF_Misc.h"
using namespace std;


//==================== proc data =====================//

//--- process data ----//
void ProcData(
	vector <vector <string> > &feat_in, DeepCNF_Model *pModel, 
	vector <vector <vector <double> > > &out_data)
{
	//read in total data
	out_data.resize(feat_in.size());
	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<(int)feat_in.size();i++)
	{
		//-> get length
		int length_s=(int)feat_in[i].size();
		//-> construct a new sequence
		DeepCNF_Seq *seq = new DeepCNF_Seq(length_s, pModel);
		seq->Read_Init_Features(feat_in[i]);
		//-> calculate prob
		vector < vector <double> > output;
		seq->MAP_Probability(output);
		//-> push_back
		out_data[i]=output;
		delete seq;
		//-> count
		fprintf(stderr,"now -> %d \r",i);
	}
}
void ProcData(
	vector <vector <vector <FeatScore> > > &feat_in, DeepCNF_Model *pModel, 
	vector <vector <vector <double> > > &out_data)
{
	//read in total data
	out_data.resize(feat_in.size());
	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<(int)feat_in.size();i++)
	{
		//-> get length
		int length_s=(int)feat_in[i].size();
		//-> construct a new sequence
		DeepCNF_Seq *seq = new DeepCNF_Seq(length_s, pModel);
		seq->Read_Init_Features(feat_in[i]);
		//-> calculate prob
		vector < vector <double> > output;
		seq->MAP_Probability(output);
		//-> push_back
		out_data[i]=output;
		delete seq;
		//-> count
		fprintf(stderr,"now -> %d \r",i);
	}
}


//------- return label ---------//
int Return_Label(vector <double> &in,double &maxval)
{
	int i;
	int size=(int)in.size();
	int label=0;
	double wsmax=0;
	for(i=0;i<size;i++)
	{
		if(in[i]>wsmax)
		{
			wsmax=in[i];
			label=i;
		}
	}
	maxval=wsmax;
	return label;
}


//------------ usage -------------//
void Usage() {
	cerr << "DeepCNF_Pred v1.4 (openmp version) [2015_09_10] \n\n";
	cerr << "Usage : \n\n";
	cerr << "./DeepCNF_Pred -i input_file -w window_str -d node_str -s state_num -l feat_num \n";
	cerr << "               -m init_model [-W label_weight] [-S feat_select] \n";
	cerr << "Options:\n\n";
	//-> required parameter
	cerr << "-i input_file :      input feature file. \n\n";
	cerr << "-w window_str :      window string for DeepCNF. e.g., '5,5' \n\n";
	cerr << "-d node_str :        node string for DeepCNF. e.g., '40,20' \n\n";
	cerr << "-s state_num :       state number. \n\n";
	cerr << "-l feat_num :        feature number at each position. \n\n";
	cerr << "-m init_model :      file for trained model. \n\n";
	cerr << "-W label_weight :    label weight. e.g., '0.1,0.9'. [default is 1 for each label] \n\n";
	cerr << "-S feat_select :     feature selection. e.g., '1-7,9' [default uses all features] \n\n";
}

//------------ main -------------//
int main(int argc, char** argv)
{
	//-- help --//
	if (argc < 2)
	{
		Usage();
		exit(0);
	}

	//---- init parameter ----//
	string input_file = "";
	string window_str = "";    //-> window string
	string node_str = "";      //-> node string
	int state_num = -1;        //-> state number
	int local_num_ori = -1;    //-> original feature number
	int local_num = -1;        //-> processed feature number
	string feat_range_str="";  //-> feature range string
	string model_file = "";
	string label_weight_str="";

	//command-line arguments process
	extern char* optarg;
	char c = 0;
	while ((c = getopt(argc, argv, "i:w:d:s:l:S:m:W:")) != EOF) {
		switch (c) {
		//-> input file
		case 'i':
			input_file = optarg;
			break;
		//-> window/node string
		case 'w':
			window_str = optarg;
			break;
		case 'd':
			node_str = optarg;
			break;
		//-> state/local number
		case 's':
			state_num = atoi(optarg);
			break;
		case 'l':
			local_num_ori = atoi(optarg);
			break;
		case 'S':
			feat_range_str = optarg;
			break;
		//-> model related
		case 'm':
			model_file = optarg;
			break;
		//-> label weight
		case 'W':
			label_weight_str = optarg;
			break;
		//-> default
		default:
			Usage();
			exit(-1);
		}
	}

	//----- check parameter -----//
	//-> check input file
	if(input_file=="")
	{
		fprintf(stderr,"input_file %s is NULL\n",input_file.c_str());
		exit(-1);
	}
	if(model_file=="")
	{
		fprintf(stderr,"model_file %s is NULL\n",model_file.c_str());
		exit(-1);
	}


	//-> check window/node string
	if(window_str=="" || node_str=="")
	{
		fprintf(stderr,"window_str %s or node_str %s is NULL\n",window_str.c_str(),node_str.c_str());
		exit(-1);
	}

	//-> check dimension
	int ws1_;
	{
		vector <int> wstmp;
		char separator=',';
		int ws1=Parse_Str(window_str,wstmp,separator);
		int ws2=Parse_Str(node_str,wstmp,separator);
		ws1_=ws1;
		//--> check window_str and node_str
		if(ws1!=ws2)
		{
			fprintf(stderr,"window_str %s dimension not equal to node_str %s \n",window_str.c_str(),node_str.c_str());
			exit(-1);
		}
	}

	//-> check state/local number
	if(state_num==-1 || local_num_ori==-1)
	{
		fprintf(stderr,"state_num %d and local_num_ori %d must be assigned\n",state_num,local_num_ori);
		exit(-1);
	}
	local_num=local_num_ori;
	//-> process range
	vector <int> range_out(local_num_ori,1);
	if(feat_range_str!="")
	{
		local_num=Parse_Feature_Range(feat_range_str,range_out);
	}

	//-> check label_weight
	vector <double> label_weight;
	if(label_weight_str!="")
	{
		char separator=',';
		int label_weight_num=Parse_Str_Double(label_weight_str,label_weight,separator);
		if(label_weight_num!=state_num)
		{
			fprintf(stderr,"label_weight_str %s dimension not equal to state_num %d \n",label_weight_str.c_str(),state_num);
			exit(-1);
		}
	}
	else
	{
		label_weight.resize(state_num);
		for(int i=0;i<state_num;i++)label_weight[i]=1;
	}

//fprintf -> start
fprintf(stderr,"model_file %s, window_str %s, node_str %s, state_num %d, feat_num %d, feat_range_str %s, label_weight_str %s \n",
	model_file.c_str(),window_str.c_str(),node_str.c_str(),state_num,local_num,feat_range_str.c_str(),label_weight_str.c_str());
//fprintf -> end

	//===================== initilize weights ================//start
	//----------- init model to get the dimension of weights ------------//
	DeepCNF_Model *m_pModel = new DeepCNF_Model(ws1_,window_str,node_str,state_num,local_num,0);
	m_pModel->cnf_model_load(model_file);
	for(int i=0;i<state_num;i++)m_pModel->label_weight[i]=label_weight[i];
	//load data
	vector <vector <string> > feat_in;
	vector <vector <int> > label_in;
	LoadData(input_file, 1, 0, local_num_ori,range_out, feat_in, label_in);

//fprintf -> start
fprintf(stderr,"load data %s done \n",input_file.c_str());
//fprintf -> end

	//proc data
	vector <vector <vector <double> > > prob_out;
	ProcData(feat_in,m_pModel,prob_out);


	//-------- prediction label ---------//
	for(int k=0;k<(int)prob_out.size();k++)
	{
		//-> output predicted result
		int i,j;
		vector <vector <double> > output=prob_out[k];
		printf("#-> %d \n",(int)output.size());
		for(i=0;i<(int)output.size();i++)
		{
			//output
			int label_real=label_in[k][i];
			double maxval;
			int label=Return_Label(output[i],maxval);
			printf("%2d -> ",label_real);
			for(j=0;j<state_num;j++)printf("%lf ",output[i][j]);
			printf("-> %lf %2d \n",maxval,label);
		}
	}

	//exit
	exit(0);
}

