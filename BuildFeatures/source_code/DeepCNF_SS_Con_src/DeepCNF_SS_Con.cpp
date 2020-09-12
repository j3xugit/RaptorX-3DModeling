#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string>
#include "getopt.h"
#include "seq.h"
#include "DeepCNF.h"
using namespace std;

//--- model data related ----# start
#include "DeepCNF_data.h"
#include "DeepCNF_model_1.h"
#include "DeepCNF_model_2.h"
#include "DeepCNF_model_3.h"
#include "DeepCNF_model_4.h"
#include "DeepCNF_model_5.h"
#include "DeepCNF_model_6.h"
#include "DeepCNF_model_7.h"
#include "DeepCNF_model_8.h"
#include "DeepCNF_model_con.h"
//--- model data related ----# end

//====== determine model ========//
void Determine_DeepCNF_Model(double * &model_weight, int &feat_num,int feat_lab)
{
	switch (feat_lab) 
	{
		case 1:
			feat_num = Feature_Model_1_size;
			model_weight = (double *)Feature_Model_1;
			break;
		case 2:
			feat_num = Feature_Model_2_size;
			model_weight = (double *)Feature_Model_2;
			break;
		case 3:
			feat_num = Feature_Model_3_size;
			model_weight = (double *)Feature_Model_3;
			break;
		case 4:
			feat_num = Feature_Model_4_size;
			model_weight = (double *)Feature_Model_4;
			break;
		case 5:
			feat_num = Feature_Model_5_size;
			model_weight = (double *)Feature_Model_5;
			break;
		case 6:
			feat_num = Feature_Model_6_size;
			model_weight = (double *)Feature_Model_6;
			break;
		case 7:
			feat_num = Feature_Model_7_size;
			model_weight = (double *)Feature_Model_7;
			break;
		case 8:
			feat_num = Feature_Model_8_size;
			model_weight = (double *)Feature_Model_8;
			break;
		case 0:
			feat_num = Feature_Model_Con_size;
			model_weight = (double *)Feature_Model_Con;
		default:
			exit(-1);
	}
};


//=========================== TGT file part =========================//

//======================= I/O related ==========================//
//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//===================== ACC_Predict ==================//
double ACC_ZY_AmiModel[20][58]={
{0.42,0.8,0.9,0.85,1.05,1.06,1.18,0.73,0.93,0.77,0.87,1.28,0.05,1,0.31,6.11,0.42,0.23,1,0.58,0.514,0.266,0.822,0.709,0.463,0.736,0.723,0.962,0.956,0.513,0.971,0.976,0.861,0.71,0.88,0.949,0.959,0.968,10,6.6,0,0,0,6.6,9,0,0,0,5,6.6,5,0,0,0,0,0,0,0},
{0.61,0.87,0.78,0.65,1.05,1.08,1.06,1.09,1.27,0.86,0.8,2.34,0.29,6.13,-1.01,10.74,0.36,0.25,0.58,1,0.814,0.566,0.605,0.825,0.601,0.82,0.859,0.405,0.389,0.959,0.579,0.465,0.798,0.879,0.809,0.658,0.687,0.417,6.6,10,0,0,0,9,6.6,0,5,0,6.6,9,9,9,0,5,5,6.6,5,0},
{2.24,0.78,1.16,0.94,0.63,1.52,1.68,1.58,0.98,1.55,1.71,1.6,0.13,2.95,-0.6,6.52,0.21,0.22,0.514,0.814,1,0.844,0.65,0.93,0.808,0.926,0.883,0.292,0.269,0.821,0.568,0.393,0.822,0.92,0.806,0.641,0.63,0.31,0,0,10,9,0,0,0,0,5,0,0,0,0,0,0,6.6,0,0,0,0},
{2.56,1.07,2.06,1.95,0.47,0.8,0.78,1.65,0.73,2.06,1.59,1.6,0.11,2.78,-0.77,2.95,0.25,0.2,0.266,0.566,0.844,1,0.431,0.766,0.932,0.724,0.689,0.068,0.053,0.571,0.311,0.152,0.604,0.754,0.6,0.404,0.403,0.092,0,0,9,10,0,0,0,0,6.6,0,0,5,0,0,0,6.6,0,0,0,0},
{0.95,0.68,0.6,0.74,1.21,1.06,1.11,0.82,0.8,0.92,1.2,1.77,0.13,2.43,1.54,6.35,0.17,0.41,0.822,0.605,0.65,0.431,1,0.777,0.511,0.829,0.822,0.73,0.711,0.549,0.856,0.802,0.852,0.75,0.825,0.892,0.845,0.737,0,0,0,0,10,0,0,0,9,0,0,0,0,5,0,9,6.6,5,5,0},
{0.6,0.78,1.04,1.38,0.83,1,1.21,0.95,1.27,0.93,0.67,1.56,0.18,3.95,-0.22,5.65,0.36,0.25,0.709,0.825,0.93,0.766,0.777,1,0.821,0.954,0.922,0.54,0.52,0.819,0.752,0.619,0.93,0.951,0.923,0.814,0.809,0.556,6.6,9,0,0,0,10,6.6,0,5,0,9,9,5,6.6,0,0,5,5,5,0},
{0.58,0.93,1.94,1.63,0.4,0.93,0.8,0.83,1.12,0.93,0.74,1.56,0.15,3.78,-0.64,3.09,0.42,0.21,0.463,0.601,0.808,0.932,0.511,0.821,1,0.739,0.724,0.314,0.303,0.586,0.49,0.369,0.727,0.792,0.719,0.576,0.594,0.335,9,6.6,0,0,0,6.6,10,0,0,0,5,9,0,0,0,0,0,0,0,0},
{1.41,1.66,1.76,1.31,0.77,0.72,2.36,1.77,1.23,0.94,1.72,0,0,0,0,6.07,0.13,0.15,0.736,0.82,0.926,0.724,0.829,0.954,0.739,1,0.935,0.565,0.546,0.789,0.793,0.665,0.929,0.939,0.91,0.845,0.825,0.576,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0},
{1.18,0.88,1.05,1.17,0.85,1.28,1.32,1.06,1.1,1.07,1.08,2.99,0.23,4.66,0.13,7.69,0.27,0.3,0.723,0.859,0.883,0.689,0.822,0.922,0.724,0.935,1,0.552,0.53,0.832,0.773,0.634,0.909,0.96,0.907,0.823,0.817,0.569,0,5,5,6.6,9,5,0,0,10,0,0,0,0,5,0,9,5,5,6.6,0},
{0.29,0.89,0.53,0.71,1.78,0.73,0.48,0.5,0.85,0.62,0.6,4.19,0.19,4,1.8,6.04,0.3,0.45,0.962,0.405,0.292,0.068,0.73,0.54,0.314,0.565,0.552,1,0.998,0.336,0.93,0.983,0.728,0.526,0.76,0.886,0.899,0.997,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,9},
{0.36,0.88,0.51,0.68,1.45,1.17,1.05,0.51,0.79,0.68,0.78,2.59,0.19,4,1.7,6.04,0.39,0.31,0.956,0.389,0.269,0.053,0.711,0.52,0.303,0.546,0.53,0.998,1,0.313,0.922,0.982,0.716,0.503,0.742,0.875,0.893,0.997,5,6.6,0,0,0,9,5,0,0,0,10,9,6.6,5,0,0,0,5,0,0},
{0.56,0.93,0.95,0.73,0.69,1.23,1.14,1.15,1.52,0.86,0.88,1.89,0.22,4.77,-0.99,9.99,0.32,0.27,0.513,0.959,0.821,0.571,0.549,0.819,0.586,0.789,0.832,0.336,0.313,1,0.523,0.388,0.727,0.873,0.785,0.584,0.614,0.347,6.6,9,0,5,0,9,9,0,0,0,9,10,5,5,0,5,5,0,0,0},
{0.42,0.82,0.52,0.72,1.35,1.11,1.07,0.65,0.95,0.67,0.77,2.35,0.22,4.43,1.23,5.71,0.38,0.32,0.971,0.579,0.568,0.311,0.856,0.752,0.49,0.793,0.773,0.93,0.922,0.523,1,0.97,0.891,0.733,0.877,0.981,0.974,0.932,5,9,0,0,0,5,0,0,0,0,6.6,5,10,6.6,0,0,0,6.6,0,0},
{0.52,1.07,0.67,0.91,1.52,1.07,0.9,0.66,0.9,0.74,0.72,2.94,0.29,5.89,1.79,5.67,0.3,0.38,0.976,0.465,0.393,0.152,0.802,0.619,0.369,0.665,0.634,0.983,0.982,0.388,0.97,1,0.8,0.6,0.799,0.937,0.938,0.984,0,9,0,0,5,6.6,0,0,5,0,5,5,6.6,10,0,5,5,9,9,0},
{1,1,1,1,1,1,1,1.84,0.63,2.15,1.05,2.67,0,2.72,0.72,6.8,0.13,0.34,0.861,0.798,0.822,0.604,0.852,0.93,0.727,0.929,0.909,0.728,0.716,0.727,0.891,0.8,1,0.898,0.927,0.941,0.94,0.74,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0},
{2.56,1.13,1.49,1.1,0.65,1.43,1.19,0.98,1.03,1.23,1.28,1.31,0.06,1.6,-0.04,5.7,0.2,0.28,0.71,0.879,0.92,0.754,0.75,0.951,0.792,0.939,0.96,0.526,0.503,0.873,0.733,0.6,0.898,1,0.923,0.792,0.789,0.542,0,5,6.6,6.6,9,0,0,0,9,0,0,5,0,5,0,10,6.6,5,6.6,0},
{2.18,1.11,1.14,1.45,0.93,1.2,0.7,1.02,1.14,0.91,1.29,3.03,0.11,2.6,0.26,5.6,0.21,0.36,0.88,0.809,0.806,0.6,0.825,0.923,0.719,0.91,0.907,0.76,0.742,0.785,0.877,0.799,0.927,0.923,1,0.897,0.908,0.773,0,5,0,0,6.6,5,0,0,5,0,0,5,0,5,0,6.6,10,5,9,0},
{0.42,1.34,0.9,0.93,1.38,0.73,0.55,0.87,0.9,0.78,0.7,3.21,0.41,8.08,2.25,5.94,0.32,0.42,0.949,0.658,0.641,0.404,0.892,0.814,0.576,0.845,0.823,0.886,0.875,0.584,0.981,0.937,0.941,0.792,0.897,1,0.983,0.888,0,6.6,0,0,5,5,0,0,5,0,5,0,6.6,9,0,5,5,10,5,0},
{0.63,1.09,0.82,0.89,1.3,1.2,0.98,0.73,1,0.73,0.77,2.94,0.3,6.47,0.96,5.66,0.25,0.41,0.959,0.687,0.63,0.403,0.845,0.809,0.594,0.825,0.817,0.899,0.893,0.614,0.974,0.938,0.94,0.789,0.908,0.983,1,0.904,0,5,0,0,5,5,0,0,6.6,0,0,0,0,9,0,6.6,9,5,10,0},
{0.26,1.06,0.66,0.99,1.7,0.68,0.46,0.61,0.87,0.61,0.77,3.67,0.14,3,1.22,6.02,0.27,0.49,0.968,0.417,0.31,0.092,0.737,0.556,0.335,0.576,0.569,0.997,0.997,0.347,0.932,0.984,0.74,0.542,0.773,0.888,0.904,1,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,10}};

//-------- process_part --------//
int AA4SUB[26] = { 1, -1,  5,  4,  7, 14,  8, 9,  10,-1, 12,11,  13, 3, -1, 15,  6,  2,  16,  17,-1, 20, 18, -1, 19, -1};


//-------- SS8<->Int --------//
int SS8_To_Int(char c)
{
	switch(c)
	{
		case 'H': return 0;
		case 'G': return 1;
		case 'I': return 2;
		case 'E': return 3;
		case 'B': return 4;
		case 'T': return 5;
		case 'S': return 6;
		case 'L': return 7;
		default: return 7;
	}
}
char Int_To_SS8(int c)
{
	switch(c)
	{
		case 0: return 'H';
		case 1: return 'G';
		case 2: return 'I';
		case 3: return 'E';
		case 4: return 'B';
		case 5: return 'T';
		case 6: return 'S';
		case 7: return 'L';
		default: return 'L';
	}
}

//------ load LABEL_file for predicted contact number ------//
int Load_LAB_File(string &cn_file,vector <int> &lab_number)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(cn_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",cn_file.c_str());
		exit(-1);
	}
	//skip
	for(int i=0;i<2;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file %s format bad!\n",cn_file.c_str());
			exit(-1);
		}
	}
	//load
	lab_number.clear();
	int count=0;
	for(int i=0;i<(int)buf.length();i++)
	{
		char c=buf[i];
		int lab=SS8_To_Int(c);
		lab_number.push_back(lab);
	}
	//return
	return (int)buf.length();
}

//-------------- for probability make -------------//
//given 1st layer feature, output 8 label probability
void Prob_Make(DeepCNF_Model *m_pModel, vector <string> &input,vector <string> &out_str)
{
	//-> construct a new sequence
	DeepCNF_Seq *seq = new DeepCNF_Seq(input.size(), m_pModel);
	seq->Read_Init_Features(input);
	//-> predict ss8 label
	vector <vector <double> > output;
	seq->MAP_Probability(output);
	//-> output
	out_str.clear();
	for(int k=0;k<(int)output.size();k++)
	{
		//----- feature -----//
		vector <double> features=output[k];
		//------ output ------//
		int featdim=(int)features.size();
		stringstream oss;
		for(int i=0;i<featdim;i++)
		{
			int wsiii=(int)features[i];
			if(wsiii!=features[i])oss << features[i] << " ";
			else oss << wsiii << " ";
		}
		string wsbuf=oss.str();
		out_str.push_back(wsbuf);
	}
	//-> delete
	delete seq;
}

//-------------- for primary feature make -------------//
//given one TGT file, output it's primary feature 
int Feature_Make(string &tgt_file,vector <string> &output,string &ami_seq)
{
	//-- read tgt ---//
	string tgt_name,tgt_root;
	getBaseName(tgt_file,tgt_name,'/','.');
	getRootName(tgt_file,tgt_root,'/');
	SEQUENCE *s=new SEQUENCE(tgt_name,tgt_root,0,1);

	//==== generate feature =====//
	output.clear();
	for(int k=0;k<s->length;k++)
	{
		//----- feature -----//
		vector <double> features;
		features.clear();
		//-- profile realted --//
		//emission score
		for(int i=0;i<20;i++)
		{
			float template_amino_Score = 1.0*s->EmissionProb[k][i];
			features.push_back(template_amino_Score);
		}
		//PSM score
		for(int i=0;i<20;i++)
		{
			double template_amino_Prob = 1.0*s->PSM[k][i];
			template_amino_Prob=1.0/(1.0+exp(-1.0*template_amino_Prob));
			features.push_back(template_amino_Prob);
		}
		//-- amino acid realted --//
		//zy ami model
		int pos=AA4SUB[s->sequence[k]-'A']-1;
		vector <int> aaind(20,0);
		if(pos<0 || pos>=20)
		{
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		else
		{
			aaind[pos]=1;
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		//-- secondary structure ---//
		for(int i=0;i<3;i++)
		{
			double psipred_reso = 1.0*s->SS2[k][i];
			features.push_back(psipred_reso);
		}

		//------ output ------//
		int featdim=(int)features.size();
		stringstream oss;
		for(int i=0;i<featdim;i++)
		{
			int wsiii=(int)features[i];
			if(wsiii!=features[i])oss << features[i] << " ";
			else oss << wsiii << " ";
		}
		string wsbuf=oss.str();
		output.push_back(wsbuf);
	}

	//delete
	ami_seq=s->sequence;
	int length=s->length;
	delete s;
	//return
	return length;
}

//------- for secondary feature make ------//
//given one TGT file, output it's secondary feature 
void Feature_Make_II(string &tgt_file,vector <string> &feat_out_ii,string &ami_seq)
{
	//--- generate primary feature ----//
	vector <string> feat_out;
	int length=Feature_Make(tgt_file,feat_out,ami_seq);
	//--- generate secondary feature ---//
	int model_layer=5;
	string window_str = "5,5,5,5,5";
	string node_str = "100,100,100,100,100";
	int state_num = 8;
	int local_num = 63;
	DeepCNF_Model *m_pModel = new DeepCNF_Model(model_layer,window_str,node_str,state_num,local_num,0);
	vector <vector <string> > prob_out_total;
	prob_out_total.clear();
	//for each DeepCNF model
	for(int s=1;s<=8;s++)
	{
		//-> load model
		double * model_weight;
		int feat_num;
		Determine_DeepCNF_Model(model_weight, feat_num,s);
		m_pModel->cnf_model_load(model_weight,feat_num);
		//-> calculate prob
		vector <string> prob_out;
		Prob_Make(m_pModel,feat_out,prob_out);
		prob_out_total.push_back(prob_out);
	}
	delete m_pModel;
	//output
	feat_out_ii.clear();
	for(int k=0;k<length;k++)
	{
		string tmp_str="";
		for(int s=0;s<(int)prob_out_total.size();s++)tmp_str=tmp_str+prob_out_total[s][k]+" ";
		feat_out_ii.push_back(tmp_str);
	}
}


//=========================== MTX+HHM part =========================//

//----------- SS definition ----------//start
#define MAXSEQLEN 50000
#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
//-> logistic 'squashing' function (output range +/- 1.0) 
#define logistic(x) ((double)1.0 / ((double)1.0 + (double)exp(-(x))))
//----------- SS definition ----------//end

//----------- NN structure -----------//start
#define IPERGRP (41)
#define WINL (-7)
#define WINR (7)
#define NUM_IN	((WINR-WINL+1)*IPERGRP)	// number of input units //
#define NUM_HID (75)			// number of hidden units //
#define NUM_OUT (3) 			// number of output units //
#define TOTAL		(NUM_IN + NUM_HID + NUM_OUT)
//----------- NN structure -----------//end

//----------- NN structure -----------//start
#define IPERGRP_ (4)
#define NUM_IN_	((WINR-WINL+1)*IPERGRP_+3+1)	// number of input units //
#define NUM_HID_ (55)			// number of hidden units //
#define TOTAL_		(NUM_IN_ + NUM_HID_ + NUM_OUT)
//----------- NN structure -----------//end


//------- select model -------//
void Determine_NN_Model(double * &model_weight, int model_tag)
{
	switch (model_tag) 
	{
		case 1:
			model_weight = (double *)NN_Data1;
			break;
		case 2:
			model_weight = (double *)NN_Data2;
			break;
		case 3:
			model_weight = (double *)NN_Data3;
			break;
		default:
			exit(-1);
	}
}

//------------ global data structure ---------//
int seqlen, nprof;
int fwt_to[TOTAL], lwt_to[TOTAL];
double activation[TOTAL], bias[TOTAL], *weight[TOTAL];
int fwt_to_[TOTAL_], lwt_to_[TOTAL_];
double activation_[TOTAL_], bias_[TOTAL_], *weight_[TOTAL_];
int profile[MAXSEQLEN][20];
double profile_[MAXSEQLEN][3];
char seq[MAXSEQLEN];
enum aacodes
{
ALA, ARG, ASN, ASP, CYS,
GLN, GLU, GLY, HIS, ILE,
LEU, LYS, MET, PHE, PRO,
SER, THR, TRP, TYR, VAL,
UNK
};


//---------- error function --------//
void err(char *s)
{
	fprintf(stderr, "%s\n", s);
}
void fail(char *s)
{
	err(s);
	exit(1);
}

//------- Convert AA letter to numeric code (0-20) ----//
int aanum(int ch)
{
	static const int aacvs[] =
	{
		999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
		20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
	};
	return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

//--------- gate function --------//
void compute_output(void)
{
	int i, j;
	double netinp;
	for (i = NUM_IN; i < TOTAL; i++)
	{
		netinp = bias[i];
		for (j = fwt_to[i]; j < lwt_to[i]; j++)netinp += activation[j] * weight[i][j];
		activation[i] = logistic(netinp);
	}
}
void compute_output_(void)
{
	int i, j;
	double netinp;
	for (i = NUM_IN_; i < TOTAL_; i++)
	{
		netinp = bias_[i];
		for (j = fwt_to_[i]; j < lwt_to_[i]; j++)netinp += activation_[j] * weight_[i][j];
		activation_[i] = logistic(netinp);
	}
}

//------- Initialize network ------//
void init(void)
{
	int i, j;
	//-- allocate memory --//
	for (i = NUM_IN; i < TOTAL; i++)
	{
		weight[i] = new double[TOTAL - NUM_OUT];
	}
	//-- Connect input units to hidden layer --//
	for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
	{
		fwt_to[i] = 0;
		lwt_to[i] = NUM_IN;
	}
	//-- Connect hidden units to output layer --//
	for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
	{
		fwt_to[i] = NUM_IN;
		lwt_to[i] = NUM_IN + NUM_HID;
	}
}
void init_(void)
{
	int i, j;
	//-- allocate memory --//
	for (i = NUM_IN_; i < TOTAL_; i++)
	{
		weight_[i] = new double[TOTAL_ - NUM_OUT];
	}
	//-- Connect input units to hidden layer --//
	for (i = NUM_IN_; i < NUM_IN_ + NUM_HID_; i++)
	{
		fwt_to_[i] = 0;
		lwt_to_[i] = NUM_IN_;
	}
	//-- Connect hidden units to output layer --//
	for (i = NUM_IN_ + NUM_HID_; i < TOTAL_; i++)
	{
		fwt_to_[i] = NUM_IN_;
		lwt_to_[i] = NUM_IN_ + NUM_HID_;
	}
}

//--------- load weight ---------//
void load_wts(double *input)
{
	int i, j;
	int count=0;
	//-- Load input units to hidden layer weights --//
	for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
	{
		for (j = fwt_to[i]; j < lwt_to[i]; j++)
		{
			weight[i][j] = input[count];
			count++;
		}
	}
	//-- Load hidden layer to output units weights --//
	for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
	{
		for (j = fwt_to[i]; j < lwt_to[i]; j++)
		{
			weight[i][j] = input[count];
			count++;
		}
	}
	//-- Load bias weights --//
	for (j = NUM_IN; j < TOTAL; j++)
	{
		bias[j] = input[count];
		count++;
	}
}
void load_wts_(double *input)
{
	int i, j;
	int count=0;
	//-- Load input units to hidden layer weights --//
	for (i = NUM_IN_; i < NUM_IN_ + NUM_HID_; i++)
	{
		for (j = fwt_to_[i]; j < lwt_to_[i]; j++)
		{
			weight_[i][j] = input[count];
			count++;
		}
	}
	//-- Load hidden layer to output units weights --//
	for (i = NUM_IN_ + NUM_HID_; i < TOTAL_; i++)
	{
		for (j = fwt_to_[i]; j < lwt_to_[i]; j++)
		{
			weight_[i][j] = input[count];
			count++;
		}
	}
	//-- Load bias weights --//
	for (j = NUM_IN_; j < TOTAL_; j++)
	{
		bias_[j] = input[count];
		count++;
	}
}

//------ Make 1st level prediction averaged over specified weight sets ---//
void predict(vector <int> &pred_out, vector <vector <double> > &pred_prob)
{
	int aa, i, j, k, n, winpos,ws;
	char fname[80], predsst[MAXSEQLEN];
	double avout[MAXSEQLEN][3], conf, confsum[MAXSEQLEN];

	//-- init --//
	for (winpos = 0; winpos < seqlen; winpos++)
	{
		avout[winpos][0] = avout[winpos][1] = avout[winpos][2] = confsum[winpos] = 0.0F;
	}

	//-- three model average --//
	for (ws=1; ws<=3; ws++)
	{
		//-> load weight
		double *input;
		Determine_NN_Model(input, ws);
		load_wts(input);
		//-> predict
		for (winpos = 0; winpos < seqlen; winpos++)
		{
			//--> init
			for (j = 0; j < NUM_IN; j++)activation[j] = 0.0;
			//--> feature collect
			for (j = WINL; j <= WINR; j++)
			{
				if (j + winpos >= 0 && j + winpos < seqlen)
				{
					for (aa=0; aa<20; aa++)activation[(j - WINL) * IPERGRP + aa] = profile[j+winpos][aa]/1000.0;
					aa = aanum(seq[j+winpos]);
					if (aa < 20)activation[(j - WINL) * IPERGRP + 20 + aa] = 1.0;
					else activation[(j - WINL) * IPERGRP + 40] = 1.0;
				}
				else
				{
					activation[(j - WINL) * IPERGRP + 40] = 1.0;
				}
			}
			//--> sigmoid
			compute_output();
			//--> output
			conf = 1.0 - MIN(MIN(activation[TOTAL - NUM_OUT], activation[TOTAL - NUM_OUT+1]), activation[TOTAL - NUM_OUT+2]);
			avout[winpos][0] += conf * activation[TOTAL - NUM_OUT];
			avout[winpos][1] += conf * activation[TOTAL - NUM_OUT+1];
			avout[winpos][2] += conf * activation[TOTAL - NUM_OUT+2];
			confsum[winpos] += conf;
		}
	}

	//-- final decision --//
	for (winpos = 0; winpos < seqlen; winpos++)
	{
		avout[winpos][0] /= confsum[winpos];
		avout[winpos][1] /= confsum[winpos];
		avout[winpos][2] /= confsum[winpos];
		if (avout[winpos][0] >= MAX(avout[winpos][1], avout[winpos][2]))predsst[winpos] = 'C';
		else if (avout[winpos][2] >= MAX(avout[winpos][0], avout[winpos][1]))predsst[winpos] = 'E';
		else predsst[winpos] = 'H';
	}

	//-- output --//
	pred_out.clear();
	pred_prob.clear();
	for (winpos = 0; winpos < seqlen; winpos++)
	{
		pred_out.push_back(predsst[winpos]);
		vector <double> pred_porb_;
		pred_porb_.push_back(avout[winpos][0]);
		pred_porb_.push_back(avout[winpos][1]);
		pred_porb_.push_back(avout[winpos][2]);
		pred_prob.push_back(pred_porb_);
	}
}

//------ Make 2nd level prediction ---//
void predict_(int niters, double dca, double dcb,
	vector <int> &pred_out, vector <vector <double> > &pred_prob)
{
	int aa, a, b, nb, i, j, k, n, winpos;
	char pred, predsst[MAXSEQLEN], lastpreds[MAXSEQLEN], *che = "CHE";
	double score_c[MAXSEQLEN], score_h[MAXSEQLEN], score_e[MAXSEQLEN], bestsc, score, conf[MAXSEQLEN], predq3, av_c, av_h, av_e;

	//---- iteration -----//
	if (niters < 1)niters = 1;
	do {
		//-> init
		memcpy(lastpreds, predsst, seqlen);
		av_c = av_h = av_e = 0.0;
		for (winpos = 0; winpos < seqlen; winpos++)
		{
			av_c += profile_[winpos][0];
			av_h += profile_[winpos][1];
			av_e += profile_[winpos][2];
		}
		av_c /= seqlen;
		av_h /= seqlen;
		av_e /= seqlen;
		//-> predict
		for (winpos = 0; winpos < seqlen; winpos++)
		{
			//--> init
			for (j = 0; j < NUM_IN_; j++)activation_[j] = 0.0;
			//--> feature collect
			activation_[(WINR - WINL + 1) * IPERGRP_] = av_c;
			activation_[(WINR - WINL + 1) * IPERGRP_ + 1] = av_h;
			activation_[(WINR - WINL + 1) * IPERGRP_ + 2] = av_e;
			activation_[(WINR - WINL + 1) * IPERGRP_ + 3] = log((double)seqlen);
			for (j = WINL; j <= WINR; j++)
			{
				if (j + winpos >= 0 && j + winpos < seqlen)
				{
					for (aa = 0; aa < 3; aa++)activation_[(j - WINL) * IPERGRP_ + aa] = profile_[j + winpos][aa];
				}
				else
				{
					activation_[(j - WINL) * IPERGRP_ + 3] = 1.0;
				}
			}
			//--> sigmoid
			compute_output_();
			//--> output
			if (activation_[TOTAL_ - NUM_OUT] > dca * activation_[TOTAL_ - NUM_OUT + 1] && activation_[TOTAL_ - NUM_OUT] > dcb * activation_[TOTAL_ - NUM_OUT + 2])pred = 'C';
			else if (dca * activation_[TOTAL_ - NUM_OUT + 1] > activation_[TOTAL_ - NUM_OUT] && dca * activation_[TOTAL_ - NUM_OUT + 1] > dcb * activation_[TOTAL_ - NUM_OUT + 2])pred = 'H';
			else pred = 'E';
			predsst[winpos] = pred;
			score_c[winpos] = activation_[TOTAL_ - NUM_OUT];
			score_h[winpos] = activation_[TOTAL_ - NUM_OUT + 1];
			score_e[winpos] = activation_[TOTAL_ - NUM_OUT + 2];
		}
		//-> result
		for (winpos = 0; winpos < seqlen; winpos++)
		{
			profile_[winpos][0] = score_c[winpos];
			profile_[winpos][1] = score_h[winpos];
			profile_[winpos][2] = score_e[winpos];
		}
	} while (memcmp(predsst, lastpreds, seqlen) && --niters);

	//----- Filter remaining singleton helix/strand assignments ----//
	//-> calculate confident score
	for (winpos = 0; winpos < seqlen; winpos++)
	{
		conf[winpos] = (2*MAX(MAX(score_c[winpos], score_h[winpos]), score_e[winpos])-(score_c[winpos]+score_h[winpos]+score_e[winpos])+MIN(MIN(score_c[winpos], score_h[winpos]), score_e[winpos]));
	}
	//-> first round filter
	for (winpos = 0; winpos < seqlen; winpos++)
	{
		if (winpos && winpos < seqlen - 1 && predsst[winpos - 1] == predsst[winpos + 1] && conf[winpos] < 0.5*(conf[winpos-1]+conf[winpos+1]))
			predsst[winpos] = predsst[winpos - 1];
	}
	//-> second round filter
	for (winpos = 0; winpos < seqlen; winpos++)
	{
		if (winpos && winpos < seqlen - 1 && predsst[winpos - 1] == 'C' && predsst[winpos] != predsst[winpos + 1])
			predsst[winpos] = 'C';
		if (winpos && winpos < seqlen - 1 && predsst[winpos + 1] == 'C' && predsst[winpos] != predsst[winpos - 1])
			predsst[winpos] = 'C';
	}

	//------- output --------//
	pred_out.clear();
	pred_prob.clear();
	for (winpos=0; winpos<seqlen; winpos++)
	{
		pred_out.push_back(predsst[winpos]);
		vector <double> pred_porb_;
		pred_porb_.push_back(score_h[winpos]);
		pred_porb_.push_back(score_e[winpos]);
		pred_porb_.push_back(score_c[winpos]);
		pred_prob.push_back(pred_porb_);
	}
}

//----- Read PSI AA frequency data -----//
int getmtx(FILE *lfil)
{
	int aa, i, j, naa;
	char buf[256], *p;

	//-- init check --//
	if (fscanf(lfil, "%d", &naa) != 1)fail("Bad mtx file - no sequence length!");
	if (naa > MAXSEQLEN)fail("Input sequence too long!");
	if (fscanf(lfil, "%s", seq) != 1)fail("Bad mtx file - no sequence!");

	//-- load file --//
	while (!feof(lfil))
	{
		if (!fgets(buf, 65536, lfil))fail("Bad mtx file!");
		if (!strncmp(buf, "-32768 ", 7))
		{
			for (j=0; j<naa; j++)
			{
				if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)fail("Bad mtx format!");
				aa = aanum(seq[j]);
				if (aa < 20)profile[j][aa] += 0000;
				if (!fgets(buf, 65536, lfil))break;
			}
		}
	}
	//return
	return naa;
}

//----- Read PSI AA frequency data -----//
int getss(vector <vector <double> >&input)
{
	int i;
	int size=(int)input.size();
	for(i=0;i<size;i++)
	{
		profile_[i][0] = input[i][0];
		profile_[i][1] = input[i][1];
		profile_[i][2] = input[i][2];
	}
	//return
	return size;
}

//------------ main process -----------//
void ss_first_layer_predict(string &mtx_file,vector <int> &pred_out2,vector <vector <double> > &pred_prob2)
{
	//first round process
	FILE *ifp = fopen(mtx_file.c_str(), "r");
	if (!ifp)exit(-1);
	seqlen = getmtx(ifp);
	fclose(ifp);
	init();
	vector <int> pred_out1;
	vector <vector <double> > pred_prob1;
	predict(pred_out1,pred_prob1);
	//second round process
	init_();
	getss(pred_prob1);
	load_wts_((double *)NN_Data_2nd);
	int niters=1;
	double dca=0.98;
	double dcb=1.09;
	predict_(niters,dca,dcb,pred_out2,pred_prob2);
}



//======================= second part prediction =====================//

//---------- read in MTX file, output sigmoid feature like ICML 2014 ----------//
//file format
/*
111
WRYRVDVTLSGKKVTGHVLVSLFGNKGNSRQYEIFQGTLKPDNTYSNEFDSDVEVGDLEKVKFIWYNNVINLTLPKVGASKITVERNDGSVFNFCSEETVREDVLLTLTAC
2.670000e-03
4.100000e-02
-3.194183e+00
1.400000e-01
2.670000e-03
5.597866e-02
-2.882785e+00
1.400000e-01
3.176060e-03
1.339561e-01
-2.010243e+00
4.012145e-01
-32768  -684  -32768  -706  -844  -714  -228  -747  -428  -660  -708  -597  -582  -758  -808  -634  -687  -690  -670  -664  1126  -100  789  -32768  -32768  -399  -32768  -32768
-32768  -560  -32768  -745  -556  -424  -624  -640  892  -739  -304  -670  -242  -402  -639  -170  649  -156  -545  -705  -686  -100  -376  -32768  -32768  -399  -32768  -32768
......
*/
int MTX_Proc(string &infile,vector < vector <double> > &output)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"infile %s not found!\n",infile.c_str());
		exit(-1);
	}
	//skip
	for(int i=0;i<14;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"infile %s format bad at %s\n",infile.c_str(),buf.c_str());
			exit(-1);
		}
	}
	//load
	output.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;  // we skip '#' comment
		//process
		istringstream www(buf);
		vector <double> val_raw;
		val_raw.clear();
		for(int j=0;j<23;j++)
		{
			int value_;
			www >> value_;
			double value=1.0/(1.0+exp(-0.01*value_)); //-> sigmoid function
			val_raw.push_back(value);
		}
		//transform
		vector <double> val_rec;
		val_rec.clear();
		{
			val_rec.push_back(val_raw[1]);
			val_rec.push_back(val_raw[16]);
			val_rec.push_back(val_raw[13]);
			val_rec.push_back(val_raw[4]);
			val_rec.push_back(val_raw[3]);
			val_rec.push_back(val_raw[15]);
			val_rec.push_back(val_raw[5]);
			val_rec.push_back(val_raw[7]);
			val_rec.push_back(val_raw[8]);
			val_rec.push_back(val_raw[9]);
			val_rec.push_back(val_raw[11]);
			val_rec.push_back(val_raw[10]);
			val_rec.push_back(val_raw[12]);
			val_rec.push_back(val_raw[6]);
			val_rec.push_back(val_raw[14]);
			val_rec.push_back(val_raw[17]);
			val_rec.push_back(val_raw[18]);
			val_rec.push_back(val_raw[20]);
			val_rec.push_back(val_raw[22]);
			val_rec.push_back(val_raw[19]);
		}
		output.push_back(val_rec);
		count++;
	}
	//return
	return count;
}


//---------- read in HHM file, output EmissionProb ----------//


//------- required matrix ------//
const double gonnet[20][20]= 
{ { 1.7378,  0.870964,0.933254,0.933254, 1.12202,  0.954993, 1,        1.12202,  0.831764, 0.831764,  0.758578, 0.912011, 0.851138, 0.588844, 1.07152,  1.28825,  1.14815,   0.436516,  0.60256,  1.02329},
  { 0.870964,2.95121, 1.07152, 0.933254, 0.60256,  1.41254,  1.09648,  0.794328, 1.14815,  0.57544,   0.60256,  1.86209,  0.676083, 0.47863,  0.812831, 0.954993, 0.954993,  0.691831,  0.660693, 0.630957},
  { 0.933254,1.07152, 2.39883, 1.65959,  0.660693, 1.1749,   1.23027,  1.09648,  1.31826,  0.524807,  0.501187, 1.20226,  0.60256,  0.489779, 0.812831, 1.23027,  1.12202,   0.436516,  0.724436, 0.60256},
  { 0.933254,0.933254,1.65959, 2.95121,  0.47863,  1.23027,  1.86209,  1.02329,  1.09648,  0.416869,  0.398107, 1.12202,  0.501187, 0.354813, 0.851138, 1.12202,  1,         0.301995,  0.524807, 0.512861},
  { 1.12202, 0.60256, 0.660693,0.47863, 14.1254,   0.57544,  0.501187, 0.630957, 0.74131,  0.776247,  0.707946, 0.524807, 0.812831, 0.831764, 0.489779, 1.02329,  0.891251,  0.794328,  0.891251, 1},
  { 0.954993,1.41254, 1.1749,  1.23027,  0.57544,  1.86209,  1.47911,  0.794328, 1.31826,  0.645654,  0.691831, 1.41254,  0.794328, 0.549541, 0.954993, 1.04713,  1,         0.537032,  0.676083, 0.707946},
  { 1,       1.09648, 1.23027, 1.86209,  0.501187, 1.47911,  2.29087,  0.831764, 1.09648,  0.537032,  0.524807, 1.31826,  0.630957, 0.40738,  0.891251, 1.04713,  0.977237,  0.371535,  0.537032, 0.645654},
  { 1.12202, 0.794328,1.09648, 1.02329,  0.630957, 0.794328, 0.831764, 4.57088,  0.724436, 0.354813,  0.363078, 0.776247, 0.446684, 0.301995, 0.691831, 1.09648,  0.776247,  0.398107,  0.398107, 0.467735},
  { 0.831764,1.14815, 1.31826, 1.09648,  0.74131,  1.31826,  1.09648,  0.724436, 3.98107,  0.60256,   0.645654, 1.14815,  0.74131,  0.977237, 0.776247, 0.954993, 0.933254,  0.831764,  1.65959,  0.630957},
  { 0.831764,0.57544, 0.524807,0.416869, 0.776247, 0.645654, 0.537032, 0.354813, 0.60256,  2.51189,   1.90546,  0.616595, 1.77828,  1.25893,  0.549541, 0.660693, 0.870964,  0.660693,  0.851138, 2.04174},
  { 0.758578,0.60256, 0.501187,0.398107, 0.707946, 0.691831, 0.524807, 0.363078, 0.645654, 1.90546,   2.51189,  0.616595, 1.90546,  1.58489,  0.588844, 0.616595, 0.74131,   0.851138,  1,        1.51356},
  { 0.912011,1.86209, 1.20226, 1.12202,  0.524807, 1.41254,  1.31826,  0.776247, 1.14815,  0.616595,  0.616595, 2.0893,   0.724436, 0.467735, 0.870964, 1.02329,  1.02329,   0.446684,  0.616595, 0.676083},
  { 0.851138,0.676083,0.60256, 0.501187, 0.812831, 0.794328, 0.630957, 0.446684, 0.74131,  1.77828,   1.90546,  0.724436, 2.69153,  1.44544,  0.57544,  0.724436, 0.870964,  0.794328,  0.954993, 1.44544},
  { 0.588844,0.47863, 0.489779,0.354813, 0.831764, 0.549541, 0.40738,  0.301995, 0.977237, 1.25893,   1.58489,  0.467735, 1.44544,  5.01187,  0.416869, 0.524807, 0.60256,   2.29087,   3.23594,  1.02329},
  { 1.07152, 0.812831,0.812831,0.851138, 0.489779, 0.954993, 0.891251, 0.691831, 0.776247, 0.549541,  0.588844, 0.870964, 0.57544,  0.416869, 5.7544,   1.09648,  1.02329,   0.316228,  0.489779, 0.660693},
  { 1.28825, 0.954993,1.23027, 1.12202,  1.02329,  1.04713,  1.04713,  1.09648,  0.954993, 0.660693,  0.616595, 1.02329,  0.724436, 0.524807, 1.09648,  1.65959,  1.41254,   0.467735,  0.645654, 0.794328},
  { 1.14815, 0.954993,1.12202, 1,        0.891251, 1,        0.977237, 0.776247, 0.933254, 0.870964,  0.74131,  1.02329,  0.870964, 0.60256,  1.02329,  1.41254,  1.77828,   0.446684,  0.645654, 1},
  { 0.436516,0.691831,0.436516,0.301995, 0.794328, 0.537032, 0.371535, 0.398107, 0.831764, 0.660693,  0.851138, 0.446684, 0.794328, 2.29087,  0.316228, 0.467735, 0.446684, 26.3027,    2.5704,   0.549541},
  { 0.60256, 0.660693,0.724436,0.524807, 0.891251, 0.676083, 0.537032, 0.398107, 1.65959,  0.851138,  1,        0.616595, 0.954993, 3.23594,  0.489779, 0.645654, 0.645654,  2.5704,    6.0256,   0.776247},
  { 1.02329, 0.630957,0.60256, 0.512861, 1,        0.707946, 0.645654, 0.467735, 0.630957, 2.04174,   1.51356,  0.676083, 1.44544,  1.02329,  0.660693, 0.794328, 1,         0.549541,  0.776247, 2.18776} };

/*
//------ required data -------//
//the conversion of AA to subScript
const int AA2SUB[26]={0,20,1,2,3,4,5,6,7,20,8,9,10,11,20,12,13,14,15,16,20,17,18,20,19,20};
//AA1Coding
const int AA1Coding[21]={0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18,20};
//AA3Coding
char* const AA3Coding[26]={"ALA","XXX","CYS","ASP","GLU","PHE","GLY","HIS","ILE","XXX","LYS","LEU","MET","ASN","XXX","PRO","GLN","ARG","SER","THR","XXX","VAL","TRP","XXX","TYR","XXX"};
//------ additional data ------//
const double HMMNull[21]={3706,5728,4211,4064,4839,3729,4763,4308,4069,3323,5509,4640,4464,4937,4285,4423,3815,3783,6325,4665,0};
const double HMMNull_f[21]={3.706,5.728,4.211,4.064,4.839,3.729,4.763,4.308,4.069,3.323,5.509,4.64,4.464,4.937,4.285,4.423,3.815,3.783,6.325,4.665,0};
const double PSM_aver[21]={6.4661,13.0051,11.824,13.9727,25.1246,11.04,11.2611,16.757,14.2443,13.6784,14.1105,10.9099,12.1047,17.5259,19.3949,7.7492,7.6928,27.4168,15.3683,10.9346,0};
//to obtain the alphabetical order from the three letter code order, use AA2SUB[ ThreeLetterOrder[i]-'A' ] where i is the order by the three letter code
const char ThreeLetterOrder[21]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','Z'};
//-------- process_part --------//
int AA4SUB[26] = { 1, -1,  5,  4,  7, 14,  8, 9,  10,-1, 12,11,  13, 3, -1, 15,  6,  2,  16,  17,-1, 20, 18, -1, 19, -1};
*/

//file format
/*
>gi|289643204|ref|ZP_06475331.1|(75-163:164) blue (type 1) copper domain protein [Frankia symbiont of Datisca glomerata]   gi|289506975|gb|EFD27947.1| blue (type 1) copper domain protein [Frankia symbiont of Datisca glomerata]  E=1e-05 s/c=0.60 id=25% cov=52%
-----------------FSPSTLTARPGKIrvnLTVPAGSAPHNFV-IPEIPEARTSVASAGTSQSVTFTIDKPGSYSFLCTIH----------------
--------------------------
>gi|227498043|ref|ZP_03928216.1|(40-139:453) possible nitrite reductase (NO-forming) [Actinomyces urogenitalis DSM 15434]   gi|226832551|gb|EEH64934.1| possible nitrite reductase (NO-forming) [Actinomyces urogenitalis DSM 15434]  E=0.001 s/c=0.47 id=22% cov=69%
---------------MRFVPDTVQVQAGDRLVITLDNTSDQVH--DLVLDTGATTGRVSAGGSARLEVGlVTGPLQGWCSIagHRAQGMVMQVVVGAAAS
GTDA----------------------
#
NULL   3706     5728    4211    4064    4839    3729    4763    4308    4069    3323    5509    4640    4464    4937    4285    4423    3815    3783    6325    4665
HMM    A        C       D       E       F       G       H       I       K       L       M       N       P       Q       R       S       T       V       W       Y
       M->M     M->I    M->D    I->M    I->I    D->M    D->D    Neff    Neff_I  Neff_D
       0        *       *       0       *       0       *       *       *       *
E 1    *        *       1963    428     *       *       *       *       *       *       *       *       *       *       *       *       *       *       *       *       1
       0        *       *       *       *       *       *       3507    0       0


*/
void HHM_Proc(string &infile,int length,string &sequence,
	vector < vector <double> > &EmissionScore,vector < vector <double> > &EmissionProb)
{
	ifstream fin;
	string wbuf,buf,temp;
	string filename=infile;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"infile %s not found!\n",infile.c_str());
		exit(-1);
	}
	//skip
	for(;;)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"infile %s format bad at %s\n",infile.c_str(),buf.c_str());
			exit(-1);
		}
		istringstream www(buf);
		www>>temp;
		if(temp=="NULL")break;
	}
	//load
	EmissionScore.resize(length);
	for(int i=0;i<length;i++)EmissionScore[i].resize(21);
	EmissionProb.resize(length);
	for(int i=0;i<length;i++)EmissionProb[i].resize(21);
	sequence.clear();
	int StartPos=0;
	//skip first
	for(int i=0;i<3;i++)
	{
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM skip first] \n",filename.c_str());
			exit(-1);
		}
	}
	//process string
	for(int i=0;i<length;i++)
	{
		//get string 1
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d] \n",filename.c_str(),i);
			exit(-1);
		}
		if(wbuf.length()<20)
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d] \n",filename.c_str(),i);
			exit(-1);
		}
		istringstream _em(wbuf);
		int residue;
		for(int j=0;j<2;j++) //skip header
		{
			if(!(_em>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d header %d] \n",filename.c_str(),i,j);
				exit(-1);
			}
			if(j==0)
			{
				char c=temp[0];
				if(c<'A'||c>'Z')residue=20;
				else residue = AA1Coding[AA2SUB[c-'A']];
				sequence.push_back(c);
			}
		}
		for(int j=0;j<20;j++) //get remaining
		{
			if(!(_em>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d col %d] \n",filename.c_str(),i,j);
				exit(-1);
			}
			if(temp[0]=='*')
			{
				EmissionScore[i-StartPos][j]=-99999;
			}
			else
			{
				EmissionScore[i-StartPos][j]=-atof(temp.c_str());
			}
			EmissionProb[i-StartPos][j]=pow(2.0,1.0*EmissionScore[i-StartPos][j]/1000.);
		}
		//check
		float wssum=0;
		for(int j=0;j<20;j++)wssum+=EmissionProb[i-StartPos][j];
		if(wssum<0.9) //this is mainly due to NEFF=1
		{
			for(int j=0;j<20;j++)
			{
				EmissionProb[i-StartPos][j]=0; //this bug is discovered by Jinbo
			}
			EmissionProb[i-StartPos][residue]=1;
		}
		if(wssum>1.1) //this is wired
		{
			for(int j=0;j<20;j++)
			{
				EmissionProb[i-StartPos][j]=0; //this bug is discovered by Jinbo
			}
			EmissionProb[i-StartPos][residue]=1;
		}
		//get string 2
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d] \n",filename.c_str(),i);
			exit(-1);
		}
		if(wbuf.length()<7)
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HMM line2 %d] \n",filename.c_str(),i);
			exit(-1);
		}
		stringstream _sin(wbuf);
		//[1] M_M
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_M] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//[2] M_I
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_I] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//[3] M_D
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_D] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//[4] I_M
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans I_M] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//[5] I_I
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans I_I] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//[6] D_M
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans D_M] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//[7] D_D
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans D_D] \n",filename.c_str(),i);
			exit(-1);
		}
		if(temp[0]=='*')
		{
		}
		else
		{
		}
		//final process
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_NEFF] \n",filename.c_str(),i);
			exit(-1);
		}
		double ws_tmp_neff=atof(temp.c_str())/1000.0 - 1;
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans I_NEFF] \n",filename.c_str(),i);
			exit(-1);
		}
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans D_NEFF] \n",filename.c_str(),i);
			exit(-1);
		}
		//get string 3
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line3 %d] \n",filename.c_str(),i);
			exit(-1);
		}
		//final process
		for(int j=0;j<20;j++)
		{
			double g=0;
			int ind = AA1Coding[j];
			for(int k=0;k<20;k++)g+=EmissionProb[i-StartPos][k]*gonnet[AA1Coding[k]][ind]*pow(2.0,-1.0*HMMNull[j]/1000.0);
			EmissionScore[i-StartPos][j] = (ws_tmp_neff*EmissionProb[i-StartPos][j]+g*10)/(ws_tmp_neff+10);
		}
		for(int j=0;j<20;j++)
		{
			EmissionProb[i-StartPos][j] = EmissionScore[i-StartPos][j];
			EmissionScore[i-StartPos][j] = log(EmissionProb[i-StartPos][j])/log(2.);
		}
	}
}



//-------------- for prediction -------------//
//given mtx and hhm file, generate the primary features
int Feature_Make(string &pre_mtx_file,string &cur_mtx_file,string &hhm_file,vector <string> &output,string &ami_seq)
{
	//--- load feature ---//
	//-> load MTX file
	vector < vector <double> > pre_mtx_output;
	int length1=MTX_Proc(pre_mtx_file,pre_mtx_output);
	vector < vector <double> > cur_mtx_output;
	int length2=MTX_Proc(cur_mtx_file,cur_mtx_output);
	if(length1!=length2)
	{
		fprintf(stderr,"pre_mtx_file length %d not equal to cur_mtx_file length %d \n",length1,length2);
		exit(-1);
	}
	int length=length1;
	//-> load HHM file
	vector < vector <double> > EmissionScore;
	vector < vector <double> > EmissionProb;
	string sequence;
	HHM_Proc(hhm_file,length,sequence,EmissionScore,EmissionProb);

	//----- first layer prediction ----//
	vector <int> pred_out2;
	vector <vector <double> > pred_prob2;
	ss_first_layer_predict(pre_mtx_file,pred_out2,pred_prob2);

	//----- second layer prediction ---//
	//==== generate feature =====//
	output.clear();
	ami_seq="";
	for(int k=0;k<length;k++)
	{
		//----- feature -----//
		vector <double> features;
		features.clear();
		//-- profile realted --//
		//emission score
		for(int i=0;i<20;i++)
		{
			double template_amino_Score = 1.0*EmissionProb[k][i];
			features.push_back(template_amino_Score);
		}
		//PSM score
		for(int i=0;i<20;i++)
		{
			double template_amino_Prob = 1.0*cur_mtx_output[k][i];
			features.push_back(template_amino_Prob);
		}
		//-- amino acid sequence --//
		int pos=AA4SUB[sequence[k]-'A']-1;
		vector <int> aaind(20,0);
		if(pos<0 || pos>=20)
		{
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		else
		{
			aaind[pos]=1;
			for(int i=0;i<20;i++)
			{
				features.push_back(aaind[i]);
			}
		}
		//-- first layer prediction ---//
		for(int i=0;i<3;i++)
		{
			double psipred_reso = 1.0*pred_prob2[k][i];
			features.push_back(psipred_reso);
		}

		//------ output ------//
		int featdim=(int)features.size();
		stringstream oss;
		for(int i=0;i<featdim;i++)
		{
			int wsiii=(int)features[i];
			if(wsiii!=features[i])oss << features[i] << " ";
			else oss << wsiii << " ";
		}
		string wsbuf=oss.str();
		output.push_back(wsbuf);
		ami_seq.push_back(sequence[k]);
	}
	//return
	return length;
}

//given mtx and hhm file, generate the secondary features
void Feature_Make_II(string &pre_mtx_file,string &cur_mtx_file,string &hhm_file,vector <string> &feat_out_ii,string &ami_seq)
{
	//--- generate primary feature ----//
	vector <string> feat_out;
	int length=Feature_Make(pre_mtx_file,cur_mtx_file,hhm_file,feat_out,ami_seq);
	//--- generate secondary feature ---//
	int model_layer=5;
	string window_str = "5,5,5,5,5";
	string node_str = "100,100,100,100,100";
	int state_num = 8;
	int local_num = 63;
	DeepCNF_Model *m_pModel = new DeepCNF_Model(model_layer,window_str,node_str,state_num,local_num,0);
	vector <vector <string> > prob_out_total;
	prob_out_total.clear();
	//for each DeepCNF model
	for(int s=1;s<=8;s++)
	{
		//-> load model
		double * model_weight;
		int feat_num;
		Determine_DeepCNF_Model(model_weight, feat_num,s);
		m_pModel->cnf_model_load(model_weight,feat_num);
		//-> calculate prob
		vector <string> prob_out;
		Prob_Make(m_pModel,feat_out,prob_out);
		prob_out_total.push_back(prob_out);
	}
	delete m_pModel;
	//output
	feat_out_ii.clear();
	for(int k=0;k<length;k++)
	{
		string tmp_str="";
		for(int s=0;s<(int)prob_out_total.size();s++)tmp_str=tmp_str+prob_out_total[s][k]+" ";
		feat_out_ii.push_back(tmp_str);
	}
}




//------- return label ---------//
/*
Secondary structure labels, with the sequence of 'L', 'B', 'E', 'G', 'I', 'H', 'S', 'T'
*/
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
int SS8_To_SS3(int in)
{
	switch(in) 
	{
		case 0: return 0; // H->H
		case 1: return 2; // G->C
		case 2: return 2; // I->C
		case 3: return 1; // E->E
		case 4: return 2; // B->C
		case 5: return 2; // T->C
		case 6: return 2; // S->C
		case 7: return 2; // L->C
		default: return 2;
	}
}
void SS8_To_SS3_Prob(vector <double> &in, vector <double> &out)
{
	out.resize(3);
	out[0]=in[0];
	out[1]=in[3];
	out[2]=in[5]+in[6]+in[7]+in[1]+in[2]+in[4];
}

//--- SS to Char ---//
char SS8_To_Char(int in)
{
	switch(in) 
	{
		case 0: return 'H'; // H
		case 1: return 'G'; // G
		case 2: return 'I'; // I
		case 3: return 'E'; // E
		case 4: return 'B'; // B
		case 5: return 'T'; // T
		case 6: return 'S'; // S
		case 7: return 'L'; // L
		default: return 'L';
	}
}
char SS3_To_Char(int in)
{
	switch(in) 
	{
		case 0: return 'H'; // H
		case 1: return 'E'; // E
		case 2: return 'C'; // C
		default: return 'C';
	}
}


//------------ usage -------------//
void Usage() {
	cerr << "DeepCNF_SS v1.02 [Apr-05-2015] \n\n";
	cerr << "./DeepCNF_SS -i pre_mtx_file -I cur_mtx_file -h hhm_file [-s SS_type] \n";
	//cerr << "./DeepCNF_SS -t tgt_file [-s SS_type] \n";
	cerr << "Options:\n\n";
	//-> required parameter
	cerr << "-i pre_mtx_file :    3-iteration mtx file. \n\n";
	cerr << "-I cur_hhm_file :    5-iteration mtx file. \n\n";
	cerr << "-h hhm_file :        input hhm file. \n\n";
	//cerr << "-t tgt_file :        input tgt file. \n\n";
	cerr << "-s SS_type :         0 for SS8 and 1 for SS3. [default=0] \n\n";
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
	string pre_mtx_file = "";
	string cur_mtx_file = "";
	string input_hhm_file = "";
	string input_tgt_file = "";
	int input_type=0;  // [0] for three files; [1] for tgt file
	int SS_type=0;
	int model_layer=5;
	string window_str = "5,5,5,5,5";
	string node_str = "100,100,100,100,100";
	int state_num = 8;
	int local_num = 64; //-> consensus mode, 8*8=64 input features

	//command-line arguments process
	extern char* optarg;
	char c = 0;
	while ((c = getopt(argc, argv, "i:I:h:t:s:")) != EOF) {
		switch (c) {
		//-> input file
		case 'i':
			pre_mtx_file = optarg;
			break;
		case 'I':
			cur_mtx_file = optarg;
			break;
		case 'h':
			input_hhm_file = optarg;
			break;
		//-> input tgt file
		case 't':
			input_tgt_file = optarg;
			break;
		//-> prediction style
		case 's':
			SS_type = atoi(optarg);
			break;
		//-> default
		default:
			Usage();
			exit(-1);
		}
	}

	//----- check parameter -----//
	//-> check input file
	if(pre_mtx_file=="" || cur_mtx_file=="" || input_hhm_file=="")
	{
		if(input_tgt_file=="")
		{
			fprintf(stderr,"input file %s, %s, %s, or input tgt %s is NULL\n",
				pre_mtx_file.c_str(),cur_mtx_file.c_str(),input_hhm_file.c_str(),input_tgt_file.c_str());
			exit(-1);
		}
		pre_mtx_file="";
		cur_mtx_file="";
		input_hhm_file="";
		input_type=1;  //-> tgt input style
	}

	//-> check state/local number
	if(SS_type<0 || SS_type>1)
	{
		fprintf(stderr,"SS_type %d should be 0 or 1 \n",SS_type);
		exit(-1);
	}

	//===================== initilize weights ================//start
	//----------- init model to get the dimension of weights ------------//
	DeepCNF_Model *m_pModel;
	m_pModel = new DeepCNF_Model(model_layer,window_str,node_str,state_num,local_num,0);
	m_pModel->cnf_model_load((double *)Feature_Model_Con,Feature_Model_Con_size);
	//-> load feature
	vector <string> feat_out;
	string ami_seq;
	if(input_type==0)Feature_Make_II(pre_mtx_file,cur_mtx_file,input_hhm_file,feat_out,ami_seq);
	else Feature_Make_II(input_tgt_file,feat_out,ami_seq);

	//-> construct a new sequence
	DeepCNF_Seq *seq = new DeepCNF_Seq(feat_out.size(), m_pModel);
	seq->Read_Init_Features(feat_out);
	//-> predict ss8 label
	vector <vector <double> > output;
	seq->MAP_Probability(output);
	//-> output header
	if(SS_type==0) //--> output SS8
	{
		printf("#DeepConCNF_SS8: eight-class secondary structure prediction results \n");
		printf("#probabilities are in the order of H G I E B T S L(loops), the 8 secondary structure types used in DSSP \n\n");
	}
	else           //--> output SS3
	{
		printf("#DeepConCNF_SS3: three-class secondary structure prediction results \n");
		printf("#probabilities are in the order of H E C, the 3 secondary structure types used in DSSP \n\n");
	}
	//-> output content
//   1 E L   0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.000
//   2 N L   0.000 0.000 0.000 0.379 0.004 0.000 0.000 0.617
	for(int i=0;i<(int)output.size();i++)
	{
		if(SS_type==0) //--> output SS8
		{
			//return maximal value
			double maxval;
			int label=Return_Label(output[i],maxval);
			char lab_c=SS8_To_Char(label);
			//output
			printf("%4d %c %c   ",i+1,ami_seq[i],lab_c);
			for(int j=0;j<state_num;j++)printf("%5.3f ",output[i][j]);
			printf("\n");
		}
		else           //--> output SS3
		{
			int state_num=3;
			//transfer to 3-SSE
			vector <double> output_sse;
			SS8_To_SS3_Prob(output[i],output_sse);
			//return maximal value
			double maxval;
			int label=Return_Label(output_sse,maxval);
			char lab_c=SS3_To_Char(label);
			//output
			printf("%4d %c %c   ",i+1,ami_seq[i],lab_c);
			for(int j=0;j<state_num;j++)printf("%5.3f ",output_sse[j]);
			printf("\n");
		}
	}

	//delete
	delete seq;
	delete m_pModel;

	//exit
	exit(0);
}
