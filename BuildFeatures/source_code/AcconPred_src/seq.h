#pragma once
#include <string>
#include <stdio.h>
#include <string.h>
#include "profile.h"
using namespace std;


class SEQUENCE : public PROFILE
{
public:
	string seq_name;
	float NEFF;
	short *isMultiHIS;
	short *ACC;
	short *SS;
	float *DISO;
	float **SS2;
	float **SS8;
	float **acc_our_10_42;
	float **SS2_;
	float **SS8_;
	float **acc_our_10_42_;
	float evd[2];
	string sequence;
	string sse_seq;
	string sse_conf;
	string acc_seq;
	string acc_conf;
	//init
	void Sequence_Create_Matrix(int length);
	void Sequence_Delete_Matrix(int length);
	//function related
	void ReadFeatures_TGT(string filename);
	SEQUENCE(string targetname,string tgt_root=".",int LoadOrig=0,int LoadBlast=0);
	~SEQUENCE(void);
	//distance related
	double*** pair_dis;
	double max_energy,min_energy;
	void Read_distance_file(string filename);
	void Read_distance_file_Binary(double *input,int size,int line_);
};
