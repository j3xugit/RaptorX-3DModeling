#pragma once
#include <string>
#include <stdio.h>
#include <string.h>
#include "profile.h"
using namespace std;


class TEMPLATE : public PROFILE
{
public:
	string temp_name;
	float NEFF;
	short *isMissing;
	short *isMultiHIS;
	float **CA;
	float **CB;
	short *ACC;
	short *ACCp;
	short *SS;
	short *CLE;
	float **acc;
	float **SS2;
	float **SS8;
	float **acc_;
	float **SS2_;
	float **SS8_;
	short *contactNum;
	short *contactNum_B;
	string sequence;
	string dssp_sequence;
	int StartPos;
	int tLength;
	//init
	void Template_Create_Matrix(int length);
	void Template_Delete_Matrix(int length);
	//functions
	void ReadFeatures_TPL(string filename);		
	TEMPLATE(string templateName,string root=".",int LoadOrig=0,int LoadBlast=0);
	~TEMPLATE(void);
	float DistOf2AAs(int i, int j, int type);
	//--- epad related ---//__140810__//
	double **dis_matrix;
	void Compute_All_Distance(void);
};
