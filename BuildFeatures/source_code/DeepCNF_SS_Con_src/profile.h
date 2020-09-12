#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;


//---- definitions ----//
//secondary structure
#ifndef HELIX
#define HELIX	0
#endif
#ifndef SHEET
#define SHEET	1
#endif
#ifndef LOOP
#define LOOP	2
#endif
//solvent accessibility
#ifndef BURIED
#define BURIED		0
#endif
#ifndef INTERMEDIATE
#define INTERMEDIATE	1
#endif
#ifndef EXPOSED
#define EXPOSED		2
#endif

//---- HHpred related ----//
const int M_M = 0;
const int M_I = 1;
const int M_D = 2;
const int I_M = 3;
const int I_I = 4;
const int D_M = 5;
const int D_D = 6;
const int _NEFF = 7;
const int I_NEFF = 8;
const int D_NEFF = 9;



//============== New+Delete+Equal Array ==============//
//new
template <class A>
inline void NewArray2D_(A *** warray,int Narray1,int Narray2)
{
	*warray=new A*[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A[Narray2];
};
template <class A>
inline void NewArray3D_(A **** warray,int Narray1,int Narray2,int Narray3)
{
	*warray=new A**[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A*[Narray2];
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++)  *(*(*warray+i)+j)=new A[Narray3];                  
};
template <class A>
inline void NewArray4D_(A ***** warray,int Narray1,int Narray2,int Narray3,int Narray4)
{
	*warray=new A***[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A**[Narray2];
	for(int i=0;i<Narray1;i++)for(int j=0;j<Narray2;j++) *(*(*warray+i)+j)=new A*[Narray3];
	for(int i=0;i<Narray1;i++)for(int j=0;j<Narray2;j++)for(int k=0;k<Narray3;k++) *(*(*(*warray+i)+j)+k)=new A[Narray4];
};
//delete
template <class A>
inline void DeleteArray2D_(A *** warray,int Narray)
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray) delete [] (*warray);
	(*warray)=NULL;
};
template <class A>
inline void DeleteArray3D_(A **** warray,int Narray1,int Narray2)
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++)  if(*(*(*warray+i)+j)) delete [] *(*(*warray+i)+j);
	for(int i=0;i<Narray1;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray1) delete [] (*warray);
	(*warray)=NULL;
};
template <class A>
inline void DeleteArray4D_(A ***** warray, int Narray1,int Narray2,int Narray3 )
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++) for(int k=0;k<Narray3;k++) if(*(*(*(*warray+i)+j)+k)) delete []*(*(*(*warray+i)+j)+k);
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++) if(*(*(*warray+i)+j)) delete [] *(*(*warray+i)+j);
	for(int i=0;i<Narray1;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray1) delete [] (*warray);
	(*warray)=NULL;
};
//equal
template <class A> 
void EqualArray_(A *out,A *in,int len)
{
	for(int i=0;i<len;i++) *(out+i)=*(in+i);
};
template <class A> 
void EqualArray2D_(A **out,A **in,int len1,int len2)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++) *(*(out+i)+j)=*(*(in+i)+j);
};
template <class A> 
void EqualArray3D_(A ***out,A ***in,int len1,int len2,int len3)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++)for(int k=0;k<len3;k++) *(*(*(out+i)+j)+k)=*(*(*(in+i)+j)+k);
};
template <class A> 
void EqualArray4D_(A ****out,A ****in,int len1,int len2,int len3,int len4)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++)for(int k=0;k<len3;k++)for(int l=0;l<len4;l++) *(*(*(*(out+i)+j)+k)+l)=*(*(*(*(in+i)+j)+k)+l);
};


//------ required data -------//
//the conversion of AA to subScript
const int AA2SUB[26]={0,20,1,2,3,4,5,6,7,20,8,9,10,11,20,12,13,14,15,16,20,17,18,20,19,20};
//AA1Coding
const int AA1Coding[21]={0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18,20};
//AA3Coding
char* const AA3Coding[26]={"ALA","XXX","CYS","ASP","GLU","PHE","GLY","HIS","ILE","XXX","LYS","LEU","MET","ASN","XXX","PRO","GLN","ARG","SER","THR","XXX","VAL","TRP","XXX","TYR","XXX"};

//------ additional data ------//
const float HMMNull[21]={3706,5728,4211,4064,4839,3729,4763,4308,4069,3323,5509,4640,4464,4937,4285,4423,3815,3783,6325,4665,0};
const float HMMNull_f[21]={3.706,5.728,4.211,4.064,4.839,3.729,4.763,4.308,4.069,3.323,5.509,4.64,4.464,4.937,4.285,4.423,3.815,3.783,6.325,4.665,0};
const float PSM_aver[21]={6.4661,13.0051,11.824,13.9727,25.1246,11.04,11.2611,16.757,14.2443,13.6784,14.1105,10.9099,12.1047,17.5259,19.3949,7.7492,7.6928,27.4168,15.3683,10.9346,0};
//to obtain the alphabetical order from the three letter code order, use AA2SUB[ ThreeLetterOrder[i]-'A' ] where i is the order by the three letter code
const char ThreeLetterOrder[21]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','Z'};


//========= class PROFILE ========//
//[note]: read in PSM,PSP and HMM
class PROFILE
{
public:
	int length;
	short *residue;
	float **PSP;
	float **PSM;
	float **ProfHMM;
	float **ProfHMM_original;
	float **EmissionScore;
	float **EmissionScore_original;
	float **EmissionProb;
	float **EmissionProb_original;
	//misc
	int WantBlast;
	int WantOriginal;
	int failure;
	//init
	void Profile_Create_Matrix(int length);
	void Profile_Delete_Matrix(int length);
	//function related
	void Process_PSM(ifstream &fin,string &filename,int StartPos=0);
	void Process_PSP(ifstream &fin,string &filename,int StartPos=0);
	void Process_HMM(ifstream &fin,string &filename,int StartPos=0);
	//constructor
	PROFILE(void);
	~PROFILE(void);
};

