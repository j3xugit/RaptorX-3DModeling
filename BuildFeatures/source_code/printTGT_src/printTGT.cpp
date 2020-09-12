#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "seq.h"
using namespace std;


//--------- base_name -----------//__110830__//
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


//---------- printTGT output --------------//
//-> example
/*
mafft-1-clean
310
DNLTKVTEFLLMEFSGIWELQVLHAGLFLLIYLAVLVGNLLIIAVITLDQHLHTPMYFFLKNLSVLDLCYISVTVPKSIRNSLTRRSSISYLGCVAQVYFFSAFASAELAFLTVMSYDRYVAICHPLQYRAVMTSGGCYQMAVTTWLSCFSYAAVHTGNMFREHVCRSSVIHQFFRDIPHVLALVSCEVFFVEFLTLALSSCLVLGCFILMMISYFQIFSTVLRIPSGQSRAKAFSTCSPQLIVIMLFLTTGLFAALGPIAKALSIQDLVIALTYTVLPPFLNPIIYSLRNKEIKTAMWRLFVKIYFLQK
0.096426,0.00997931,0.0656436,0.0939615,0.0155987,0.0936845,0.0419547,0.0283569,0.0757793,0.0712961,0.0128253,0.0466415,0.0281514,0.033096,0.073963,0.0524047,0.0669068,0.0626943,0.010692,0.0199152
0.0521099,0.00927134,0.0644151,0.0519574,0.0138357,0.0575319,0.0341234,0.0198963,0.0507563,0.0379009,0.00988588,0.306041,0.0269086,0.0341805,0.0425636,0.0543555,0.0649483,0.0326685,0.00864084,0.0271795
0.0611163,0.011305,0.0384098,0.0635195,0.0269299,0.0410597,0.0265689,0.0499617,0.0454972,0.177965,0.0217659,0.028472,0.0674353,0.0601667,0.0372176,0.0710769,0.0576301,0.0813799,0.00620021,0.0248821
....
0.331566,-0.918843,0.281798,0.652214,-1.16343,0.312955,0.187975,-0.832158,0.346947,-0.487034,-0.77586,0.217757,-0.686651,0.0198013,0.527948, 0.16884,-0.0867024,-0.212521,-0.222327,-0.984989
 -0.5563,-1.02501,0.254542,-0.202526,-1.33646,-0.390495,-0.110094,-1.34336,-0.231269,-1.39863,-1.15141,  2.9318,-0.751791,0.0663187,-0.269236,0.221569,-0.129564,-1.15295,-0.529613,-0.536339
......
       0,       0,       0,       0,       0,       0,       0,       1,
   0.001,   0.004,       0,   0.023,   0.014,   0.028,   0.001,   0.929,
   0.012,   0.022,       0,    0.12,   0.038,   0.133,   0.242,   0.432,
   0.013,   0.023,       0,   0.219,   0.016,   0.118,     0.2,   0.411,
   0.043,   0.018,       0,   0.336,    0.03,   0.037,   0.187,   0.349,
   0.045,   0.003,       0,   0.347,   0.007,   0.012,   0.068,   0.517,
   0.045,   0.004,       0,   0.458,   0.004,   0.011,   0.051,   0.426,
    0.05,   0.004,       0,   0.679,   0.003,   0.005,   0.163,   0.094,
.....
*/

//------- printTGT ----------//
void PrintTGT(string &tgt_file)
{
	//-- load tgt ---//
	string tgt_name,tgt_root;
	getBaseName(tgt_file,tgt_name,'/','.');
	getRootName(tgt_file,tgt_root,'/');
	SEQUENCE *s=new SEQUENCE(tgt_name,tgt_root,1,1);
	//-- print TGT ---//
	//-> header
	cout << tgt_name << endl;
	cout << s->length << endl;
	cout << s->sequence << endl;
	//-> PSSM (prob)
	for(int k=0;k<s->length;k++)
	{
		//HHpred prob
		for(int i=0;i<20;i++)
		{
			double template_amino_Prob = 1.0*s->EmissionProb[k][i];
			if(i<19)cout <<setw(8)<< template_amino_Prob << ",";
			else cout <<setw(8)<< template_amino_Prob << endl;
		}
	}
	//-> PSSM (score)
	for(int k=0;k<s->length;k++)
	{
		//PSM score
		for(int i=0;i<20;i++)
		{
			double template_amino_Score = 1.0*s->EmissionScore[k][i]+HMMNull[i]/1000.0;
			if(i<19)cout <<setw(8)<< template_amino_Score << ",";
			else cout <<setw(8)<< template_amino_Score << endl;
		}
	}
	//-> SS8
	for(int k=0;k<s->length;k++)
	{
		for(int i=0;i<8;i++)
		{
			if(i<7)cout << setw(8) << s->SS8[k][i] << ",";
			else cout << setw(8) << s->SS8[k][i] << "," << endl;
		}
	}
	//delete
	delete s;
}


//---------- main ----------//
int main(int argc,char **argv)
{
	//---- printTGT ----//
	{
		if(argc<2)
		{
			fprintf(stderr,"printTGT <tgt> \n");
			exit(-1);
		}
		//input
		string tgt=argv[1];
		//proc
		PrintTGT(tgt);
		//exit
		exit(-1);
	}
}

