#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
using namespace std;


//=================== upper and lower case ====================//
//----------upper_case-----------//
void toUpperCase(char *buffer) 
{  
	for(int i=0;i<(int)strlen(buffer);i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
//----------lower_case-----------//
void toLowerCase(char *buffer)
{  
	for(int i=0;i<(int)strlen(buffer);i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}

//=============== Load EVfold MSA, output A3M =============//
//-> example
/*
>4hkrA/1-214
msqsgedlhSPTYLSWRKLQLSRAKLKASSKTSALLSGFAMVAMVEVQLDHDTNVPPGMLIAFAICTTLLVAVHMLALMI
STCILPNIETVSNLHSISLVHESPHERLHWYIETAWAFSTLLGLILFLLEIAILCWVKFYDLSRRAAWSATVVLIPVMII
FMAFAIHFYRSLVSHKYEVTVSGIRELEMLKEQMEQDhlehhnnirnngegeef
>UniRef90_A0A2A6C3L0/5-164
.........-------FQFDLGKAQLKNSSRISALLAGFAVIALVELIYET--NSPRWLLLMIGVLTVLLVSIHLLALMM
STCIQPYMQTSQE------VDESI-LELKFYIDLSWLFSTCVGLIILLAEIGLIFYYQFRIDGEHAGWMTLTIMIPVCFT
FVLVSCFINRIRADLSIERLETKI-------------.................
>UniRef90_A0A2A6C3L0/230-340
.........-YEMAEQYLYDLKKAQLKASSNTSALLSGFAMIALVELHYDA--DTPHWLLIMLGVVTALLVSVHLLALMI
STCIQPYMQAA------GPTQYSPHIRLKFLIDLSWFFSTCVGLILFLA-------------------------------
-------------------------------------.................
>UniRef90_A0A2A6C3L0/373-550
.........---MADKYTYVLSKAQLKASSNTSALLAGFAMVCLVELQYDE--HTPHGLMIVLGVVTALLVSVHLLALMM
STCILPYMEAS------GATQDSPHIRLKFYIDLSWFFSTCIGLILFLVEIGLIFYVKFRAIDFQAGWITTAILIPVLFV
FVIISCFIHRSRANYSFDRIDSKVSGLKKMLNSSEEN.................
>UniRef90_L9K1R8/39-248
....tsnhhSVQALSWRKLYLSRAKLKASSRTSALLSGFAMVAMVEVQLETQYQYPRPLLIAFSACTTVLVAVHLFALLI
STCILPNVEAVSNIHNLNSISESPHERMHPYIELAWGFSTVLGILLFLAEVVLLCWIKFLPVDARAALVSTIIMVPVGLI
FVVFTIHFYRSLVRHKTERHNREIEELHKLKVQLDGHe................
*/

//---------------- change '.' to '-' -----------------//
void Change_Dot_to_Dash(string &in)
{
	int i;
	int len=(int)in.length();
	for(i=0;i<len;i++)
	{
		if(in[i]=='.')in[i]='-';
	}
}


//---------------- multi fasta input -----------------//
int Multi_FASTA_Input(string &multi_fasta,vector <string> &nam_list,vector <string> &fasta_list)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		exit(-1);
	}
	//load
	int first=1;
	int count=0;
	int number=0;
	string name;
	string seq;
	nam_list.clear();
	fasta_list.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf.length()>=1 && buf[0]=='>')
		{
			//get header
			temp=buf.substr(1,buf.length()-1);
			nam_list.push_back(temp);
			count++;
			if(first!=1)
			{
				//transfer to upper
				toUpperCase(seq);
				Change_Dot_to_Dash(seq);
				//add
				fasta_list.push_back(seq);
				number++;
			}
			first=0;
			seq="";
		}
		else
		{
			if(first!=1)seq+=buf;
		}
	}
	//final
	if(first!=1)
	{
		//transfer to upper
		toUpperCase(seq);
		Change_Dot_to_Dash(seq);
		//add
		fasta_list.push_back(seq);
		number++;
	}
	//check
	if(number!=count)
	{
		fprintf(stderr,"%s -> num %d != count %d \n",multi_fasta.c_str(),number,count);
		exit(-1);
	}
	return count;
}

//--------- main ----------//
int main(int argc,char **argv)
{
	//---- EVfold_MSA_Trans ----//
	{
		if(argc<2)
		{
			fprintf(stderr,"./EVfold_MSA_Trans <evfold_a2m> \n");
			exit(-1);
		}
		//input
		string evfold_a2m=argv[1];

		//read
		int retv;
		vector <string> nam_list;
		vector <string> fasta_list;
		retv=Multi_FASTA_Input(evfold_a2m,nam_list,fasta_list);

		//output
		for(int i=0;i<retv;i++)
		{
			printf(">%s\n",nam_list[i].c_str());
			printf("%s\n",fasta_list[i].c_str());
		}

		//exit
		exit(0);
	}
}
