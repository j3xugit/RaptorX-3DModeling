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

//---------- read FASTA file ----------//
void remove_blank(string &in_str,string &out_str)
{
	int i;
	int len=(int)in_str.length();
	out_str="";
	for(i=0;i<len;i++)
	{
		if(in_str[i]!=' ')out_str+=in_str[i];
	}
}
int Read_FASTA_SEQRES(string &fasta_file,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fasta_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",fasta_file.c_str());
		return -1;
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file format bad! %s \n",fasta_file.c_str());
			return -1;
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	remove_blank(temp,seqres);
	return (int)seqres.length();
}

//------------- sequence to SEQRES ------------//
const char* One2Three_III(char c)
{
	//switch
	switch(c)
	{
		case 'A':return "ALA";
		case 'R':return "ARG";
		case 'N':return "ASN";
		case 'D':return "ASP";
		case 'C':return "CYS";
		case 'Q':return "GLN";
		case 'E':return "GLU";
		case 'G':return "GLY";
		case 'H':return "HIS";
		case 'I':return "ILE";
		case 'L':return "LEU";
		case 'K':return "LYS";
		case 'M':return "MET";
		case 'F':return "PHE";
		case 'P':return "PRO";
		case 'S':return "SER";
		case 'T':return "THR";
		case 'W':return "TRP";
		case 'Y':return "TYR";
		case 'V':return "VAL";
		default:return "UNK";
	}
}

//------ output format -----//
/*
GLU ASN ILE GLU VAL HIS MET LEU ASN LYS GLY ALA GLU GLY ALA MET 
VAL PHE GLU PRO ALA TYR ILE LYS ALA ASN PRO GLY ASP THR VAL THR 
PHE ILE PRO VAL ASP LYS GLY HIS ASN VAL GLU SER ILE LYS ASP MET 
ILE PRO GLU GLY ALA GLU LYS PHE LYS SER LYS ILE ASN GLU ASN TYR 
VAL LEU THR VAL THR GLN PRO GLY ALA TYR LEU VAL LYS CYS THR PRO 
HIS TYR ALA MET GLY MET ILE ALA LEU ILE ALA VAL GLY ASP SER PRO 
ALA ASN LEU ASP GLN ILE VAL SER ALA LYS LYS PRO LYS ILE VAL GLN 
PRO ARG LEU GLU LYS VAL ILE ALA SER ALA LYS 
*/

//------------ main -------------//
int main(int argc, char** argv)
{
	//------- SEQ_To_Input -----//
	{
		if(argc<3)
		{
			fprintf(stderr,"SEQ_To_Input <input_fasta> <output_file> \n");
			exit(-1);
		}
		string input_fasta=argv[1];
		string output_file=argv[2];
		//load input seqres
		string seqres;
		int seqlen=Read_FASTA_SEQRES(input_fasta,seqres);
		//output three-letter file
		int i;
		FILE *fp=fopen(output_file.c_str(),"wb");
		for(i=0;i<seqlen;i++)
		{
			if(i%16==0 && i>0)fprintf(fp,"\n");
			string str=One2Three_III(seqres[i]);
			fprintf(fp,"%s ",str.c_str());
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
}

