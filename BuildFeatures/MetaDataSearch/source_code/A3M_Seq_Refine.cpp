#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <string.h>
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

//----- get upper case -----//
int getUpperCase(char *buffer)
{
	int count=0;
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=65 && buffer[i]<=90) count++;
	return count;
}
int getUpperCase(string &buffer)
{
	int count=0;
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=65 && buffer[i]<=90) count++;
	return count;
}
//----- get lower case -----//
int getLowerCase(char *buffer)
{
	int count=0;
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=97 && buffer[i]<=122) count++;
	return count;
}
int getLowerCase(string &buffer)
{
	int count=0;
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=97 && buffer[i]<=122) count++;
	return count;
}

//----- trim white space -----//
int Trim_White_Space(string &in)
{
	int count=0;
	string out="";
	for(int i=0;i<(int)in.length();i++)
	if(in[i]!=' ')out+=in[i],count++;
	in=out;
	return count;
}


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


//-------- read in MSA in a3m format (i.e., normal FASTA with upper/lower) ------------//
//[note]: we set the first sequence as the query sequence,
//        that is to say, all the following sequences should be longer than the first
int Multi_FASTA_Input(string &multi_fasta,vector <string> &nam_list,vector <string> &fasta_list)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		return -1;
	}
	//load
	int relfirst=1;
	int firstlen;
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
			name=buf.substr(1,buf.length()-1);
			nam_list.push_back(name);
			count++;
			if(first!=1)
			{
				fasta_list.push_back(seq);
				number++;
				if(number==1)
				{
					firstlen=(int)seq.length();
				}
				else
				{
					Trim_White_Space(seq);
					int lowlen=getLowerCase(seq);
					int curlen=(int)seq.length()-lowlen;
					if(curlen!=firstlen)
					{
						fprintf(stderr,"length not equal at %s, [%d!=%d] \n",buf.c_str(),curlen,firstlen);
						return -1;
					}
				}
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
		fasta_list.push_back(seq);
		number++;
		if(number==1)
		{
			firstlen=(int)seq.length();
		}
		else
		{
			Trim_White_Space(seq);
			int lowlen=getLowerCase(seq);
			int curlen=(int)seq.length()-lowlen;
			if(curlen!=firstlen)
			{
				fprintf(stderr,"length not equal at %s, [%d!=%d] \n",buf.c_str(),curlen,firstlen);
				return -1;
			}
		}
	}
	//check
	if(number!=count)
	{
		fprintf(stderr,"num %d != count %d \n",number,count);
		return -1;
	}
	return count;
}


//--------- read FASTA sequence --------//
int Read_FASTA_SEQRES(string &seqfile,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(seqfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",seqfile.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",seqfile.c_str());
			exit(-1);
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	seqres=temp;
	//return
	return (int)seqres.length();
}



//----- eliminate lower case -------//
void Eliminate_LowerCase(string &instr,string &outstr)
{
	int i;
	outstr.clear();
	for(i=0;i<(int)instr.length();i++)
	{
		if(instr[i]>='a' && instr[i]<='z') continue;
		outstr.push_back(instr[i]);
	}
}


//========= validate sequence ==========//
int Ori_AA_Map[26]=
{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};
// A B C D E F G H I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
// 0 1 2 3 4 5 6 7 8  9 10 11 12 14 14 15 16 17 18 19 20 21 22 23 24 25

void Validate_Sequence(string &instr,string &outstr)
{
	int i;
	int len=(int)instr.length();
	outstr=instr;
	for(i=0;i<len;i++)
	{
		if(instr[i]=='-')continue;
		char a=instr[i];
		if(a<'A' || a>='Z')
		{
			outstr[i]='X';
			continue;
		}
		int retv=Ori_AA_Map[a-'A'];
		if(retv==20)
		{
			outstr[i]='X';
			continue;
		}
	}
}

//-------- main -------//
int main(int argc,char **argv)
{
	//------ A3M_Seq_Refine -------//
	{
		if(argc<4)
		{
			fprintf(stderr,"A3M_Seq_Refine <a3m_input> <seq_input> <a3m_output> \n");
			exit(-1);
		}
		string a3m_input=argv[1];
		string seq_input=argv[2];
		string a3m_output=argv[3];
		//-> get base name
		string seq_name;
		getBaseName(seq_input,seq_name,'/','.');
		//-> load A3M
		vector <string> nam_list;
		vector <string> fasta_list;
		int totnum=Multi_FASTA_Input(a3m_input,nam_list,fasta_list);
		if(totnum<=0)exit(-1);
		//-> load SEQ
		string seqres;
		int retv=Read_FASTA_SEQRES(seq_input,seqres);
		if(retv<=0)exit(-1);
		//-> check
		size_t found = seqres.find(fasta_list[0]);
		if (found==string::npos)
		{
			fprintf(stderr,"seq NOT found in a3m \n");
			exit(-1);
		}
		//-> get head and tail
		int head=found;
		int tail=seqres.length()-fasta_list[0].length()-head;
		string head_str(head, '-');
		string tail_str(tail, '-');
		//-> output
		FILE *fp=fopen(a3m_output.c_str(),"wb");
		fprintf(fp,">%s\n",seq_name.c_str());
		fprintf(fp,"%s\n",seqres.c_str());
		for(int i=1;i<totnum;i++)
		{
			fprintf(fp,">%s\n",nam_list[i].c_str());
			fprintf(fp,"%s%s%s\n",
				head_str.c_str(),fasta_list[i].c_str(),tail_str.c_str());
		}
		fclose(fp);
		//exit
		exit(0);
	}
}

