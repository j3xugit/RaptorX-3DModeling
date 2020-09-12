#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
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

//------- int to string -------//
string IntToString(int Number)
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}

//--------- Parse_Str_Str --------//
int Parse_Str_Str(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		out.push_back(buf);
		count++;
	}
	return count;
}
//--------- Parse_Str_Str (automatic) --------//
int Parse_Str_Str(string &in,vector <string> &out)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(! (www>>buf) )break;
		out.push_back(buf);
		count++;
	}
	return count;
}

//------ remove CHAR --------//
void remove_CHAR(string &in, string &out, char c)
{
	out="";
	for(int i=0;i<(int)in.length();i++)
	{
		if(in[i]!=c)out.push_back(in[i]);
	}
}


//=============== Load MSA, output to FASTA =============//
//->
/*
115422923|10_157|1_163
115526545|16_162|1_171
121594756|65_214|45_222
121603993|35_186|15_194
13476734|19_171|1_191
*/
int Multi_FASTA_Input(string &multi_fasta,
	vector <string> &nam_list,vector <string> &fasta_list,
	vector <string> &header_list)
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
	header_list.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf.length()>=1 && buf[0]=='>')
		{
			if(first!=1)
			{
				fasta_list.push_back(seq);
				number++;
			}
			temp=buf.substr(1,buf.length()-1);
			header_list.push_back(temp);
			//-> the following foramt is based on uniprot_human_nr
			//vector <string> out_str;
			//int retv=Parse_Str_Str(temp,out_str);
			//-> get name string
			//name=out_str[0];
			string int_nam = IntToString(number);
			nam_list.push_back(int_nam);
			count++;
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
	//---- MSA_Output ----//__2014_10_30__//
	{
		if(argc<3)
		{
			fprintf(stderr,"Version: 1.00 \n");
			fprintf(stderr,"./MSA_To_SEQ <input_multi_seq> <out_seq> \n");
			exit(-1);
		}
		//input
		string uniprot_fasta_msa=argv[1];
		string out_seq=argv[2];

		//read
		int retv;
		vector <string> nam_list;
		vector <string> fasta_list;
		vector <string> header_list;
		retv=Multi_FASTA_Input(uniprot_fasta_msa,nam_list,fasta_list,header_list);

		//output string
		FILE *fp=fopen(out_seq.c_str(),"wb");
		for(int i=0;i<retv;i++)
		{
			string str1,str2;
			remove_CHAR(fasta_list[i],str1,'.');
			remove_CHAR(str1,str2,'-');
			toUpperCase(str2);
			fprintf(fp,">%s\n",header_list[i].c_str());
			fprintf(fp,"%s\n",str2.c_str());
		}
		fclose(fp);

		//exit
		exit(0);
	}
}

