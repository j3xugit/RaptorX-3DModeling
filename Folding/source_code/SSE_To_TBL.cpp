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

//---- parse feat_num ---//
int Parse_FeatNum(string &in_str)
{
	int feat_num=0;
	istringstream www(in_str);
	for(;;)
	{
		string temp;
		if(! (www>>temp) )break;
		feat_num++;
	}
	return feat_num;
}

//================ load secondary structure string ================//
//-> example
/*
>106952
CCCHHHHHHHHHHHHHCCCCCCCCCCEEECCCCCCHHHHHHHHHHHHHHCCCCCCCEEEEECCCCCCCHHHCCCCCCCHHHHHHHHCCCCCCCCCCCHHHHHHHHHHHHCCCCCCCCCCEEEECCCCCCCCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHHCCCCCCEEEEECHHHHHHHHHHHCCCCEEEEEEEEECCCCCCCCCCCCCHHHHHHHHHHHHHHHHCCHHHCCCCCCC
*/
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

//================ SSE string to Start,End block ============//
int SSE_To_StartEnd(string &sse_str, vector < pair<int,int> > & start_end_block, char TargChar)
{
	int i;
	int length=(int)sse_str.length();
	start_end_block.clear();
	int start,end,len;
	int first=1;
	int count=0;
	for(i=0;i<length;i++)
	{
		if(sse_str[i]==TargChar)
		{
			if(first==1)
			{
				first=0;
				start=i;
				len=1;
			}
			else
			{
				len++;
			}
		}
		else
		{
			if(first==0)
			{
				first=1;
				end=start+len-1;
				//push_back
				start_end_block.push_back(pair<int,int>(start,end));
				count++;
			}
		}
	}
	//final assign
	if(first==0)
	{
		first=1;
		end=start+len-1;
		//push_back
		start_end_block.push_back(pair<int,int>(start,end));
		count++;
	}
	//return
	return count;
}


//================ output ssnoe.tbl =============//->SS_noe
//-> example
/*
assign (resid   4 and name  C) (resid   9 and name  C) 8.16 0.10 0.10 !helix
assign (resid   4 and name  C) (resid   7 and name  C) 4.87 0.10 0.10 !helix
assign (resid   4 and name  C) (resid   8 and name  C) 6.09 0.10 0.10 !helix
assign (resid   4 and name CA) (resid   9 and name CA) 8.63 0.10 0.10 !helix
assign (resid   4 and name CA) (resid   7 and name CA) 5.13 0.10 0.10 !helix
assign (resid   4 and name CA) (resid   8 and name CA) 6.16 0.10 0.10 !helix
assign (resid   4 and name  N) (resid   9 and name  N) 8.07 0.10 0.10 !helix
assign (resid   4 and name  N) (resid   7 and name  N) 4.84 0.10 0.10 !helix
assign (resid   4 and name  N) (resid   8 and name  N) 6.10 0.10 0.10 !helix
assign (resid   4 and name  O) (resid   9 and name  O) 8.40 0.10 0.10 !helix
assign (resid   4 and name  O) (resid   7 and name  O) 4.99 0.10 0.10 !helix
assign (resid   4 and name  O) (resid   8 and name  O) 6.12 0.10 0.10 !helix
assign (resid   5 and name  C) (resid  10 and name  C) 8.16 0.10 0.10 !helix
assign (resid   5 and name  C) (resid   8 and name  C) 4.87 0.10 0.10 !helix
assign (resid   5 and name  C) (resid   9 and name  C) 6.09 0.10 0.10 !helix
assign (resid   5 and name CA) (resid  10 and name CA) 8.63 0.10 0.10 !helix
assign (resid   5 and name CA) (resid   8 and name CA) 5.13 0.10 0.10 !helix
assign (resid   5 and name CA) (resid   9 and name CA) 6.16 0.10 0.10 !helix
assign (resid   5 and name  N) (resid  10 and name  N) 8.07 0.10 0.10 !helix
assign (resid   5 and name  N) (resid   8 and name  N) 4.84 0.10 0.10 !helix
assign (resid   5 and name  N) (resid   9 and name  N) 6.10 0.10 0.10 !helix
assign (resid   5 and name  O) (resid  10 and name  O) 8.40 0.10 0.10 !helix
assign (resid   5 and name  O) (resid   8 and name  O) 4.99 0.10 0.10 !helix
assign (resid   5 and name  O) (resid   9 and name  O) 6.12 0.10 0.10 !helix
....
assign (resid 228 and name  C) (resid 231 and name  C) 4.87 0.10 0.10 !helix
assign (resid 228 and name  C) (resid 232 and name  C) 6.09 0.10 0.10 !helix
assign (resid 228 and name CA) (resid 231 and name CA) 5.13 0.10 0.10 !helix
assign (resid 228 and name CA) (resid 232 and name CA) 6.16 0.10 0.10 !helix
assign (resid 228 and name  N) (resid 231 and name  N) 4.84 0.10 0.10 !helix
assign (resid 228 and name  N) (resid 232 and name  N) 6.10 0.10 0.10 !helix
assign (resid 228 and name  O) (resid 231 and name  O) 4.99 0.10 0.10 !helix
assign (resid 228 and name  O) (resid 232 and name  O) 6.12 0.10 0.10 !helix
....
assign (resid  27 and name  O) (resid  28 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid  28 and name  O) (resid  29 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid  57 and name  O) (resid  58 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid  58 and name  O) (resid  59 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid  59 and name  O) (resid  60 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid  60 and name  O) (resid  61 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 120 and name  O) (resid 121 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 121 and name  O) (resid 122 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 122 and name  O) (resid 123 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 174 and name  O) (resid 175 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 175 and name  O) (resid 176 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 176 and name  O) (resid 177 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 177 and name  O) (resid 178 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 195 and name  O) (resid 196 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 196 and name  O) (resid 197 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 197 and name  O) (resid 198 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 198 and name  O) (resid 199 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 199 and name  O) (resid 200 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 200 and name  O) (resid 201 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 201 and name  O) (resid 202 and name  O) 4.57 0.10 0.10 !unpaired E residue
assign (resid 202 and name  O) (resid 203 and name  O) 4.57 0.10 0.10 !unpaired E residue
*/

//--- rules for SSNOE ---//
//-> 1. helix
/*
head part contains 3x restraints, tail part contains 2x restraints,
where each restraint includes C,CA,N,O

here 3x is i+3,i+4,i+5
2x is i+3,i+4

[note]: the minimal length for SSNOE_helix is five (actually, this is questionable. Maybe we can try length four later)
*/
//-> 2. sheet
/*
K length sheet will have (K-1) restraints
*/

//----------- SSnoe helix -----------//
void Output_SSnoe_Helix(vector < pair<int,int> > & start_end_block, vector <string> &output)
{
	int i,j;
	int size=(int)start_end_block.size();
	char command[30000];
	string temp;
	for(i=0;i<size;i++)
	{
		int start=start_end_block[i].first;
		int end=start_end_block[i].second;
		int len=end-start+1;
		//length check
		if(len<5)continue;
		//assign head
		for(j=0;j<len-5;j++)
		{
			//--- triplet ---//
			//-> C
			sprintf(command,"assign (resid %3d and name  C) (resid %3d and name  C) 8.16 0.10 0.10 !helix",start+j+1,start+j+1+5);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  C) (resid %3d and name  C) 4.87 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  C) (resid %3d and name  C) 6.09 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
			//-> CA
			sprintf(command,"assign (resid %3d and name CA) (resid %3d and name CA) 8.63 0.10 0.10 !helix",start+j+1,start+j+1+5);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name CA) (resid %3d and name CA) 5.13 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name CA) (resid %3d and name CA) 6.16 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
			//-> N
			sprintf(command,"assign (resid %3d and name  N) (resid %3d and name  N) 8.07 0.10 0.10 !helix",start+j+1,start+j+1+5);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  N) (resid %3d and name  N) 4.84 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  N) (resid %3d and name  N) 6.10 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
			//-> O
			sprintf(command,"assign (resid %3d and name  O) (resid %3d and name  O) 8.40 0.10 0.10 !helix",start+j+1,start+j+1+5);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  O) (resid %3d and name  O) 4.99 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  O) (resid %3d and name  O) 6.12 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
		}
		//assign tail
		j=len-5;
		{
			//--- doublet ---//
			//-> C
			sprintf(command,"assign (resid %3d and name  C) (resid %3d and name  C) 4.87 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  C) (resid %3d and name  C) 6.09 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
			//-> CA
			sprintf(command,"assign (resid %3d and name CA) (resid %3d and name CA) 5.13 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name CA) (resid %3d and name CA) 6.16 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
			//-> N
			sprintf(command,"assign (resid %3d and name  N) (resid %3d and name  N) 4.84 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  N) (resid %3d and name  N) 6.10 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
			//-> O
			sprintf(command,"assign (resid %3d and name  O) (resid %3d and name  O) 4.99 0.10 0.10 !helix",start+j+1,start+j+1+3);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name  O) (resid %3d and name  O) 6.12 0.10 0.10 !helix",start+j+1,start+j+1+4);
			temp=command;
			output.push_back(temp);
		}
	}
}
//----------- SSnoe sheet -----------//
void Output_SSnoe_Sheet(vector < pair<int,int> > & start_end_block, vector <string> &output)
{
	int i,j;
	int size=(int)start_end_block.size();
	char command[30000];
	string temp;
	for(i=0;i<size;i++)
	{
		int start=start_end_block[i].first;
		int end=start_end_block[i].second;
		int len=end-start+1;
		//length check
		if(len<2)continue;
		//assign head
		for(j=0;j<len-1;j++)
		{
			//-> O
			sprintf(command,"assign (resid %3d and name  O) (resid %3d and name  O) 4.57 0.10 0.10 !unpaired E residue",start+j+1,start+j+1+1);
			temp=command;
			output.push_back(temp);
		}
	}
}


//================ output dihedral.tbl =============//->Dihedral Angle
//---- helix ----//(minimal length is one, and different for 'head' and 'tail')
//-> five helix
/*
assign (resid  98 and name c) (resid  99 and name  n) (resid  99 and name ca) (resid  99 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  99 and name n) (resid  99 and name ca) (resid  99 and name  c) (resid 100 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  99 and name c) (resid 100 and name  n) (resid 100 and name ca) (resid 100 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid 100 and name n) (resid 100 and name ca) (resid 100 and name  c) (resid 101 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid 100 and name c) (resid 101 and name  n) (resid 101 and name ca) (resid 101 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid 101 and name n) (resid 101 and name ca) (resid 101 and name  c) (resid 102 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid 101 and name c) (resid 102 and name  n) (resid 102 and name ca) (resid 102 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid 102 and name n) (resid 102 and name ca) (resid 102 and name  c) (resid 103 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid 102 and name c) (resid 103 and name  n) (resid 103 and name ca) (resid 103 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid 103 and name n) (resid 103 and name ca) (resid 103 and name  c) (resid 104 and name n) 5.0  -41.51   3.936 2 !helix psi
*/
//-> four helix
/*
assign (resid  25 and name c) (resid  26 and name  n) (resid  26 and name ca) (resid  26 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  26 and name n) (resid  26 and name ca) (resid  26 and name  c) (resid  27 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  26 and name c) (resid  27 and name  n) (resid  27 and name ca) (resid  27 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  27 and name n) (resid  27 and name ca) (resid  27 and name  c) (resid  28 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  27 and name c) (resid  28 and name  n) (resid  28 and name ca) (resid  28 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  28 and name n) (resid  28 and name ca) (resid  28 and name  c) (resid  29 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  28 and name c) (resid  29 and name  n) (resid  29 and name ca) (resid  29 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  29 and name n) (resid  29 and name ca) (resid  29 and name  c) (resid  30 and name n) 5.0  -41.51   3.936 2 !helix psi
*/
//-> three helix
/*
assign (resid  25 and name c) (resid  26 and name  n) (resid  26 and name ca) (resid  26 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  26 and name n) (resid  26 and name ca) (resid  26 and name  c) (resid  27 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  26 and name c) (resid  27 and name  n) (resid  27 and name ca) (resid  27 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  27 and name n) (resid  27 and name ca) (resid  27 and name  c) (resid  28 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  27 and name c) (resid  28 and name  n) (resid  28 and name ca) (resid  28 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  28 and name n) (resid  28 and name ca) (resid  28 and name  c) (resid  29 and name n) 5.0  -41.51   3.936 2 !helix psi
*/
//-> two helix
/*
assign (resid  25 and name c) (resid  26 and name  n) (resid  26 and name ca) (resid  26 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  26 and name n) (resid  26 and name ca) (resid  26 and name  c) (resid  27 and name n) 5.0  -41.51   3.936 2 !helix psi
assign (resid  26 and name c) (resid  27 and name  n) (resid  27 and name ca) (resid  27 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  27 and name n) (resid  27 and name ca) (resid  27 and name  c) (resid  28 and name n) 5.0  -41.51   3.936 2 !helix psi
*/
//-> one helix
/*
assign (resid  25 and name c) (resid  26 and name  n) (resid  26 and name ca) (resid  26 and name c) 5.0  -63.47    3.68 2 !helix phi
assign (resid  26 and name n) (resid  26 and name ca) (resid  26 and name  c) (resid  27 and name n) 5.0  -41.51   3.936 2 !helix psi
*/
//-> head helix
/*
assign (resid   1 and name n) (resid   1 and name ca) (resid   1 and name  c) (resid   2 and name n) 5.0  -41.51   3.936 2 !helix psi
*/
//-> tail helix
/*
assign (resid 122 and name c) (resid 123 and name  n) (resid 123 and name ca) (resid 123 and name c) 5.0  -63.47    3.68 2 !helix phi
*/

//------ sheet -----// (minimal length is two)
//-> long sheet
/*
assign (resid   3 and name n) (resid   3 and name ca) (resid   3 and name  c) (resid   4 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   3 and name c) (resid   4 and name  n) (resid   4 and name ca) (resid   4 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
assign (resid   4 and name n) (resid   4 and name ca) (resid   4 and name  c) (resid   5 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   4 and name c) (resid   5 and name  n) (resid   5 and name ca) (resid   5 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
assign (resid   5 and name n) (resid   5 and name ca) (resid   5 and name  c) (resid   6 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   5 and name c) (resid   6 and name  n) (resid   6 and name ca) (resid   6 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
assign (resid   6 and name n) (resid   6 and name ca) (resid   6 and name  c) (resid   7 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   6 and name c) (resid   7 and name  n) (resid   7 and name ca) (resid   7 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
assign (resid   7 and name n) (resid   7 and name ca) (resid   7 and name  c) (resid   8 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   7 and name c) (resid   8 and name  n) (resid   8 and name ca) (resid   8 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
assign (resid   8 and name n) (resid   8 and name ca) (resid   8 and name  c) (resid   9 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   8 and name c) (resid   9 and name  n) (resid   9 and name ca) (resid   9 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
assign (resid   9 and name n) (resid   9 and name ca) (resid   9 and name  c) (resid  10 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid   9 and name c) (resid  10 and name  n) (resid  10 and name ca) (resid  10 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
*/
//-> two sheet
/*
assign (resid  17 and name n) (resid  17 and name ca) (resid  17 and name  c) (resid  18 and name n) 5.0  134.95    7.06 2 !unpaired E residue psi
assign (resid  17 and name c) (resid  18 and name  n) (resid  18 and name ca) (resid  18 and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi
*/

//----------- dihedral helix -----------//
void Output_Dihedral_Helix(string &seqres,
	vector < pair<int,int> > & start_end_block, vector <string> &output)
{
	int i,j;
	int length=(int)seqres.length();
	int size=(int)start_end_block.size();
	char command[30000];
	string temp;
	for(i=0;i<size;i++)
	{
		int start=start_end_block[i].first;
		int end=start_end_block[i].second;
		int len=end-start+1;
		//---- note that head and tail are special cases ----//
		for(j=0;j<len;j++)
		{
			if(start+j+0>0)
			{
				sprintf(command,"assign (resid %3d and name c) (resid %3d and name  n) (resid %3d and name ca) (resid %3d and name c) 5.0  -63.47    3.68 2 !helix phi",
					start+j+0,start+j+1,start+j+1,start+j+1);
				temp=command;
				output.push_back(temp);
			}
			if(start+j+2<=length)
			{
				sprintf(command,"assign (resid %3d and name n) (resid %3d and name ca) (resid %3d and name  c) (resid %3d and name n) 5.0  -41.51   3.936 2 !helix psi",
					start+j+1,start+j+1,start+j+1,start+j+2);
				temp=command;
				output.push_back(temp);
			}
		}
	}
}

//----------- dihedral sheet -----------//
void Output_Dihedral_Sheet(string &seqres,
	vector < pair<int,int> > & start_end_block, vector <string> &output)
{
	int i,j;
	int length=(int)seqres.length();
	int size=(int)start_end_block.size();
	char command[30000];
	string temp;
	for(i=0;i<size;i++)
	{
		int start=start_end_block[i].first;
		int end=start_end_block[i].second;
		int len=end-start+1;
		//--- length check ---//
		if(len<2)continue;
		for(j=0;j<len-1;j++)
		{
			sprintf(command,"assign (resid %3d and name n) (resid %3d and name ca) (resid %3d and name  c) (resid %3d and name n) 5.0  134.95    7.06 2 !unpaired E residue psi",
				start+j+1,start+j+1,start+j+1,start+j+2);
			temp=command;
			output.push_back(temp);
			sprintf(command,"assign (resid %3d and name c) (resid %3d and name  n) (resid %3d and name ca) (resid %3d and name c) 5.0 -118.91   8.692 2 !unpaired E residue phi",
				start+j+1,start+j+2,start+j+2,start+j+2);
			temp=command;
			output.push_back(temp);
		}
	}
}

//================ output hbond.tbl =============//->SS_hbond
//->helix only with at least five residues
/*
HHHEEEEEEECCCCCCEECCCEEEEHHHHEEEEEECCCCCEEEECCCCCCCCCCCCECCCCCEEEEECCCCEEEEEECCCCCCCCCEEEEEECCCCCCHHHHHCCCCCHHHHHHHHHHHHHHH
*/
/*
assign (resid  99 and name O) (resid 103 and name H) 1.99  0.1  0.1 !helix
assign (resid 109 and name O) (resid 113 and name H) 1.99  0.1  0.1 !helix
assign (resid 110 and name O) (resid 114 and name H) 1.99  0.1  0.1 !helix
assign (resid 111 and name O) (resid 115 and name H) 1.99  0.1  0.1 !helix
assign (resid 112 and name O) (resid 116 and name H) 1.99  0.1  0.1 !helix
assign (resid 113 and name O) (resid 117 and name H) 1.99  0.1  0.1 !helix
assign (resid 114 and name O) (resid 118 and name H) 1.99  0.1  0.1 !helix
assign (resid 115 and name O) (resid 119 and name H) 1.99  0.1  0.1 !helix
assign (resid 116 and name O) (resid 120 and name H) 1.99  0.1  0.1 !helix
assign (resid 117 and name O) (resid 121 and name H) 1.99  0.1  0.1 !helix
assign (resid 118 and name O) (resid 122 and name H) 1.99  0.1  0.1 !helix
assign (resid 119 and name O) (resid 123 and name H) 1.99  0.1  0.1 !helix
*/

//----------- hbond helix -----------//
void Output_HBond_Helix(string &seqres,
	vector < pair<int,int> > & start_end_block, vector <string> &output)
{
	int i,j;
	int size=(int)start_end_block.size();
	char command[30000];
	string temp;
	for(i=0;i<size;i++)
	{
		int start=start_end_block[i].first;
		int end=start_end_block[i].second;
		int len=end-start+1;
		//--- length check ---//
		if(len<5)continue;
		for(j=0;j<len-4;j++)
		{
			if(seqres[start+j+4]=='P')
			{
				sprintf(command,"assign (resid %3d and name O) (resid %3d and name N) 2.99  0.1  0.1 !helix",
					start+j+1,start+j+1+4);
			}
			else
			{
				sprintf(command,"assign (resid %3d and name O) (resid %3d and name H) 1.99  0.1  0.1 !helix",
					start+j+1,start+j+1+4);
			}
			temp=command;
			output.push_back(temp);
		}
	}
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//------- SSE_To_TBL -----//
	{
		if(argc<6)
		{
			fprintf(stderr,"Version: 1.00 (2017-01-17) \n");
			fprintf(stderr,"SSE_To_TBL <seqres_str> <sse_string>  \n");
			fprintf(stderr,"           <noe_file> <angle_file> <hbond_file> \n");
			fprintf(stderr,"[note]: input file <seqres_str> should be FASTA format \n");
			fprintf(stderr,"                    for input Amino Acid sequence. \n");
			fprintf(stderr,"                   <sse_string> should be FASTA format \n");
			fprintf(stderr,"                    for Secondary Structure string. \n");
			fprintf(stderr,"        output files <noe_file> for 'ssnoe.tbl' \n");
			fprintf(stderr,"                     <angle_file> for 'dihedral.tbl' \n");
			fprintf(stderr,"                     <hbond_file> for 'hbond.tbl' \n");
			exit(-1);
		}
		string seqres_file=argv[1];
		string sse_string=argv[2];
		string noe_file=argv[3];
		string angle_file=argv[4];
		string hbond_file=argv[5];

		//input files
		//-> load seqres_str
		string seqres_str;
		int l1=Read_FASTA_SEQRES(seqres_file,seqres_str);
		//-> load sse_string
		string sse_str;
		int l2=Read_FASTA_SEQRES(sse_string,sse_str);
		//-> length check
		if(l1!=l2)
		{
			fprintf(stderr,"seqres length %d not equal to sse length %d \n",l1,l2);
			exit(-1);
		}

		//process sse_string
		//-> get helix blocks
		vector < pair<int,int> > Helix_start_end_block;
		SSE_To_StartEnd(sse_str, Helix_start_end_block, 'H');
		//-> get sheet blocks
		vector < pair<int,int> > Sheet_start_end_block;
		SSE_To_StartEnd(sse_str, Sheet_start_end_block, 'E');

		//output TBL files
		FILE *fp;
		vector <string> out_str;
		//-> output noe_file
		out_str.clear();
		Output_SSnoe_Helix(Helix_start_end_block, out_str);
		Output_SSnoe_Sheet(Sheet_start_end_block, out_str);
		fp=fopen(noe_file.c_str(),"wb");
		for(int i=0;i<(int)out_str.size();i++)fprintf(fp,"%s\n",out_str[i].c_str());
		fclose(fp);
		//-> output angle_file
		out_str.clear();
		Output_Dihedral_Helix(sse_str,Helix_start_end_block, out_str);
		Output_Dihedral_Sheet(sse_str,Sheet_start_end_block, out_str);
		fp=fopen(angle_file.c_str(),"wb");
		for(int i=0;i<(int)out_str.size();i++)fprintf(fp,"%s\n",out_str[i].c_str());
		fclose(fp);
		//-> output hbond_file
		out_str.clear();
		Output_HBond_Helix(seqres_str,Helix_start_end_block, out_str);
		fp=fopen(hbond_file.c_str(),"wb");
		for(int i=0;i<(int)out_str.size();i++)fprintf(fp,"%s\n",out_str[i].c_str());
		fclose(fp);

		//exit
		exit(0);
	}
}
