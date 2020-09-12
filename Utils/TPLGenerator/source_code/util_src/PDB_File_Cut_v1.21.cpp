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

//----- just utility ----//
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++) 
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}

//---------- dynamic programming ----------//
int WWW_Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <int> D[3];      // the path (directions) matrix
	vector <double> M[3];   // the current scores (values) matrix
	//resize(m,n)
	for (i = 0; i < 3; ++i) 
	{
		D[i].resize(m*n);
		M[i].resize(m*n);
	}
	//init()
	double MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = MIN;
	M[_V_][0*DP_maximal+ 0] = MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = MIN;
		M[_H_][i*DP_maximal+ 0] = MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		for (j = 1; j < n; j++) 
		{
			//condition upper
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			v1 = M[_V_][(i-1)*DP_maximal+ j] + gap_ext;
			v2 = M[_S_][(i-1)*DP_maximal+ j] + gap_open;
			v3 = M[_H_][(i-1)*DP_maximal+ j] + gap_open;
			M[_V_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_V_][i*DP_maximal+ j] == v1) D[_V_][i*DP_maximal+ j] = _V_;
			else if(M[_V_][i*DP_maximal+ j] == v2) D[_V_][i*DP_maximal+ j] = _S_;
			else D[_V_][i*DP_maximal+ j] = _H_;
			//condition left
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			v1 = M[_H_][i*DP_maximal+ j-1] + gap_ext;
			v2 = M[_S_][i*DP_maximal+ j-1] + gap_open;
			v3 = M[_V_][i*DP_maximal+ j-1] + gap_open;
			M[_H_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_H_][i*DP_maximal+ j] == v1) D[_H_][i*DP_maximal+ j] = _H_;
			else if(M[_H_][i*DP_maximal+ j] == v2) D[_H_][i*DP_maximal+ j] = _S_;
			else D[_H_][i*DP_maximal+ j] = _V_;
			//condition diag
			dist = score.at((i-1)*IN_maximal+ j-1);  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = M[_V_][(i-1)*DP_maximal+ j-1] + dist;
			v2 = M[_H_][(i-1)*DP_maximal+ j-1] + dist;
			v3 = M[_S_][(i-1)*DP_maximal+ j-1] + dist;
			M[_S_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_S_][i*DP_maximal+ j] == v3) D[_S_][i*DP_maximal+ j] = _S_;
			else if (M[_S_][i*DP_maximal+ j] == v1) D[_S_][i*DP_maximal+ j] = _V_;
			else D[_S_][i*DP_maximal+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=M[_V_][i*DP_maximal+ j];
	v2=M[_H_][i*DP_maximal+ j];
	v3=M[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		pre_case=D[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<int,int>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<int,int>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<int,int>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<int,int>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<int,int>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
}
int Seqres_DynaProg(string &seqres,string &ami_,vector <int> &mapping)
{
	//--[0]check
	int len=(int)seqres.length();
	int totnum=(int)ami_.length();

	//--[1]dynamic_programming
	int i,j;
	int head=0;
	int n1=len;    //SEQRES
	int n2=totnum; //ATOM
	vector <double> score;
	score.resize(len*totnum);
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			if(seqres[i]==ami_[j+head])score.at(i*n2+j)=10;
			else score.at(i*n2+j)=-15;
			{
//				if(seqres[i]=='X'||seqres[i]=='Z'||seqres[i]=='.')score.at(i*n2+j)=0;
//				else if(ami_[j+head]=='X'||ami_[j+head]=='Z'||ami_[j+head]=='.')score.at(i*n2+j)=0;
//				else score.at(i*n2+j)=-15;
			}
		}
	}
	double sco;
	int matchs;
	vector<pair<int,int> > WWW_alignment;
	matchs=WWW_Advance_Align_Dyna_Prog_Double(n1,n2,score,-11,-1,-110,-10,0,0,-110,-110,
		WWW_alignment,sco);
	int lcmp=(int)WWW_alignment.size();
	
	//extract
	mapping.resize(len);
	for(i=0;i<len;i++)mapping[i]=0; //default: NO
	int first,second;
	int retv=1;
	for(i=0;i<lcmp;i++)
	{
		first=WWW_alignment[i].first;
		second=WWW_alignment[i].second;
		if(first<=0)
		{
			if(second>0)
			{
//				fprintf(stderr,"Alignment BAD!!!\n");
//				return -1;
				retv=-1;
				continue;
			}
		}
		if(first>0 && second>0)
		{
			mapping[first-1]=1;
		}
	}
	return retv;
}

//========== process PDB ============//
char WWW_Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(1.0*26,1.0*i);
	//switch
	switch(result)
	{
		case 286:return 'A';
		case 4498:return 'R';
		case 9256:return 'N';
		case 10608:return 'D';
		case 12794:return 'C';
		case 9080:return 'Q';
		case 13812:return 'E';
		case 16516:return 'G';
		case 12383:return 'H';
		case 2998:return 'I';
		case 13635:return 'L';
		case 12803:return 'K';
		case 12960:return 'M';
		case 2901:return 'F';
		case 9921:return 'P';
		case 11614:return 'S';
		case 11693:return 'T';
		case 10601:return 'W';
		case 12135:return 'Y';
		case 7457:return 'V';
		default:return 'X';
	}
}
void Read_PDB_SEQRES(string &mapfile,string &seqres,map<string, int > &_mapping) //->from .pdb file
{
	//--- list for mapping ---//
	map<string, int>::iterator iter;
	_mapping.clear();
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(mapfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",mapfile.c_str());
		exit(-1);
	}
	int len;
	int count=0;
	char cur_chain;
	char rel_chain;
	int first=1;
	char chain='_';
	seqres="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM" && temp!="HETA")continue;
		//check CA
		temp=buf.substr(13,2);
		if(temp!="CA")continue;
		//check DNA
//		if(buf[17]==' ' || buf[18]==' ' || buf[19]==' ')continue;
		//chain
		cur_chain=buf[21];
		if(chain!='_')
		{
			if(first==1)
			{
				if(cur_chain!=chain)continue;
				else
				{
					first=0;
					rel_chain=cur_chain;
				}
			}
			else
			{
				if(cur_chain!=rel_chain)break;
			}
		}
		else
		{
			if(first==1)
			{
				first=0;
				rel_chain=cur_chain;
			}
			else
			{
				if(cur_chain!=rel_chain)break;
			}
		}
		//record name
		name=buf.substr(21,6);
		iter = _mapping.find(name);
		if(iter != _mapping.end())continue;
		count++;
		_mapping.insert(map < string, int >::value_type(name, count));
		//final
		temp=buf.substr(17,3);
		char c=WWW_Three2One_III(temp.c_str());
		seqres=seqres+c;
	}
}
void Output_PDB_Full(string &pdb,map<string, int > &_mapping,vector <int> &_rec,FILE *fp)
{
	//read
	ifstream fin;
	string buf,temp;
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		exit(-1);
	}
	//output
	int len;
	int key;
	string name;
	int first=1;
	char cur_chain;
	char rel_chain;
	char chain='_';
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM" && temp!="HETA")continue;
		//check CA
//		temp=buf.substr(13,2);
//		if(temp!="CA")continue;
		//check DNA
//		if(buf[17]==' ' || buf[18]==' ' || buf[19]==' ')continue;
		//chain
		cur_chain=buf[21];
		if(chain!='_')
		{
			if(first==1)
			{
				if(cur_chain!=chain)continue;
				else
				{
					first=0;
					rel_chain=cur_chain;
				}
			}
			else
			{
				if(cur_chain!=rel_chain)break;
			}
		}
		else
		{
			if(first==1)
			{
				first=0;
				rel_chain=cur_chain;
			}
			else
			{
				if(cur_chain!=rel_chain)break;
			}
		}
		//check name
		name=buf.substr(21,6);
		key=_mapping[name];
		if(key<0)continue;
		key=key-1;
		if(_rec[key]==0)continue; //gapped
		//output
		fprintf(fp,"%s\n",buf.c_str());
	}
}


//=============== PDB_File_Cut =============//
void Read_FASTA_SEQRES(string &mapfile,string &seqres,int skip) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(mapfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		printf("no such file! %s \n",mapfile.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			printf("file bad! %s \n",mapfile.c_str());
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
}
void PDB_File_Cut(string &pdbfile,string &seqfile,string &outfile,int skip)
{
	//-> read fasta sequence
	string fasta_seqres;
	Read_FASTA_SEQRES(seqfile,fasta_seqres,skip);
	//-> read pdb sequence
	string pdb_seqres;
	map<string, int > _mapping;
	Read_PDB_SEQRES(pdbfile,pdb_seqres,_mapping);
	//-> do dynamic programming
	int pdb_length=(int)pdb_seqres.length();
	int seq_length=(int)fasta_seqres.length();
	if(seq_length>pdb_length)
	{
		fprintf(stderr,"over range bad! %s \n",pdbfile.c_str());
//		exit(-1);
	}
	vector <int> mapping;
	int retv=Seqres_DynaProg(pdb_seqres,fasta_seqres,mapping);
	if(retv!=1)
	{
		fprintf(stderr,"mapping bad! %s \n",pdbfile.c_str());
//		exit(-1);
	}
	//-> output cut file
	FILE *fp=fopen(outfile.c_str(),"wb");
	Output_PDB_Full(pdbfile,_mapping,mapping,fp);
	fclose(fp);
}


//---------- main ----------//
int main(int argc,char **argv)
{
	//---- process HMM template ----//
	{
		if(argc<5)
		{
			printf("Version: 1.21\n");
			printf("PDB_File_Cut <pdbfile> <fasta> <outfile> <skip>\n");
			exit(-1);
		}
		string infile=argv[1];
		string fasta=argv[2];
		string outfile=argv[3];
		int skip=atoi(argv[4]);
		PDB_File_Cut(infile,fasta,outfile,skip);
		exit(0);
	}
}
