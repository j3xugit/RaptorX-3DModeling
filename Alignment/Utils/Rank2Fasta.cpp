#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>         // -> openmp related !!
using namespace std;


//============== Head_Tail_Mapping =================//
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


//--------- FASTA I/O ------------//
//FASTA
int ReadToFile_FASTA(string &fn,vector<pair<int, int> > &alignment,
					  string &nam1_content,string &nam2_content,
					  string &nam1_full,string &nam2_full,
					  string &nam1,string &nam2)
{
	int i;
	int cur1=0;
	int cur2=0;
	int len;
	int len1,len2;
	alignment.clear();
	//init
	string seq="";  //sequence
	string tmp="";  //template
	//load
	ifstream fin;
	string buf,temp;
	fin.open(fn.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"alignment file not found [%s] !!!\n",fn.c_str());
		return -1;
	}
	//read tmp
	for(;;)
	{
		if(!getline(fin,buf,'\n'))goto badend;
		len=(int)buf.length();
		if(len>1)
		{
			if(buf[0]=='>')
			{
				istringstream www(buf);
				www>>temp;
				len=(int)temp.length();
				nam1=temp.substr(1,len-1);
				break;
			}
		}
	}
	for(;;)
	{
		if(!getline(fin,buf,'\n'))goto badend;
		len=(int)buf.length();
		if(len==0)continue;
		if(len>1)
		{
			if(buf[0]=='>')
			{
				istringstream www(buf);
				www>>temp;
				len=(int)temp.length();
				nam2=temp.substr(1,len-1);
				break;
			}
		}
		tmp+=buf;
	}
	//read seq
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len==0)continue;
		seq+=buf;
	}
	//process
	len1=(int)seq.length();
	len2=(int)tmp.length();
	if(len1!=len2)
	{
		fprintf(stderr,"alignment len not equal [%s] !!!\n",fn.c_str());
		return -1;
	}
	len=len1;
	nam1_content.clear();
	nam2_content.clear();
	for(i=0;i<len;i++)
	{
		if(tmp[i]!='-' && seq[i]!='-') //match
		{
			nam1_content.push_back(tmp[i]);
			nam2_content.push_back(seq[i]);
			cur1++;
			cur2++;
			alignment.push_back(pair<int,int>(cur1,cur2));
		}
		else
		{
			if(tmp[i]!='-') //Ix
			{
				nam1_content.push_back(tmp[i]);
				cur1++;
				alignment.push_back(pair<int,int>(cur1,-cur2));
			}
			if(seq[i]!='-') //Iy
			{
				nam2_content.push_back(seq[i]);
				cur2++;
				alignment.push_back(pair<int,int>(-cur1,cur2));
			}
		}
	}
	//return
	nam1_full=tmp;
	nam2_full=seq;
	return 1; //success

badend:
	fprintf(stderr,"alignment file format bad [%s] !!!\n",fn.c_str());
	return -1;
}

//------ read fasta sequence -----//
void WS_Read_FASTA_SEQRES(string &mapfile,string &seqres,int skip=1) //->from .fasta file
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
	double WS_MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = WS_MIN;
	M[_V_][0*DP_maximal+ 0] = WS_MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = WS_MIN;
		M[_H_][i*DP_maximal+ 0] = WS_MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = WS_MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = WS_MIN;
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
int WS_Seqres_DynaProg(string &seqres,string &ami_,int *mapping)
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
			else
			{
				if(seqres[i]=='X'||seqres[i]=='Z'||seqres[i]=='.')score.at(i*n2+j)=0;
				else if(ami_[j+head]=='X'||ami_[j+head]=='Z'||ami_[j+head]=='.')score.at(i*n2+j)=0;
				else score.at(i*n2+j)=-15;
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
	for(i=0;i<len;i++)mapping[i]=-1; //default: NO
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
				retv=-1;
				continue;
			}
		}
		if(first>0 && second>0)
		{
			mapping[first-1]=second-1;
		}
	}
	return retv;
}

//------ aliX to alignment ------//
void AliX_To_Alignment(int *ali1,int *ali2,int moln1,int moln2,vector<pair<int, int> > &alignment)
{
	//init
	alignment.clear();
	int length_x=moln1;
	int length_y=moln2;
	//start
	int i;
	int j;
	int wlen;
	int ii,jj;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=length_x;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
			      pre_ii++;
			      alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
			      pre_jj++;
			      alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=length_x;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=length_y;i++)alignment.push_back (pair<int,int>(-length_x, i));  //Iy
}

//======== head_tail mapping ========//
//-> fasta is a cutted one, sequence is a full one
void WS_Head_Tail_Mapping(string &fasta_file,string &sequence_file,string &out_file,int skip)
{
	int i;
	int retv;
	//load alignment file
	string nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2;
	vector<pair<int, int> > alignment;
	retv=ReadToFile_FASTA(fasta_file,alignment,nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2);
	//load sequence file
	string sequence_fasta;
	WS_Read_FASTA_SEQRES(sequence_file,sequence_fasta,skip);
	//do dynamic programming
	int map_len1=(int)sequence_fasta.length();
	int map_len2=(int)nam2_content.length();
	int *mapping=new int[map_len1+map_len2];
	retv=WS_Seqres_DynaProg(sequence_fasta,nam2_content,mapping);
	//create ali
	int moln1=(int)nam1_content.length();
	int moln2=(int)nam2_content.length();
	int moln3=(int)sequence_fasta.length();
	int *ali1=new int[moln1];
	int *ali2=new int[moln2];
	int *ali3=new int[moln3];
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)ali2[i]=-1;
	for(i=0;i<moln3;i++)ali3[i]=-1;
	//get initial alignment
	int ii,jj;
	int size=(int)alignment.size();
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			ali1[ii-1]=jj-1;
			ali2[jj-1]=ii-1;
		}
	}
	//do mapping
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln3;i++)ali3[i]=-1;
	for(i=0;i<moln3;i++)
	{
		if(mapping[i]!=-1)
		{
			jj=mapping[i];
			ii=ali2[jj];
			if(ii>=0)
			{
				ali1[ii]=i;
				ali3[i]=ii;
			}
		}
	}
	//output
	AliX_To_Alignment(ali1,ali3,moln1,moln3,alignment);
	size=(int)alignment.size();
	string out1="";
	string out2="";
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			out1=out1+nam1_content[ii-1];
			out2=out2+sequence_fasta[jj-1];
		}
		else
		{
			if(ii>0)
			{
				out1=out1+nam1_content[ii-1];
				out2=out2+'-';
			}
			if(jj>0)
			{
				out1=out1+'-';
				out2=out2+sequence_fasta[jj-1];
			}
		}
	}
	FILE *fp=fopen(out_file.c_str(),"wb");
	fprintf(fp,">%s\n",nam1.c_str());
	fprintf(fp,"%s\n",out1.c_str());
	fprintf(fp,">%s\n",nam2.c_str());
	fprintf(fp,"%s\n",out2.c_str());
	fclose(fp);
	//delete
	delete [] mapping;
	delete [] ali1;
	delete [] ali2;
	delete [] ali3;
}


//================ build model ====================//
//================= mapping_series ===================//
//---- DynaProg ----//
int Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
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
	double WS_MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = WS_MIN;
	M[_V_][0*DP_maximal+ 0] = WS_MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = WS_MIN;
		M[_H_][i*DP_maximal+ 0] = WS_MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = WS_MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = WS_MIN;
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

//======== PDB_alignment ============//
int analyse_seperate(int i1,char c1,int i2,char c2,int &sep) 
{
	if(c1==' ')
	{
		if(c2==' ')sep=i2-i1;
		else
		{
			if(i1==i2)
			{
				sep=c2-'A';
				return -1;
			}
			else
			{
				sep=0;
				return -1;
			}
		}
	}
	else
	{
		if(c2!=' ')
		{
			if(c1==c2)sep=i2-i1;
			else
			{
				if(i1==i2)sep=abs(c2-c1);
				else
				{
					sep=0;
					return -1;
				}
			}
		}
		else
		{
			sep=0;
			return -1;
		}
	}
	if(sep<0)return -1;
	if(i1*i2<0)
	{
		sep--;
		return 2;  // represent:"0"	
	}
	else return 1;  // normal_return
}
int process_oriami_record(char *seqres,char *ami_,int *int_,char *ins_,char *tag_,
	vector<pair<int,int> > &WWW_alignment,vector <int> &wali1,vector <int> &wali2)
{
	string out;
	string out1;
	int i,j;
	int head=0;
	int len;
	int totnum;
	int ii,jj;
	int seperate;
	int ret_val;
	int ws_rec_num;
	int n1,n2;

	//--[0]check
	len=(int)strlen(seqres);
	totnum=(int)strlen(ami_);

	//--[1]dynamic_programming	
	n1=len;    //SEQRES
	n2=totnum; //ATOM
	vector <double> score;
	score.resize(len*totnum);
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			if(seqres[i]==ami_[j+head])score.at(i*n2+j)=10;
			else
			{
				if(seqres[i]=='X'||seqres[i]=='Z'||seqres[i]=='.')score.at(i*n2+j)=0;
				else if(ami_[j+head]=='X'||ami_[j+head]=='Z'||ami_[j+head]=='.')score.at(i*n2+j)=0;
				else score.at(i*n2+j)=-15;
			}
		}
	}
	double sco;
	int matchs;
	matchs=Advance_Align_Dyna_Prog_Double(n1,n2,score,-11,-1,-110,-10,0,0,0,0,
		WWW_alignment,sco);
	int lcmp=(int)WWW_alignment.size();

	//--Output_DP_Result
	//init
	wali1.resize(n1);
	wali2.resize(n2);
	for(j=0;j<n1;j++)wali1[j]=-1;
	for(j=0;j<n2;j++)wali2[j]=-1;

	//record_neo
	i=0;
	ii=0;
	int first,second; //first SEQUES, second ATOM
	int IsInsert=1;
	for(j=0;j<lcmp;j++)
	{
		first=WWW_alignment.at(j).first;
		second=WWW_alignment.at(j).second;
		if(first<=0)
		{
			if(second>0)
			{
				IsInsert=-1;
				tag_[i+head]='i';
				i++;
			}
		}
		else
		{
			if(second>0)
			{
				wali1.at(ii)=i;  //seqres_ami
				wali2.at(i)=ii;  //atom_ami
				if(seqres[ii]!=ami_[i+head])tag_[i+head]='/';
				i++;
			}
			ii++;
		}
	}
	//bad process
	if(IsInsert==-1)return IsInsert;
	if(matchs<2)return -1;

	//head_tail_tag	
	ii=wali2.at(totnum-1);
	if(ii!=-1)for(i=ii+1;i<len;i++)wali1.at(i)=-2; //tail_tag
	ii=wali2.at(0);
	if(ii!=-1)for(i=0;i<ii;i++)wali1.at(i)=-2; //head_tag

	//analyse_main_backword
	ws_rec_num=0;
	for(i=totnum-1;i>=1;i--)
	{
		if(tag_[i+head]=='i')continue; //__Found_Insert__//__080326__//
		ret_val=analyse_seperate(int_[i-1+head],ins_[i-1+head],int_[i+head],ins_[i+head],seperate);
		ii=wali2.at(i)-seperate; //expected_position
		jj=wali2.at(i-1);        //current_position
		if(ii!=jj && ret_val!=-1 && ii>=0 && ii<len) //error
		{
			if(ws_rec_num>=8)tag_[i+head]*=-1;  //solid!!
			ws_rec_num=0;
			for(j=0;j<ret_val;j++)
			{
				if(ii-j<0)break;
				if(wali1.at(ii-j)==-1)
				{
					if(ami_[i-1+head]==seqres[ii-j])
					{
						if(tag_[i-1+head]=='/')tag_[i-1+head]=' ';
					}
					else if(tag_[i-1+head]!='/')continue;
					if(jj>=0 && jj<n1)wali1.at(jj)=-1;
					wali1.at(ii-j)=i-1;
					wali2.at(i-1)=ii-j;				
					if(seperate==1)tag_[i-1+head]*=-1;  //solid!!
					break;
				}
			}
		}
		else ws_rec_num++;
	}

	//analyse_main_forward
	for(i=0;i<totnum-1;i++)
	{
		if(tag_[i+head]=='i')continue; //__Found_Insert__//__080326__//
		ret_val=analyse_seperate(int_[i+head],ins_[i+head],int_[i+1+head],ins_[i+1+head],seperate);
		ii=wali2.at(i)+seperate; //expected_position
		jj=wali2.at(i+1);        //current_position
		if(ii!=jj && ret_val!=-1 && ii>=0 && ii<len) //error
		{
			if(seperate!=1 && tag_[i+1+head]<0)continue;
			for(j=0;j<ret_val;j++)
			{
				if(ii+j>=len)break;
				if(wali1.at(ii+j)==-1)
				{
					if(ami_[i+1+head]==seqres[ii+j])
					{
						if(tag_[i+1+head]=='/')tag_[i+1+head]=' ';
					}
					else if(tag_[i+1+head]!='/')continue;
					if(jj>=0 && jj<n1)wali1.at(jj)=-1;
					wali1.at(ii+j)=i+1;
					wali2.at(i+1)=ii+j;
					break;
				}
			}
		}
	}

	//[final correction]
	int cur;
	//head_correct
	cur=0;
	ii=wali2.at(cur);     //current
	jj=wali2.at(cur+1)-1; //mapping
	if(ii!=jj && jj>=0 && jj<len)
	{
		if(wali1.at(jj)==-1)
		{
			if(ami_[cur+head]==seqres[jj])
			{
				if(tag_[cur+head]=='/')tag_[cur+head]=' ';
				wali1.at(ii)=-1;
				wali1.at(jj)=cur;
				wali2.at(cur)=jj;
			}
		}
	}
	//tail_correct
	cur=n2-1;
	ii=wali2.at(cur);     //current
	jj=wali2.at(cur-1)+1; //mapping
	if(ii!=jj && jj>=0 && jj<len)
	{
		if(wali1.at(jj)==-1)
		{
			if(ami_[cur+head]==seqres[jj])
			{
				if(tag_[cur+head]=='/')tag_[cur+head]=' ';
				wali1.at(ii)=-1;
				wali1.at(jj)=cur;
				wali2.at(cur)=jj;
			}
		}
	}

	//return
	return 1;
}
//--------- DSSP_Process --------//
char Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(26.0,1.0*i);
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
// -> updated!! //__120330__//
int WS_Process_PDB(string &file,char *ami_,int *int_,char *ins_)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)return -1;
	//record
	int len;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp=="TER " || temp=="END " || temp=="ENDM")break; //this might be modified in the future
		if(temp!="ATOM" && temp!="HETA")continue;
		if(len<55)continue;
		temp=buf.substr(13,2);
		if(temp!="CA")continue;
		temp=buf.substr(17,3);
		ami_[count]=Three2One_III(temp.c_str());
		temp=buf.substr(22,4);
		int_[count]=atoi(temp.c_str());
		ins_[count]=buf[26];
		count++;
	}
	ami_[count]='\0';
	ins_[count]='\0';
	return count;
}
//---- PDB_Mapping ----//
int WS_PDB_Mapping(char *seqres,string &pdbfile,string &seqws,string &dssp)
{
	int maxlen=30000;
	char *ami_=new char[maxlen];
	int *int_=new int[maxlen];
	char *ins_=new char[maxlen];
	char *tag_=new char[maxlen];
	//load pdb
	int pdb_len=WS_Process_PDB(pdbfile,ami_,int_,ins_);
	if(pdb_len<0)
	{
		printf("PDB file error! [%s]\n",pdbfile.c_str());
		delete [] ami_;
		delete [] int_;
		delete [] ins_;
		delete [] tag_;
		return -1;
	}
	//alignment
	int k;
	for(k=0;k<pdb_len;k++)tag_[k]=' ';
	tag_[k]='\0';
	vector<pair<int,int> > WWW_alignment;
	vector <int> wali1;
	vector <int> wali2;
	int IsInsert=process_oriami_record(seqres,ami_,int_,ins_,tag_,WWW_alignment,wali1,wali2);
	if(IsInsert!=1)
	{
//		printf("Mapping Insert Bad!!\n");
//		delete [] ami_;
//		delete [] int_;
//		delete [] ins_;
//		delete [] tag_;
//		exit(-1);
	}
	//output
	int i,j;
	int ii;
	int count;
	int first,second; //first SEQUES, second ATOM
	int lcmp=(int)WWW_alignment.size();
	char wstemp[60000];
	//[1]
	count=0;
	for(j=0;j<lcmp;j++)
	{
		first=WWW_alignment.at(j).first;
		if(first<=0)continue;
		ii=first-1;
		wstemp[count]=seqres[ii];
		count++;
	}
	wstemp[count]='\n';
	seqws=wstemp;
	//[2]
	count=0;
	for(j=0;j<lcmp;j++)
	{
		first=WWW_alignment.at(j).first;
		second=WWW_alignment.at(j).second;
		if(first<=0)continue;
		ii=first-1;
		if(wali1.at(ii)<=-1)
		{
			wstemp[count]='-';
			count++;
		}
		else
		{
			i=wali1.at(ii);
			wstemp[count]=ami_[i];
			if(ami_[i]=='X')wstemp[count]='.';
			count++;
		}
	}
	wstemp[count]='\n';
	dssp=wstemp;
	//final
	delete [] ami_;
	delete [] int_;
	delete [] ins_;
	delete [] tag_;
	//return
	return 1;
}
//PIR
void SaveToFile_PIR(FILE *fp,string &tempnam,string &targnam,vector<pair<int, int> > &alignment,
					string &template_content,string &sequence_content,int cutnum=60)
{
	string nam1=targnam;
	string nam2=tempnam;
	string chainName="";
	if(tempnam.length()==5)chainName=chainName+tempnam[4]; //tempnam must be like 1pdbA
	else chainName="A";
//	chainName=chainName+tempnam[4];

	//WS_Modification//__110510__//
	fprintf(fp,">P1;%s_%s\n",nam2.c_str(),nam1.c_str());
	fprintf(fp,"sequence:%s_%s::::::::",nam2.c_str(),nam1.c_str());
	//output sequence content	
	for(int i=0; i<(int)alignment.size(); i++){
		if (i%cutnum==0) fprintf(fp,"\n");
		if (alignment[i].second <=0) fprintf(fp,"-");
		else fprintf(fp,"%c",sequence_content[alignment[i].second-1]);
	}
	//* is the ending mark 
	fprintf(fp,"*\n\n\n");

	//WS_Modification//__110510__//
	fprintf(fp,">P1;%s\n",nam2.c_str());
	fprintf(fp,"structure:%s::%s::%s::::",nam2.c_str(), chainName.c_str(), chainName.c_str());
	//output template content
	for(int i=0; i<(int)alignment.size(); i++){
		if (i%cutnum==0) fprintf(fp,"\n");
		if (alignment[i].first <=0) fprintf(fp,"-");
		else fprintf(fp,"%c",template_content[alignment[i].first-1]);
	}
	//* is the ending mark 
	fprintf(fp,"*\n\n");
}

//PROCESS
void WS_ALI_To_PIR_Including_MAP(string &tempnam,string &targnam,string &ali_root,string &pdb_file)
{
	string file;
	vector<pair<int, int> > alignment;
	string temp_content,targ_content;
	string temp_full,targ_full;
	string seqres,dssp;
	string nam1,nam2;
	//process
	file=ali_root+"/"+tempnam+"-"+targnam+".ws_tmp_fasta";
	ReadToFile_FASTA(file,alignment,temp_content,targ_content,temp_full,targ_full,nam1,nam2);
	char seqws[60000];
	strcpy(seqws,temp_content.c_str());
	WS_PDB_Mapping(seqws,pdb_file,seqres,dssp);
	//output
	file=tempnam+"-"+targnam+".pir";
	FILE *fp=fopen(file.c_str(),"wb");
	SaveToFile_PIR(fp,tempnam,targnam,alignment,dssp,targ_content);
	fclose(fp);
}

//=============== WS_Full_Pair_Mod ==============//__110730__//
//------------ WS_Python_Write --------------//
void WS_Pair_Modeler_Python_Write(FILE *fp,const char *nam1,const char *nam2,string &alifile,string &pdb_root,int MOL_NUM=1)
{
	//previous items
	fprintf(fp,"# Homology modeling with multiple templates\n");
	fprintf(fp,"from modeller import *              # Load standard Modeller classes\n");
	fprintf(fp,"from modeller.automodel import *    # Load the automodel class\n");
	fprintf(fp,"\n");
	//create new environment
	fprintf(fp,"log.verbose()    # request verbose output\n");
	fprintf(fp,"env = environ()  # create a new MODELLER environment to build this model in\n");
	fprintf(fp,"\n");
	//PDB directory
	fprintf(fp,"# directories for input atom files\n");
	fprintf(fp,"env.io.atom_files_directory = ['%s']\n",pdb_root.c_str());
	fprintf(fp,"\n");
	//MultiStruAlign alifile
	fprintf(fp,"a = automodel(env,\n");
	fprintf(fp,"              alnfile  = '%s', # alignment filename\n",alifile.c_str());
	//structures
	fprintf(fp,"              knowns   = (");
	fprintf(fp,"'%s', ",nam1);
	fprintf(fp,"),     # codes of the templates\n");
	//sequence
	fprintf(fp,"              sequence = ");
	fprintf(fp,"'%s', ",nam2);
	fprintf(fp,")               # code of the target\n");
	//final items
	fprintf(fp,"a.starting_model= 1                 # index of the first model\n");
	fprintf(fp,"a.ending_model  = %d                # index of the last model\n",MOL_NUM);
	fprintf(fp,"                                    # (determines how many models to calculate)\n");
	fprintf(fp,"a.make()                            # do the actual homology modeling\n");
}

//============= WS_Simple_Mod ============//
//only build model, don't output TMscore (for real CASP case)
void WS_Pair_Mod_Single_Fasta_Simp(string &nam1,string &nam2,
							 string &mod_bin,string &pdb_root,string &pir_root,
							 int MOL_NUM=1)
{
	string name;
	string python;
	string ssstemp;
	string command;
	FILE *fp;

	//alignment file
	name=pdb_root+"/"+nam1+".pdb";  //-> this is template PDB file
	WS_ALI_To_PIR_Including_MAP(nam1,nam2,pir_root,name);
//	return;
	name=nam1+"-"+nam2+".pir"; //this is alignment PIR file
	ssstemp=nam1+"_"+nam2;
	python=ssstemp+".py";
	fp=fopen(python.c_str(),"wb");
	WS_Pair_Modeler_Python_Write(fp,nam1.c_str(),ssstemp.c_str(),name,pdb_root,MOL_NUM);
	fclose(fp);
//	return;
	//modeller
	command=mod_bin+" "+python;
	system(command.c_str());

	//--------- multi_model -------//__110730__//
	int i,k;
	char model_nam[5];
	char wscommand[30000];
	//each num process
	for(i=1;i<=MOL_NUM;i++)
	{
		//name_create
		sprintf(model_nam,"%4d",i);
		for(k=0;k<(int)strlen(model_nam);k++)if(model_nam[k]==' ')model_nam[k]='0';
		//delete
		sprintf(wscommand,"rm -f %s.D0000%s",ssstemp.c_str(),model_nam);
		system(wscommand);
		sprintf(wscommand,"rm -f %s.V9999%s",ssstemp.c_str(),model_nam);
		system(wscommand);
	}

	//final delete
	command="rm -f "+ssstemp+".py";
	system(command.c_str());
	command="rm -f "+ssstemp+".log";
	system(command.c_str());
	command="rm -f "+ssstemp+".ini";
	system(command.c_str());
	command="rm -f "+ssstemp+".rsr";
	system(command.c_str());
	command="rm -f "+ssstemp+".sch";
	system(command.c_str());
	//delete temp pir file
	sprintf(wscommand,"rm -f %s",name.c_str());
	system(wscommand);
}

//---------- final build3Dmodel -----------//
int build3Dmodel(string &fasta_file,string &query_name,string &pdb_root,string &mod_bin,int mod_num)
{
	//read fasta
	string nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2;
	vector<pair<int, int> > alignment;
	int retv=ReadToFile_FASTA(fasta_file,alignment,nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2);
	//check fasta
	if(retv!=1)
	{
		fprintf(stderr,"fasta_file %s not found\n",fasta_file.c_str());
		return -1;
	}
	//check first or second
	int swap=0;
	string tmpnam;
	if(nam1==query_name) //first is query (must swap!!)
	{
		swap=1;
		tmpnam=nam2+"-"+nam1+".ws_tmp_fasta";
		FILE *fp=fopen(tmpnam.c_str(),"wb");
		fprintf(fp,">%s\n",nam2.c_str());
		fprintf(fp,"%s\n",nam2_full.c_str());
		fprintf(fp,">%s\n",nam1.c_str());
		fprintf(fp,"%s\n",nam1_full.c_str());
		fclose(fp);
		//assign
		nam1=nam2;
		nam2=query_name;
	}
	else if(nam2==query_name) //second is query
	{
		swap=0;
		tmpnam=nam1+"-"+nam2+".ws_tmp_fasta";
		FILE *fp=fopen(tmpnam.c_str(),"wb");
		fprintf(fp,">%s\n",nam1.c_str());
		fprintf(fp,"%s\n",nam1_full.c_str());
		fprintf(fp,">%s\n",nam2.c_str());
		fprintf(fp,"%s\n",nam2_full.c_str());
		fclose(fp);
		//don't need assign
	}
	else //none is query
	{
		fprintf(stderr,"query_name %s couldn't be found in fasta_file %s \n",query_name.c_str(),fasta_file.c_str());
		return -1;
	}
	//start modeller
	string pir_root=".";
	WS_Pair_Mod_Single_Fasta_Simp(nam1,nam2,mod_bin,pdb_root,pir_root,mod_num);
	//delete tmp files
	string command;
	command="rm -f "+tmpnam;
	system(command.c_str());
	//return
	return 1;
}




//===================== file process =====================//
//----- process alignment -----//
//note: the input file must in post-process format
int WS_Proc_Alignment(string &infile,int topk,string &targ_nam_out,
	vector <string> &temp_rec, vector <string> &targ_rec)
{
	ifstream fin;
	string buf,temp;
	//load
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!!\n",infile.c_str());
		exit(-1);
	}
	string temp_nam,targ_nam;
	string temp_cur,targ_cur;
	string temp_ful,targ_ful;
	string temp_ali,targ_ali;
	int count=0;
	int first=0;
	temp_nam="";
	targ_nam="";
	temp_ful="";
	targ_ful="";
	temp_rec.clear();
	targ_rec.clear();
	//process
	for(;;)
	{
		//readin
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp>>temp_cur>>temp_ali;
		if(!getline(fin,buf,'\n'))break;
		istringstream sss(buf);
		sss>>temp>>targ_cur>>targ_ali;
		//judge
		if( (temp_cur != temp_nam || targ_cur != targ_nam)  && first==1 )
		{
			//print old
			string outnam=temp_nam+"-"+targ_nam_out+".fasta";
			FILE *fp=fopen(outnam.c_str(),"wb");
			fprintf(fp,">%s\n",temp_nam.c_str());
			fprintf(fp,"%s\n",temp_ful.c_str());
			fprintf(fp,">%s\n",targ_nam_out.c_str());
			fprintf(fp,"%s\n",targ_ful.c_str());
			fclose(fp);
			//record
			temp_rec.push_back(temp_nam);
			targ_rec.push_back(targ_nam);
			//clear
			temp_ful="";
			targ_ful="";
			count++;
		}
		//termi
		if(count>=topk)break;
		temp_nam=temp_cur;
		targ_nam=targ_cur;
		temp_ful+=temp_ali;
		targ_ful+=targ_ali;
		first=1;
	}
	//final check
	if(count<topk && first==1) //need to output final
	{
		//print old
		string outnam=temp_cur+"-"+targ_nam_out+".fasta";
		FILE *fp=fopen(outnam.c_str(),"wb");
		fprintf(fp,">%s\n",temp_cur.c_str());
		fprintf(fp,"%s\n",temp_ful.c_str());
		fprintf(fp,">%s\n",targ_nam_out.c_str());
		fprintf(fp,"%s\n",targ_ful.c_str());
		fclose(fp);
		//record
		temp_rec.push_back(temp_cur);
		targ_rec.push_back(targ_cur);
		//clear
		temp_ful="";
		targ_ful="";
		count++;
	}
	//return
	return count;
}

//----- proc seqname ----//
int WS_Proc_SeqName(string &infile,string &seq_name)
{
	ifstream fin;
	string buf,temp;
	//load
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!!\n",infile.c_str());
		exit(-1);
	}
	int len;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<13)continue;
		temp=buf.substr(0,13);
		if(temp=="Query Name = ")
		{
			seq_name=buf.substr(13,len-13);
			return 1;
		}
	}
	return 0;
}

//----- proc sequence ---//
int WS_Proc_Sequence(string &infile,string &seq_fasta)
{
	ifstream fin;
	string buf,temp;
	//load
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!!\n",infile.c_str());
		exit(-1);
	}
	int len;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<17)continue;
		temp=buf.substr(0,17);
		if(temp=="Query Sequence = ")
		{
			seq_fasta=buf.substr(17,len-17);
			return 1;
		}
	}
	return 0;
}


//----- proc list ----//
int WS_Proc_List(string &infile,int topk,
	vector <string> &temp_rec, vector <string> &targ_rec)
{
	int retv;
	//proc seqname
	string query_name;
	getBaseName(infile,query_name,'/','.');
	//proc seqfasta
	string seq_fasta;
	retv=WS_Proc_Sequence(infile,seq_fasta);
	if(retv==0)
	{
		fprintf(stderr,"file %s contains no Query Sequence !!\n",infile.c_str());
		exit(-1);
	}
	//output sequence
	string sequence_file=query_name+".seqfasta";
	FILE *fp=fopen(sequence_file.c_str(),"wb");
	fprintf(fp,">%s\n",query_name.c_str());
	fprintf(fp,"%s\n",seq_fasta.c_str());
	fclose(fp);
	//generate rank_tmp_ali
	char ws_command[30000];
	sprintf(ws_command,"awk '{ if(NF==6){ if($1==\"T\" || $1==\"S\") {print $1\" \"$2\" \"$4}  } }' %s > %s.rank_tmp_ali ",
		infile.c_str(),query_name.c_str());
	system(ws_command);
	//pre-proc alignment
	string proc_file=query_name+".rank_tmp_ali";
	retv=WS_Proc_Alignment(proc_file,topk,query_name,temp_rec,targ_rec);
	if(retv==0)
	{
		fprintf(stderr,"file %s contains no alignments !!\n",infile.c_str());
		exit(-1);
	}
	//post-proc alignment
	for(int i=0;i<retv;i++)
	{
		string fasta_file=temp_rec[i]+"-"+query_name+".fasta";
		string out_file=temp_rec[i]+"-"+query_name+".fasta_out";
		WS_Head_Tail_Mapping(fasta_file,sequence_file,out_file,1);
	}
	//delete rank_tmp_ali
	sprintf(ws_command,"rm -f %s",proc_file.c_str());
	system(ws_command);
	sprintf(ws_command,"rm -f %s",sequence_file.c_str());
	system(ws_command);
	//return 
	return retv;
}


//===================== build model =====================//
//----- multi build model -----//
void multi_build_model(string &out_nam,vector <string> &temp_list, string &targ_nam_out,
	string &fasta_root,string &pdb_root,string &fasta_suffix,string &mod_bin)
{
	int i;
	int size=(int)temp_list.size();
	int mod_num=1;
	#pragma omp parallel for schedule(dynamic)
	for(i=0;i<size;i++)
	{
		string fasta_file=fasta_root+"/"+temp_list[i]+"-"+targ_nam_out+"."+fasta_suffix;
		int retv=build3Dmodel(fasta_file,targ_nam_out,pdb_root,mod_bin,mod_num);
		char ws_command[30000];
		sprintf(ws_command,"mv %s_%s.B99990001.pdb %s-RPX-m%d-%s.pdb",
			temp_list[i].c_str(),targ_nam_out.c_str(),out_nam.c_str(),i+1,temp_list[i].c_str());
		system(ws_command);
	}
}


//================== usage =============//
void Usage(char *arg)
{
	printf("Version: 1.02 \n");
	printf("Usage: \n");
	printf("===================================================================\n");
//	printf("%s -i rank_file [-k TopK] [-d pdb_root] [-m mod_bin] \n",arg);
	printf("%s -i rank_file [-k TopK] \n",arg);
	printf("------------------------------------------------------------------ \n");
	printf("  -i rank_file: the rank file generated by CNFsearch. \n");
	printf(" [-k TopK]:     build models for templates among TopK in rank file.\n");
	printf("                (default=1) \n");
//	printf(" [-d pdb_root]: the folder containing the PDB file of the template.\n");
//	printf("                (default = databases/pdb_BC100/ ) \n");
//	printf(" [-m mod_bin]:  the MODELLER executable file, \n");
//	printf("                (default = ~/bin/modeller9v8/bin/mod9v8 ) \n");
//	printf(" [-n mod_num]:  the number of 3D models to be generated \n");
//	printf("                from the alignment (default=1) \n");
	printf("===================================================================\n");
}


//---- parameter editor ----//
static option long_options[] =
{
	{"input",  required_argument, NULL, 'i'},
	{"ktop",   no_argument, NULL, 'k'},
	{"data",   no_argument, NULL, 'd'},
	{"modbin",  no_argument, NULL, 'm'},
//	{"number",  no_argument, NULL, 'n'},
	{0, 0, 0, 0}
};

//--------- main ---------//
int main(int argc,char **argv)
{
	//----- multiple process build model -----//
	{
		if(argc==1)
		{
			Usage(argv[0]);
			exit(-1);
		}
		string rank_file="";
		int TopK=1;
		string pdb_root="databases/pdb_BC100/";
		string mod_bin="~/bin/modeller9v8/bin/mod9v8 ";

		//--- get parameter --//
		char c = 0;
		int option_index=0;
		while ((c = getopt_long(argc, argv, "i:k:d:m:",long_options,&option_index)) != EOF)
		{
			switch (c)
			{
				case 'i':
					rank_file = optarg;
					break;
				case 'k':
					TopK = atoi(optarg);
					break;
				case 'd':
					pdb_root = optarg;
					break;
				case 'm':
					mod_bin = optarg;
					break;
//				case 'n':
//					mod_num = atoi(optarg);
//					break;
				default:
					Usage(argv[0]);
					exit(-1);
			}
		}

		//--- check input ---//
		if(rank_file=="" || TopK<=0)
		{
			Usage(argv[0]);
			exit(-1);
		}
		string out_nam;
		getBaseName(rank_file,out_nam,'/','.');

		//----- list process ----//
		vector <string> temp_rec;
		vector <string> targ_rec;
		int retv=WS_Proc_List(rank_file,TopK,temp_rec,targ_rec);

		//----- build model -----//
//		string fasta_root=".";
//		string fasta_suffix="fasta_out";
//		multi_build_model(out_nam,temp_rec,out_nam,fasta_root,pdb_root,fasta_suffix,mod_bin);

		//======== delete exit =====//
/*
		char ws_command[30000];
		for(int i=0;i<retv;i++)
		{
			sprintf(ws_command,"rm -f %s-%s.fasta",
				temp_rec[i].c_str(),targ_rec[i].c_str());
			system(ws_command);
			sprintf(ws_command,"rm -f %s-%s.fasta_out",
				temp_rec[i].c_str(),targ_rec[i].c_str());
			system(ws_command);
		}
*/
		exit(0);
	}
}
