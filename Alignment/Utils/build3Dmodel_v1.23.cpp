#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <getopt.h>
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

//================= mapping_series ===================//
//DynaProg related
vector<pair<int,int> > WWW_alignment;
vector <int> wali1;
vector <int> wali2;

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
int process_oriami_record(char *seqres,char *ami_,int *int_,char *ins_,char *tag_)
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
		if(len<54)continue;
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
void WS_PDB_Mapping(char *seqres,string &pdbfile,string &seqws,string &dssp)
{
	int maxlen=3000;
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
		exit(-1);
	}
	//alignment
	int k;
	for(k=0;k<pdb_len;k++)tag_[k]=' ';
	tag_[k]='\0';
	int IsInsert=process_oriami_record(seqres,ami_,int_,ins_,tag_);
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
	char wstemp[6000];
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
}

//MAP
/*
void WS_Read_MAP(string &mapfile,string &seqres,string &dssp)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(mapfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		printf("no map file! %s \n",mapfile.c_str());
		exit(-1);
	}
	//process
	if(!getline(fin,buf,'\n'))
	{
		printf("map file bad! %s \n",mapfile.c_str());
		exit(-1);
	}
	if(!getline(fin,buf,'\n'))
	{
		printf("map file bad! %s \n",mapfile.c_str());
		exit(-1);
	}
	seqres="";
	seqres=buf;
	if(!getline(fin,buf,'\n'))
	{
		printf("map file bad! %s \n",mapfile.c_str());
		exit(-1);
	}
	if(!getline(fin,buf,'\n'))
	{
		printf("map file bad! %s \n",mapfile.c_str());
		exit(-1);
	}
	dssp="";
	dssp=buf;
}
*/

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
	char seqws[6000];
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
	string &mod_bin,string &pdb_root,string &pir_root,int MOL_NUM=1,int KEEP_FILE=0)
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
	char wscommand[3000];
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
	if(KEEP_FILE==0)
	{
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
}




//================== usage =============//
void Usage(char *arg)
{
	printf("Version: 1.23 \n");
	printf("Usage: \n");
	printf("===================================================================\n");
	printf("%s -i fasta_file -q query_name \n",arg);
	printf("        [ -d pdb_root] [ -m mod_bin ] [ -n mod_num ] [ -k keep_file ] \n");
	printf("-------------------------------------------------------------------\n");
	printf(" -i fasta_file:  the input alignment file in fasta format. \n");
	printf(" -q query_name:  the name of the query protein. \n");
	printf(" [-d pdb_root]:  the folder containing the PDB file of the template.\n");
	printf("                 (default = databases/pdb_BC100/ ) \n");
	printf(" [-m mod_bin]:   the MODELLER executable file, \n");
	printf("                 (default = ~/bin/modeller9v8/bin/mod9v8 ) \n");
	printf(" [-n mod_num]:   the number of 3D models to be generated \n");
	printf("                 from the alignment (default=1) \n");
	printf(" [-k keep_file]: keep temporary files (default=0, not keep) \n");
	printf("===================================================================\n");
//	printf("-----------------------------------------------------------\n");
//	printf(" the output model should be temp_targ.B9999000X.pdb, where \n");
//	printf(" temp is template name, targ is query name and X is the \n");
//	printf(" X-th model \n");
}
//---- parameter editor ----//
static option long_options[] =
{
	{"input",  required_argument, NULL, 'i'},
	{"query",  required_argument, NULL, 'q'},
	{"data",    no_argument, NULL, 'd'},
	{"modbin",  no_argument, NULL, 'm'},
	{"number",  no_argument, NULL, 'n'},
	{"keep",    no_argument, NULL, 'k'},
	{0, 0, 0, 0}
};


//---------- main ----------//
int main(int argc,char **argv)
{
	//------ mod script ----------//
	{
		if(argc==1)
		{
			Usage(argv[0]);
			exit(-1);
		}
		string fasta_file="";
		string query_name="";
		string pdb_root="databases/pdb_BC100/";
//		string mod_bin="/home/";
//		char *username=getlogin();
//		mod_bin=mod_bin+username+"/bin/modeller9v8/bin/mod9v8 ";
//		string mod_bin="mod9v8 ";
		string mod_bin="~/bin/modeller9v8/bin/mod9v8 ";
		int mod_num=1;
		int keep_file=0;

		//--- get parameter --//
		char c = 0;
		int option_index=0;
		while ((c = getopt_long(argc, argv, "i:q:d:m:n:k:",long_options,&option_index)) != EOF)
		{
			switch (c)
			{
				case 'i':
					fasta_file = optarg;
					break;
				case 'q':
					query_name = optarg;
					break;
				case 'd':
					pdb_root = optarg;
					break;
				case 'm':
					mod_bin = optarg;
					break;
				case 'n':
					mod_num = atoi(optarg);
					break;
				case 'k':
					keep_file = atoi(optarg);
					break;
				default:
				Usage(argv[0]);
				exit(-1);
			}
		}
		//--- check input ---//
		if(fasta_file=="" || query_name=="")
		{
			Usage(argv[0]);
			exit(-1);
		}
		//read fasta
		string nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2;
		vector<pair<int, int> > alignment;
		int retv=ReadToFile_FASTA(fasta_file,alignment,nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2);
		//check fasta
		if(retv!=1)
		{
			fprintf(stderr,"fasta_file %s not found\n",fasta_file.c_str());
			exit(-1);
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
			exit(-1);
		}
		//start modeller
		string pir_root=".";
		WS_Pair_Mod_Single_Fasta_Simp(nam1,nam2,mod_bin,pdb_root,pir_root,mod_num,keep_file);
		//delete tmp files
		string command;
		command="rm -f "+tmpnam;
		system(command.c_str());
		//exit
		exit(0);
	}
}
