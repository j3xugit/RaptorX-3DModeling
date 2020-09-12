#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include "getopt.h"
#include "cdhit_seqdb.h"
#include "blast_pssm.h"
using namespace std;


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

//---- remove tail blank ----//
void Remove_Tail_Blank(string &in,string &out)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]!=' ')break;
	}
	out=in.substr(0,i+1);
}

//-------- read in MSA in PSI format (i.e., normal FASTA with upper case) ------------//
//[note]: we set the first sequence as the query sequence,
//        that is to say, all the following sequences should be longer than the first
int Multi_FASTA_Input_PSI(string &multi_fasta,vector <string> &nam_list,vector <string> &fasta_list)
{
	ifstream fin;
	string buf,temp;

	//-> 1. read
	fin.open(multi_fasta.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",multi_fasta.c_str());
		return -1;
	}

	//-> 2. init
	int firstlen=0;
	int first=1;
	int count=0;
	string name;
	string seq;
	nam_list.clear();
	fasta_list.clear();

	//-> 3. process
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		temp=buf.substr(0,31);
		Remove_Tail_Blank(temp,name);
		seq=buf.substr(33,buf.length()-33);
		nam_list.push_back(name);
		fasta_list.push_back(seq);
		count++;
		if(first==1)
		{
			firstlen=(int)seq.length();
			first=0;
		}
		else
		{
			int curlen=(int)seq.length();
			if(curlen!=firstlen)
			{
				fprintf(stderr,"length not equal at %s, [%d!=%d] \n",buf.c_str(),curlen,firstlen);
				return -1;
			}
		}
	}

	//-> 4. return
	return count;
}

//====================== Check the given MSA ===================//
//-> MSA_Check
int MSA_Check(
	vector <string> &MSA_in,                                 //-> input, MSA
	vector <string> &MSA_out,                                //-> output, valid MSA
	string &query_rec)                                       //-> output, query record
{
	int i,k;
	int totnum=(int)MSA_in.size();
	int totlen=(int)MSA_in[0].length();
	MSA_out.clear();

	//-> 1. check query
	query_rec=MSA_in[0];
	for(k=0;k<totlen;k++)
	{
		int kResidue = BLAST_AA_To_ARND_Mapping(MSA_in[0][k]); //-> must return 0->19
		if(kResidue==-1)
		{
			fprintf(stderr,"ERROR: bad residue %c at query seq pos %d \r",MSA_in[0][k],k);
			MSA_in[0][k]='A';
//			return -1;
		}
	}
	MSA_out.push_back(MSA_in[0]);

	//-> 2. check other sequences
	int count=1;
	for(i=1;i<totnum;i++)
	{
		//check length
		if((int)MSA_in[i].length() != totlen) //length must equal !
		{
			fprintf(stderr,"ERROR: length not equal at %d sequence, [%d != %d]\n",
				i,(int)MSA_in[i].length(), totlen);
			return -1;
		}
		//check residue
		int valid=0;
		for(k=0;k<totlen;k++)
		{
			char pos=MSA_in[i][k];
			if(pos=='-')continue;  //skip gap
			if(pos>='a' && pos<='z') //turn lower case to gap
			{
				MSA_in[i][k]='-';
				continue;
			}
			int kResidue = BLAST_AA_To_ARND_Mapping(pos); //-> must return 0->19
			if(kResidue==-1)  //turn bad residue to gap
			{
				MSA_in[i][k]='-';
			}
			else valid=1;
		}
		//add this sequence
		if(valid==1)
		{
			MSA_out.push_back(MSA_in[i]);
			count++;
		}
	}

	//-> 3. return
	return count;
}

//============ MSA_To_Fragments ===========//
//given an input MSA, extracts all ungapped fragments, 
// and record their origin id and range.
int MSA_To_Fragments(
	vector <string> &MSA_in,                            //-> input, MSA
	int Minimal_UnGap_Len,                              //-> parameter, default should be 11
	vector <string> &out_nam,                           //-> output, OriginID/Range
	vector <string> &out_seq)                           //-> output, ungapped fragments
{
	int i,k;
	int totnum=(int)MSA_in.size();
	int totlen=(int)MSA_in[0].length();
	out_nam.clear();
	out_seq.clear();

	//process
	char command[30000];
	string nam_rec;
	string seq_rec;
	int count=0;
	for(i=0;i<totnum;i++)
	{
		int len=0;
		int start=-1;
		for(k=0;k<totlen;k++)
		{
			if(MSA_in[i][k]=='-')
			{
				//fragment process
				if(len>Minimal_UnGap_Len)
				{
					sprintf(command,"%d:%d-%d",i,start,start+len-1);
					nam_rec=command;
					out_nam.push_back(nam_rec);
					seq_rec=MSA_in[i].substr(start,len);
					out_seq.push_back(seq_rec);
					count++;
				}
				len=0;
				start=-1;
			}
			else
			{
				if(start==-1)start=k;
				len++;
			}
		}
		//terminal process
		if(len>Minimal_UnGap_Len)
		{
			sprintf(command,"%d:%d-%d",i,start,start+len-1);
			nam_rec=command;
			out_nam.push_back(nam_rec);
			seq_rec=MSA_in[i].substr(start,len);
			out_seq.push_back(seq_rec);
			count++;
		}
	}
	//return
	return count;
}



//============= CD-HIT process ===============//
int CD_HIT_StrTrans(string &input,string &output)
{
	int len=(int)input.length();
	if(len<=4)
	{
		fprintf(stderr,"input string %s bad ! \n",input.c_str());
		return -1;
	}
	if(input[0]!='>')
	{
		fprintf(stderr,"input string %s bad ! \n",input.c_str());
		return -1;
	}
	output=input.substr(1,len-4);
	return 1; //success
}

//[note]: sim_thres should be set to 0.95
//        len_thres should be set to 0.95
//[output]: all clusters with more than one members
int CD_HIT_Process_New(	
	vector <string> &out_nam,                           //-> input, OriginID/Range
	vector <string> &out_seq,                           //-> input, ungapped fragments
	double sim_thres,double len_thres,
	vector <vector <string> > & output, vector <int> & center )
{
	//-> 0. run CD-HIT
//	int sys_retv;
//	char command[30000];
//	sprintf(command,"./cd-hit -i %s -o %s -c %lf -s %lf -M 0 -T 0 -d 0 1> ws1 2> ws2",
//		infile.c_str(),outname.c_str(),sim_thres,len_thres);
//	sys_retv=system(command);

	//-> 1. initial CDHIT
	Options options;
	SequenceDB seq_db;
	options.cluster_thd=sim_thres;   // -c sim_thres
	options.diff_cutoff=len_thres;   // -s len_thres
	options.max_memory=0;            // -M 0
	options.des_len=0;               // -d 0
#ifndef NO_OPENMP
	int cpu = omp_get_num_procs();   // -T 0
	options.threads = cpu;
#endif
	options.Validate();
	InitNAA( MAX_UAA );
	options.NAAN = NAAN_array[options.NAA];
	seq_db.NAAN = NAAN_array[options.NAA];

	//-> 2. run CDHIT to get all clusters with more than one members
	seq_db.Read_From_String( options, out_nam, out_seq );
	seq_db.SortDivide( options );
	seq_db.DoClustering( options );
	seq_db.WriteExtra1D_To_String( options, output,center);

	//return
	CDHIT_CleanUpTempFiles();
	return (int)center.size();
}


//============= CD-HIT Purge ===============//
int CD_HIT_StrProc(string &in_str,int &start,int &end,
	int totnum,int totlen)
{
	int i;
	int cur;
	int seq_id;
	string temp;
	int len=(int)in_str.length();
	//get seq_id
	for(i=0;i<len;i++)
	{
		if(in_str[i]==':')break;
	}
	temp=in_str.substr(0,i);
	seq_id=atoi(temp.c_str());
	if(seq_id<0 || seq_id>=totnum)return -1;
	//get start
	cur=i+1;
	for(i=cur;i<len;i++)
	{
		if(in_str[i]=='-')break;
	}
	temp=in_str.substr(cur,i-cur);
	start=atoi(temp.c_str());
	if(start<0 || start>=totlen)return -1;
	//get end
	cur=i+1;
	temp=in_str.substr(cur,len-cur);
	end=atoi(temp.c_str());
	if(end<0 || end>=totlen)return -1;
	//return
	return seq_id;
}

//[note]: thres should be set to 0.95
int CD_HIT_CalcOverlap(
	int cent_seqid,int cent_start,int cent_end,
	int memb_seqid,int memb_start,int memb_end,
	int &start, int &end,double thres)
{
	//init
	int cent_len=cent_end-cent_start+1;
	int memb_len=memb_end-memb_start+1;
	int minlen=cent_len<memb_len?cent_len:memb_len;
	start=memb_start>cent_start?memb_start:cent_start;
	end=memb_end<cent_end?memb_end:cent_end;
	int len=end-start+1;
	//check
	if(1.0*len/minlen>thres)return 1;  //succeed
	else return -1;  //failed
}

//[note]: thres should be set to 0.95
void CD_HIT_Purge(
	vector <string> &MSA_in,                    //-> input/output, MSA
	vector <vector <string> > & purge_string,   //-> input, purge sequence
	vector <int> & purge_center,                //-> input, purge center
	double overlap_thres)                       //-> input, overlap threshold
{
	//init
	int totnum=(int)MSA_in.size();
	int totlen=(int)MSA_in[0].length();
	//process
	int i,k,l;
	int retv;
	int size=(int)purge_center.size();
	int center;
	int member;
	string cent_str;
	string memb_str;
	int cent_seqid,cent_start,cent_end;
	int memb_seqid,memb_start,memb_end;
	int start,end;
	for(i=0;i<size;i++)
	{
		//get center string
		center=purge_center[i];
		cent_str=purge_string[i][center];
		cent_seqid=CD_HIT_StrProc(cent_str,cent_start,cent_end,totnum,totlen);
		if(cent_seqid==-1)
		{
			fprintf(stderr,"cent_str %s overalp ! \n",cent_str.c_str());
			continue;
		}
		//process member
		member=(int)purge_string[i].size();
		for(k=0;k<member;k++)
		{
			if(k==center)continue;
			memb_str=purge_string[i][k];
			memb_seqid=CD_HIT_StrProc(memb_str,memb_start,memb_end,totnum,totlen);
			if(cent_seqid==-1)
			{
				fprintf(stderr,"memb_str %s overalp ! \n",memb_str.c_str());
				continue;
			}
			//check overlap
			retv=CD_HIT_CalcOverlap(cent_seqid,cent_start,cent_end,
				memb_seqid,memb_start,memb_end,start,end,overlap_thres);
			if(retv==1 && memb_seqid!=0) //purge this fragment
			{
				for(l=start;l<=end;l++)MSA_in[memb_seqid][l]='-';
			}
		}
	}
}

//-> MSA_Purge_Check
int MSA_Purge_Check(
	vector <string> &MSA_in,
	vector <string> &MSA_out)
{
	MSA_out.clear();
	//init
	int totnum=(int)MSA_in.size();
	if(totnum<=0)return totnum;
	int totlen=(int)MSA_in[0].length();
	//process
	int i,k;
	int valid;
	int count=1;
	MSA_out.push_back(MSA_in[0]); //query sequence
	for(i=1;i<totnum;i++)
	{
		valid=0;
		for(k=0;k<totlen;k++)
		{
			if(MSA_in[i][k]!='-')
			{
				valid=1;
				break;
			}
		}
		if(valid==1)
		{
			MSA_out.push_back(MSA_in[i]);
			count++;
		}
	}
	return count;
}



//=============== Purge Pairwise Alignment =============//
void MSA_Purge_Pair_Alignments(
	vector <string> &MSA_in, vector <int> &MSA_valid,
	int seq_index1,int seq_index2,
	double max_percent_identity)
{
	//-> 1. init check 
	if ( seq_index1 == seq_index2 ) return;
	if( MSA_valid[seq_index1]==0 || MSA_valid[seq_index2]==0 )return;

	//-> 2. init
	int k;
	int p;
	int query_len=(int)MSA_in[0].length();
	int len=0;
	int match=0;
	int start=-1;

	//-> 3. process
	for (p = 0; p < query_len; p++) 
	{  	
		// if neither position is aligned
		char pos1=MSA_in[seq_index1][p];
		char pos2=MSA_in[seq_index2][p];
		int kPos1Aligned = pos1=='-'?0:1;
		if(seq_index1==0)kPos1Aligned=0;
		int kPos2Aligned = pos2=='-'?0:1;

		// judge aligned or not
		if( !kPos1Aligned && !kPos2Aligned )
		{
			if(len>0)
			{
				if(1.0*match/len>=max_percent_identity)
				{
					//purge alignment on seq_index2
					for(k=start;k<start+len;k++)MSA_in[seq_index2][k]='-';
				}
			}
			len=0;
			match=0;
			start=-1;
		}
		else
		{
			if(start==-1)start=p;
			if(pos1==pos2)match++;
			len++;
		}
	}

	//-> 4. terminal process
	if(len>0)
	{
		if(1.0*match/len>=max_percent_identity)
		{
			//purge alignment on seq_index2
			for(k=start;k<start+len;k++)MSA_in[seq_index2][k]='-';
		}
	}

	//-> 5. complete purge
	int isaligned=0;
	for(k=0;k<query_len;k++)
	{
		if(MSA_in[seq_index2][k]!='-')
		{
			isaligned=1;
			break;
		}
	}
	if(isaligned==0)
	{
		MSA_valid[seq_index2]=0;
	}
}

//-> MSA_Purge (just self-hits)
int MSA_Purge(
	vector <string> &MSA_in,
	vector <string> &MSA_out,
	int CDHIT_or_NOT)
{
	//-> 1. init
	int i;
	int totnum=(int)MSA_in.size();
	vector <int> MSA_valid(totnum,1); //default is OK

	//-> 2. purge self hits
	double kPSIIdentical=1.0;
	for(i=1;i<totnum;i++)
	{
		MSA_Purge_Pair_Alignments(MSA_in,MSA_valid,0,i,kPSIIdentical);
	}

	//-> 3. purge pairwise
	if(CDHIT_or_NOT==0) //no cdhit applied !!
	{
		int j;
		double kPSINearIdentical=0.94;
		for (i = 1; i < totnum; i++) 
		{ 
			for (j = 1; (i + j) < totnum; j++) 
			{
				// N.B.: The order of comparison of sequence pairs is deliberate,
				// tests on real data indicated that this approach allowed more
				// sequences to be purged 
				MSA_Purge_Pair_Alignments(MSA_in,MSA_valid, j, (i + j), kPSINearIdentical);
			}
		}
	}

	//-> 4. output
	MSA_out.clear();
	MSA_out.push_back(MSA_in[0]);
	int count=1;
	for(i=1;i<totnum;i++)
	{
		//add this sequence
		if(MSA_valid[i]==1)
		{
			MSA_out.push_back(MSA_in[i]);
			count++;
		}
	}

	//-> 5. return
	return count;
}



//=========== BLAST PSSM File Output =============//
/* an example is as follows,

Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts
            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    1 E    -2   1   0   3  -4   2   5  -2  -1  -4  -3   2  -3  -4  -2  -1  -1  -4  -3  -3    0   7   0  18   0   4  54   0   0   0   0  16   0   0   0   0   0   0   0   0  0.65 0.10
    2 N    -2  -2   2   3  -3   0   3  -2  -1  -3  -2  -1  -2  -3  -2   0   3  -4  -3  -2    0   0  13  21   0   0  29   0   0   0   5   0   0   0   0   0  31   0   0   0  0.48 0.12
...

                      K         Lambda
PSI Ungapped         0.1319     0.3166
PSI Gapped           0.0404     0.2670

*/
void BLAST_PSSM_File_Output(FILE *fp,
	string &query,                                      //-> input, query string
	int *pssm_in, double *match_weights,                //-> input, pssm and match_weights
	double *freq_entro, double *column_weights,         //-> input, entropy and column_weights
	double gap_lambda,double gap_K,double gap_H,        //-> input, gap KARLIN parameter
	double ungap_lambda,double ungap_K,double ungap_H)  //-> input, ungap KARLIN parameter
{
	//output header
	fprintf(fp,"\nLast position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n");
	fprintf(fp,"            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V\n");
	//output others
	int i,k;
	int length=(int)query.length();
	for(i=0;i<length;i++)
	{
		//head
		fprintf(fp,"%5d %c  ",i+1,query[i]);
		//pssm
		for(k=0;k<20;k++)fprintf(fp," %3d",pssm_in[i*20+k]);
		fprintf(fp," ");
		//match_weights
		for(k=0;k<20;k++)
		{
			double value=100;
			value*=match_weights[i*20+k];
			int value_int = (int)(value + (value > 0. ? 0.5 : -0.5));
			fprintf(fp," %3d",value_int);
		}
		fprintf(fp," ");
		//entropy
		fprintf(fp," %4.2f ",freq_entro[i]);
		//column_weights
		fprintf(fp,"%4.2f\n",column_weights[i]);
	}
	//output tailor
	fprintf(fp,"\n                      K         Lambda\n");
	fprintf(fp,"PSI Ungapped         %6.4f     %6.4f\n",ungap_K,ungap_lambda);
	fprintf(fp,"PSI Gapped           %6.4f     %6.4f\n",gap_K,gap_lambda);
}



//=========== BLAST Checkpoint File Output in ASN.1 format ===========//
/* an example for query sequence is as follows,
        seq-data iupacaa "TITQDTPINQIFTDTALAEKMKTVLGKTNVTDTVSQTDLDQVTTLQADRLGI
KSIDGVEYLNNLTQINFSNNQLTDITPLKNLTKLVDILMNNNQIADITPLANLTNLTGLTLFNNQITDIDPLKNLTNL
NRLELSSNTISDISALSGLTSLQQLSFSSNQVTDLKPLANLTTLERLDISSNKVSDISVLAKLTNLESLIATNNQISD
ITPLGILTNLDELSLNGNQLKDIGTLASLTNLTDLDLANNQISNLAPLSGLTKLTELKLGANQISNISPLAGLTALTN
LELNENQLEDISPISNLKNLTYLTLYFNNISDISPVSSLTKLQRLFFYNNKVSDVSSLANLTNINWLSAGHNQISDLT
PLANLTRITQLGLNDQAWTNAPVNYKANVSIPNTVKNVTGALIAPATISDGGSYTEPDITWNLPSYTNEVSYTFSQPV
TIGKGTTTFSGTVTQPLK"
*/
void BLAST_Checkpoint_File_StrProc(string &in_str,vector <string> &out_str)
{
	//init
	out_str.clear();
	//init proc
	int i;
	int start=52;
	int cutlen=78;
	int length=(int)in_str.length();
	if(length<start)
	{
		out_str.push_back(in_str);
		return;
	}
	//real proc
	string header=in_str.substr(0,start);
	out_str.push_back(header);
	int divide=(length-start)/cutlen;
	for(i=0;i<divide;i++)
	{
		string curstr=in_str.substr(start+i*cutlen,cutlen);
		out_str.push_back(curstr);
	}
	string tailer=in_str.substr(start+divide*cutlen,length-start+divide*cutlen);
	if(tailer!="")out_str.push_back(tailer);
}

// double (or float) trasfer to ASN.1 string 
// for example, { 470521826525912, 10, -16 }, -> 1.2345678123456780e-01
// by default, prec = 14
void Double_To_ASN1(double input, string &output,int prec=14)
{
	output="";
	//get temp string
	stringstream ss;
	ss.precision(prec);
	ss << scientific << input;
	string temp=ss.str();
	//get final string
	string head;
	string midd;
	string tail;
	int i;
	int len=(int)temp.length();
	for(i=0;i<len;i++)if(temp[i]=='e')break;
	head=temp.substr(0,1);
	midd=temp.substr(2,i-2);
	tail=temp.substr(i+1,len-i-1);
	//process
	string finhead=head+midd;
	temp=tail.substr(1,(int)tail.length()-1);
	int sign=1;
	if(tail[0]=='-')sign=-1;
	int digit=sign*atoi(temp.c_str())-prec;
	stringstream ww;
	ww << digit;
	string findigit=ww.str();
	//final
	output=finhead + ", 10, "+findigit;
}

//---- output ASN.1 format BLAST checkpoint file //-> for new BLAST+ (say, psiblast)
void BLAST_Checkpoint_File_Output(FILE *fp,
	string &query, string &query_nam, string &matnam,   //-> input, query string
	double *freq_in,int *pssm_in,                       //-> input, frequency and pssm
	double gap_lambda,double gap_K,double gap_H,        //-> input, gap KARLIN parameter
	double ungap_lambda,double ungap_K,double ungap_H)  //-> input, ungap KARLIN parameter
{
	//init setup
	int digit=28;
	int count;
	//get length
	int i,k;
	int length=(int)query.length();
	//output header
	fprintf(fp,"PssmWithParameters ::= {\n");
	fprintf(fp,"  pssm {\n");
	fprintf(fp,"    isProtein TRUE,\n");
	fprintf(fp,"    numRows %d,\n",digit);
	fprintf(fp,"    numColumns %d,\n",length);
	fprintf(fp,"    byRow FALSE,\n");
	fprintf(fp,"    query seq {\n");
	fprintf(fp,"      id {\n");
	fprintf(fp,"        local str \"%s\"\n",query_nam.c_str());
	fprintf(fp,"      },\n");
	fprintf(fp,"      inst {\n");
	fprintf(fp,"        repr raw,\n");
	fprintf(fp,"        mol aa,\n");
	fprintf(fp,"        length %d,\n",length);
	//output query string
	{
		vector <string> out_str;
		BLAST_Checkpoint_File_StrProc(query,out_str);
		int size=out_str.size();
		if(size>1)
		{
			fprintf(fp,"        seq-data iupacaa \"%s\n",out_str[0].c_str());
			for(i=1;i<size-1;i++)fprintf(fp,"%s\n",out_str[i].c_str());
			fprintf(fp,"%s\"\n",out_str[size-1].c_str());
		}
		else
		{
			fprintf(fp,"        seq-data iupacaa \"%s\"\n",out_str[0].c_str());
		}
		fprintf(fp,"      }\n");
		fprintf(fp,"    },\n");
	}
	//output frequency
	fprintf(fp,"    intermediateData {\n");
	fprintf(fp,"      freqRatios {\n");
	count=0;
	for(i=0;i<length;i++)
	{
		for(k=0;k<digit;k++) //28
		{
			int pos=BLAST_Digit_To_ARND_Mapping(k);
			if(pos<0)
			{
				if(count<length*digit-1)fprintf(fp,"        { 0, 10, 0 },\n");
				else fprintf(fp,"        { 0, 10, 0 }\n"); 
			}
			else
			{
				string output;
				double input=freq_in[i*20+pos];
				Double_To_ASN1(input, output);
				if(count<length*digit-1)fprintf(fp,"        { %s },\n",output.c_str());
				else fprintf(fp,"        { %s }\n",output.c_str());
			}
			count++;
		}
	}
	fprintf(fp,"      }\n");
	fprintf(fp,"    },\n");
	//output pssm
	fprintf(fp,"    finalData {\n");
	fprintf(fp,"      scores {\n");
	count=0;
	for(i=0;i<length;i++)
	{
		for(k=0;k<digit;k++) //28
		{
			if(k==21)
			{
				fprintf(fp,"        -1,\n");
			}
			else if(k==25)
			{
				fprintf(fp,"        -4,\n");
			}
			else
			{
				int pos=BLAST_Digit_To_ARND_Mapping(k);
				if(pos<0)
				{
					if(count<length*digit-1)fprintf(fp,"        -32768,\n");
					else fprintf(fp,"        -32768\n");
				}
				else
				{
					int pssm_out=pssm_in[i*20+pos];
					if(count<length*digit-1)fprintf(fp,"        %d,\n",pssm_out);
					else fprintf(fp,"        %d\n",pssm_out);
				}
			}
			count++;
		}
	}
	fprintf(fp,"      },\n");
	//output KARLIN parameter
	{
		string output;
		//-> gap 
		Double_To_ASN1(gap_lambda, output);
		fprintf(fp,"      lambda { %s },\n",output.c_str());
		Double_To_ASN1(gap_K, output);
		fprintf(fp,"      kappa { %s },\n",output.c_str());
		Double_To_ASN1(gap_H, output);
		fprintf(fp,"      h { %s },\n",output.c_str());
		//-> ungap
		Double_To_ASN1(ungap_lambda, output);
		fprintf(fp,"      lambdaUngapped { %s },\n",output.c_str());
		Double_To_ASN1(ungap_K, output);
		fprintf(fp,"      kappaUngapped { %s },\n",output.c_str());
		Double_To_ASN1(ungap_H, output);
		fprintf(fp,"      hUngapped { %s }\n",output.c_str());
	}
	//final output
	fprintf(fp,"    }\n");
	fprintf(fp,"  },\n");
	fprintf(fp,"  params {\n");
	fprintf(fp,"    pseudocount 0,\n");
	fprintf(fp,"    rpsdbparams {\n");
	fprintf(fp,"      matrixName \"%s\"\n",matnam.c_str());
	fprintf(fp,"    }\n");
	fprintf(fp,"  }\n");
	fprintf(fp,"}\n");
}
//---- output ASN.1 format BLAST checkpoint file //-> for old BLAST (say, blastpgp)
void BLAST_Checkpoint_File_Output(FILE *fp,
	string &query, string &query_nam, double *freq_in)   //-> input, query string
{
	//init setup
	int digit=28;
	int count;
	//get length
	int i,k;
	int length=(int)query.length();
	//output header
	fprintf(fp,"PssmWithParameters ::= {\n");
	fprintf(fp,"  pssm {\n");
	fprintf(fp,"    isProtein TRUE ,\n");
	fprintf(fp,"    numRows %d ,\n",digit);
	fprintf(fp,"    numColumns %d ,\n",length);
	fprintf(fp,"    byRow FALSE ,\n");
	fprintf(fp,"    query\n");
	fprintf(fp,"      seq {\n");
	fprintf(fp,"        id {\n");
	fprintf(fp,"          local\n");
	fprintf(fp,"            str \"%s\" } ,\n",query_nam.c_str());
	fprintf(fp,"        inst {\n");
	fprintf(fp,"          repr raw ,\n");
	fprintf(fp,"          mol aa ,\n");
	fprintf(fp,"          length %d ,\n",length);
	fprintf(fp,"          seq-data\n");
	//output query string
	{
		vector <string> out_str;
		BLAST_Checkpoint_File_StrProc(query,out_str);
		int size=out_str.size();
		if(size>1)
		{
			fprintf(fp,"            ncbieaa \"%s\n",out_str[0].c_str());
			for(i=1;i<size-1;i++)fprintf(fp,"%s\n",out_str[i].c_str());
			fprintf(fp,"%s\" } } ,\n",out_str[size-1].c_str());
		}
		else
		{
			fprintf(fp,"            ncbieaa \"%s\" } } ,\n",out_str[0].c_str());
		}
	}
	//output frequency
	fprintf(fp,"    intermediateData {\n");
	fprintf(fp,"      freqRatios {\n");
	count=0;
	for(i=0;i<length;i++)
	{
		for(k=0;k<digit;k++) //28
		{
			int pos=BLAST_Digit_To_ARND_Mapping(k);
			if(pos<0)
			{
				if(count<length*digit-1)fprintf(fp,"        { 0, 10, 0 } ,\n");
				else fprintf(fp,"        { 0, 10, 0 } } } } ,\n"); 
			}
			else
			{
				string output;
				double input=freq_in[i*20+pos];
				Double_To_ASN1(input, output);
				if(count<length*digit-1)fprintf(fp,"        { %s } ,\n",output.c_str());
				else fprintf(fp,"        { %s } } } } ,\n",output.c_str());
			}
			count++;
		}
	}
	//final output
	fprintf(fp,"  params {\n");
	fprintf(fp,"    pseudocount 0 ,\n");
	fprintf(fp,"    rpsdbparams {\n");
	fprintf(fp,"      matrixName \"BLOSUM62\" ,\n");
	fprintf(fp,"      gapOpen 11 ,\n");
	fprintf(fp,"      gapExtend 1 } } }\n");
}

//=========== BLAST MTX Format File Output ===========//
void BLAST_MTX_Output(FILE *fp,string &query,           //-> input, query length
	double *freq_in,double *std_prob,double ideal_lambda, //-> input, frequency and std_prob
	double origap_lambda,double origap_K,double origap_H, //-> input, origap KARLIN parameter
	double gap_lambda,double gap_K,double gap_H,          //-> input, gap KARLIN parameter
	double ungap_lambda,double ungap_K,double ungap_H)    //-> input, ungap KARLIN parameter
{
	int i,k;
	int digit=28;
	int len=(int)query.length();
	//sequence output
	fprintf(fp,"%d\n",len);
	fprintf(fp,"%s\n",query.c_str());
	//KARLIN output
	//-> origap
	fprintf(fp,"%e\n",origap_lambda/100.0);
	fprintf(fp,"%e\n",origap_K);
	fprintf(fp,"%e\n",log(origap_K));
	fprintf(fp,"%e\n",origap_H);
	//-> gap
	fprintf(fp,"%e\n",gap_lambda/100.0);
	fprintf(fp,"%e\n",gap_K);
	fprintf(fp,"%e\n",log(gap_K));
	fprintf(fp,"%e\n",gap_H);
	//-> ungap
	fprintf(fp,"%e\n",ungap_lambda/100.0);
	fprintf(fp,"%e\n",ungap_K);
	fprintf(fp,"%e\n",log(ungap_K));
	fprintf(fp,"%e\n",ungap_H);
	//output MTX
	for(i=0;i<len;i++)
	{
		for(k=0;k<digit;k++) //28
		{
			if(k==21)
			{
				fprintf(fp,"-100  ");
			}
			else if(k==25)
			{
				fprintf(fp,"-396  ");
			}
			else
			{
				int pos=BLAST_Digit_To_ARND_Mapping(k);
				if(pos<0)
				{
					fprintf(fp,"-32768  ");
				}
				else
				{
					double qOverPEstimate = freq_in[i*20+pos] / std_prob[pos];
					int score = BLAST_Nint (100.0*log(qOverPEstimate)/ideal_lambda);
					fprintf(fp,"%d  ",score);
				}
			}
		}//end of digit
		fprintf(fp,"\n");
	}//end of len
}


//---------- usage ---------//
void Usage() 
{
	fprintf(stderr,"Version: 1.05 \n");
	fprintf(stderr,"MSA_To_PSSM -i psi_input [-l block_len] [-o pssm_out] [-m mtx_out]  \n");
	fprintf(stderr,"          [-t chk_out_new] [-T chk_out_old] [-c cut_num] \n\n");
	fprintf(stderr,"Usage : \n\n");
	fprintf(stderr,"-i psi_input :         Input MSA file in PSI format. \n\n");
	fprintf(stderr,"-l block_len :         Block length for blasgpgp operation. \n");
	fprintf(stderr,"                       (by default, block_len = 5) \n\n");
	fprintf(stderr,"-o pssm_out :          Output PSSM file in ARND style. \n\n");
	fprintf(stderr,"-m mtx_out :           Output blast MTX format file. \n\n");
	fprintf(stderr,"-t chk_out_new :       Output blast checkpoint file for BLAST+. \n\n");
	fprintf(stderr,"-T chk_out_old :       Output blast checkpoint file for blastpgp. \n\n");
	fprintf(stderr,"-c cut_num :           If MSA number is larger than cut_num, \n");
	fprintf(stderr,"                       then cdhit will be applied during purge. \n");
	fprintf(stderr,"                       (by default, cdhit won't be called.) \n");
	fprintf(stderr,"                       [if set, please use -c 20] \n\n");
}

//-------- main -------//
int main(int argc,char **argv)
{
	//------ MSA_To_PSSM -------//
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		string psi_input="";
		string pssm_output="";
		string chk_output="";
		string mtx_output="";
		int block_len=5;
		int CDHIT_or_NOT=0; //default: DONOT use cdhit to purge near identical sequence
		int pssm_key=0;
		int chk_key=0;
		int mtx_key=0;
		int cut_num=0;

		//command-line arguments process
		extern char* optarg;
		char c = 0;
		while ((c = getopt(argc, argv, "i:l:o:m:t:T:c:")) != EOF) {
			switch (c) {
			case 'i':
				psi_input = optarg;
				break;
			case 'l':
				block_len = atoi(optarg);
				break;
			case 'o':
				pssm_output = optarg;
				pssm_key=1;
				break;
			case 'm':
				mtx_output = optarg;
				mtx_key=1;
				break;
			case 't':
				chk_output = optarg;
				chk_key=1;
				break;
			case 'T':
				chk_output = optarg;
				chk_key=2;
				break;
			case 'c':
				cut_num = atoi(optarg);
				CDHIT_or_NOT = 1;
				break;
			default:
				Usage();
				exit(-1);
			}
		}

		//check arguments
		if(psi_input=="")
		{
			fprintf(stderr,"psi_input should be specified \n");
			exit(-1);
		}
		if( (pssm_key+chk_key+mtx_key)==0  )
		{
			fprintf(stderr,"at least one of the output should be specified \n");
			exit(-1);
		}

		// PSI input
		int retv;
		vector <string> nam_list;
		vector <string> fasta_list;
		retv=Multi_FASTA_Input_PSI(psi_input,nam_list,fasta_list);
		if(retv<=0)
		{
			fprintf(stderr,"Multi_FASTA_Input_PSI error with return code %d \n",retv);
			exit(-1);
		}
		// MSA check
		vector <string> MSA_init;
		string query_rec;
		retv=MSA_Check(fasta_list,MSA_init,query_rec);
		if(retv<=0)
		{
			fprintf(stderr,"MSA_Check error with return code %d \n",retv);
			exit(-1);
		}


//--- we add a checkpoint here, if input sequence is relatively small, then there is no need to call cd-hit -----//
if(retv<cut_num)CDHIT_or_NOT=0;
//--- we add a checkpoint here ---//over

		// MSA purge
		vector <string> MSA_fin_;
		retv=MSA_Purge(MSA_init,MSA_fin_,CDHIT_or_NOT);
		if(retv<=0)
		{
			fprintf(stderr,"MSA_Purge 1nd stage error with return code %d \n",retv);
			exit(-1);
		}

		//============ CDHIT purge ===========//
		vector <string> MSA_fin;
		if(CDHIT_or_NOT==1) //cdhit process
		{
			//-> 1. MSA fragment extract
			int Minimal_UnGap_Len=11;
			vector <string> Frag_nam;
			vector <string> Frag_seq;
			MSA_To_Fragments(MSA_fin_,Minimal_UnGap_Len,Frag_nam,Frag_seq);
			
			//-> 2. CDHIT Purge main process ============//
			double sim_threshold=0.93;
			double len_threshold=0.93;
			vector <vector <string> > purge_string;
			vector <int> purge_center;
			retv=CD_HIT_Process_New(Frag_nam,Frag_seq,sim_threshold,len_threshold,
				purge_string,purge_center);
			if(retv>0)CD_HIT_Purge(MSA_fin_,purge_string,purge_center,sim_threshold);

			//-> 3. Get final MSA
			retv=MSA_Purge_Check(MSA_fin_,MSA_fin);
			if(retv<=0)
			{
				fprintf(stderr,"MSA_Purge 2nd stage error with return code %d \n",retv);
				exit(-1);
			}
		}
		else                //no cdhit, then just copy MSA_fin_ to MSA_fin
		{
			MSA_fin.clear();
			for(int i=0;i<(int)MSA_fin_.size();i++)MSA_fin.push_back(MSA_fin_[i]);
		}

		// PSSM process
		int *PSSM_Out;
		double *Weight_Out;
		double *Freq_Out;
		double *Std_Prob;
		double *Freq_Entro;
		double *Colum_Weights_P;
		int query_len=(int)MSA_fin[0].length();
		string matnam;
		double idea_lambda;
		double origap_lambda,origap_K,origap_H;
		double gap_lambda,gap_K,gap_H;
		double ungap_lambda,ungap_K,ungap_H;
		{
			BLAST_PSSM blast_pssm;
			//init matrix
			string input_matrix="BLOSUM62";
			int gap_open=11;
			int gap_extend=1;
			blast_pssm.Init_BLAST_Stat(input_matrix,gap_open,gap_extend);    //init matrix related data
			blast_pssm.Init_Pseudo_Expno();              //init pseudo count
			blast_pssm.Create_DataStructure(query_len);  //init PSSM related data
			retv=blast_pssm.MSA_To_PSSM(MSA_fin,block_len);
			if(retv!=0)
			{
				fprintf(stderr,"blast_pssm.MSA_To_PSSM error with return code %d \n", retv);
				exit(-1);
			}
			//assign output
			matnam=input_matrix;
			PSSM_Out=new int[query_len*20];
			Weight_Out=new double[query_len*20];
			Freq_Out=new double[query_len*20];
			Std_Prob=new double[20];
			Freq_Entro=new double[query_len];
			Colum_Weights_P=new double[query_len];
			for(int i=0;i<query_len*20;i++)PSSM_Out[i]=blast_pssm.PSSM_Fin[i];
			for(int i=0;i<query_len*20;i++)Weight_Out[i]=blast_pssm.Match_Weights[i];
			for(int i=0;i<query_len*20;i++)Freq_Out[i]=blast_pssm.Freq_Ratio[i];
			for(int i=0;i<20;i++)Std_Prob[i]=blast_pssm.Freq_Std[i];
			for(int i=0;i<query_len;i++)Freq_Entro[i]=blast_pssm.Freq_Entro[i];
			for(int i=0;i<query_len;i++)Colum_Weights_P[i]=blast_pssm.Colum_Weights_P[i];
			//assign KARLIN parameter
			idea_lambda=blast_pssm.Ideal_Lambda;
			origap_lambda=blast_pssm.Gap_Lambda;
			origap_K=blast_pssm.Gap_K;
			origap_H=blast_pssm.Gap_H;
			gap_lambda=blast_pssm.Gap_Lambda;
			gap_K=blast_pssm.Cur_Gap_K;
			gap_H=blast_pssm.Gap_H;
			ungap_lambda=blast_pssm.Cur_Lambda;
			ungap_K=blast_pssm.Cur_K;
			ungap_H=blast_pssm.Cur_H;

/*
//testout
{
int count=0;
for(int i=0;i<query_len;i++)
{
	for(int k=0;k<20;k++)
	{
		printf("%6.4f ",Freq_Out[count]);
		count++;
	}
	printf("\n");
}
}
*/
		}

		// PSSM output in ARNDCQEGHILKMFPSTWYVZ code
		FILE *fp;
		MSA_fin[0]=query_rec;
		if(pssm_output!="")
		{
			fp=fopen(pssm_output.c_str(),"wb");
			string query=MSA_fin[0];
			BLAST_PSSM_File_Output(fp,query,PSSM_Out,Weight_Out,Freq_Entro,Colum_Weights_P,
				gap_lambda,gap_K,gap_H,ungap_lambda,ungap_K,ungap_H);
			fclose(fp);
		}

		// Checkpoint file output in BLAST ASN.1 format
		if(chk_output!="")
		{
			fp=fopen(chk_output.c_str(),"wb");
			string query=MSA_fin[0];
			string query_nam="QUERY";
			if(chk_key==1) //output new format for BLAST+
			{
				BLAST_Checkpoint_File_Output(fp,query,query_nam,matnam,Freq_Out,PSSM_Out,
					gap_lambda,gap_K,gap_H,ungap_lambda,ungap_K,ungap_H);
			}
			else           //output old format for BLAST+
			{
				BLAST_Checkpoint_File_Output(fp,query,query_nam,Freq_Out);
			}
			fclose(fp);
		}

		// BLAST MTX file by makemat
		if(mtx_output!="")
		{
			fp=fopen(mtx_output.c_str(),"wb");
			string query=MSA_fin[0];
			BLAST_MTX_Output(fp,query,Freq_Out,Std_Prob,idea_lambda,
				origap_lambda,origap_K,origap_H,gap_lambda,gap_K,gap_H,ungap_lambda,ungap_K,ungap_H);
			fclose(fp);
		}

		//========= exit =========//
		delete [] PSSM_Out;
		delete [] Weight_Out;
		delete [] Freq_Out;
		delete [] Std_Prob;
		delete [] Freq_Entro;
		delete [] Colum_Weights_P;
		exit(0);
	}
}
