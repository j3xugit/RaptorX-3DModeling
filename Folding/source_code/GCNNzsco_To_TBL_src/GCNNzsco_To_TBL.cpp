#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include "Fast_Sort.h"
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
int Read_FASTA_SEQRES(string &fasta_file,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fasta_file.c_str(), ios::in);
	if(fin.fail()!=0)
		{
		fprintf(stderr,"no such file! %s \n",fasta_file.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",fasta_file.c_str());
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

//-------- sequence cut to pieces -----------//
void Sequence_Cut_To_Pieces(string &input,vector <string> &output,int cut)
{
	int i;
	int len=(int)input.length();
	int start=0;
	output.clear();
	for(i=0;i<len;i++)
	{
		if(i>0 && i%cut==0)
		{
			string sub=input.substr(start,i-start);
			output.push_back(sub);
			start=i;
		}
	}
	//final
	string sub=input.substr(start,len-start);
	output.push_back(sub);
}


//=========== MAT_To_PSICOV ===========//
//-> input:  a symmetric contact prediction matrix, with higher value the better prediction
//-> output: PSICOV format, or other constraint format

//---------- input a string, output a vector -----//
int String_To_Vector(string &input,vector <double> &output)
{
	istringstream www(input);
	output.clear();
	int count=0;
	double value;
	for(;;)
	{
		if(! (www>>value) )break;
		output.push_back(value);
		count++;
	}
	return count;
}

//---------- input a symmetric matrix --------//
int Input_Symmetric_Matrix(string &input_file, vector < vector < double > > &output_mat)
{
	//start
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		return -1;
	}
	//load
	int count=0;
	int colnum;
	int colcur;
	int first=1;
	output_mat.clear();
	vector <double> tmp_rec;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		colcur=String_To_Vector(buf,tmp_rec);
		if(first==1)
		{
			first=0;
			colnum=colcur;
		}
		else
		{
			if(colcur!=colnum)
			{
				fprintf(stderr,"current column number %d not equal to the first column number %d \n",
					colcur,colnum);
				return -1;
			}
		}
		output_mat.push_back(tmp_rec);
		count++;
	}
	//final check
	if(count!=colnum)
	{
		fprintf(stderr,"row number %d not equal to column number %d \n",
			count,colnum);
		return -1;
	}
	//symmetric check
	double epsilu=1.e-9;
	for(int i=0;i<count;i++)
	{
		for(int j=i+1;j<count;j++)
		{
			double value=fabs(output_mat[i][j]-output_mat[j][i]);
			if(value>epsilu)
			{
				fprintf(stderr,"WARNING: symmetric condition broken at [%d,%d] -> %lf not equal to %lf \n",
					i,j,output_mat[i][j],output_mat[j][i]);
			}
		}
	}
	//return
	return count;
}

//---- calculate matrix Z_score ------//
void Calc_Mat_Zscore(vector < vector < double > > &input_mat, double &mean, double &vari)
{
	int i,j;
	int size=(int)input_mat.size();
	int count;
	//mean
	mean=0;
	count=0;
	for(i=0;i<size;i++)
	{
		for(j=i+6;j<size;j++)
		{
			mean+=input_mat[i][j];
			count++;
		}
	}
	mean/=count;
	//vari
	vari=0;
	count=0;
	for(i=0;i<size;i++)
	{
		for(j=i+6;j<size;j++)
		{
			vari+=(input_mat[i][j]-mean)*(input_mat[i][j]-mean);
			count++;
		}
	}		
	vari=1.0*sqrt(1.0*vari/count);
}

//--------- symmetric matrix to column stype --------//
//-> note: posi starts from 0 !!
int SymmMat_To_Column(vector < vector < double > > &input_mat,
	vector<pair<int, int> > &posi, vector <double> &value)
{
	int i,j;
	int size=(int)input_mat.size();
	posi.clear();
	value.clear();
	int count=0;
	for(i=0;i<size;i++)
	{
		for(j=i;j<size;j++)
		{
			posi.push_back(pair<int,int>(i,j));
			value.push_back(0.5*(input_mat[i][j]+input_mat[j][i]));
			count++;
		}
	}
	//return
	return count;
}

//----------- vector sort ----------//
//-> UPorDOWN: UP for ascending order sort, DOWN for descending order sort
void Vector_Sort(vector <double> &input, vector <int> &order,int UPorDOWN=0)
{
	order.clear();
	//prepare Fast_Sort
	Fast_Sort <double> fast_sort;
	int size=(int)input.size();
	double *temp_score=new double[size];
	int *temp_index=new int[size];
	for(int i=0;i<size;i++)temp_score[i]=input[i];
	//Fast_Sort
	if(UPorDOWN==1)fast_sort.fast_sort_1up(temp_score,temp_index,size);
	else fast_sort.fast_sort_1(temp_score,temp_index,size);
	//output
	order.resize(size);
	for(int i=0;i<size;i++)order[i]=temp_index[i];
	//delete
	delete [] temp_score;
	delete [] temp_index;
}

//============= shrink column list ===========//
//-> given separation start and end, as well as TopK ratio value
int Shrink_Column_List(int length,
	vector<pair<int, int> > &posi, vector <double> &value, vector <int> &order,
	vector<pair<int, int> > &posi_out, vector <double> &value_out,
	int sepa_start, int sepa_end, double topk_ratio)
{
	int i;
	int size=(int)value.size();
	int count=0;
	int topk=(int)(topk_ratio*length); //topk_ratio should be 30 as default
	posi_out.clear();
	value_out.clear();
	int long_range=0;
	int med_range=0;
	int short_range=0;
	for(i=0;i<size;i++)
	{
		int index=order[i];
		int ii=posi[index].first;
		int jj=posi[index].second;
		if(sepa_start >=0 && sepa_end>=0)
		{
			if(abs(ii-jj)<sepa_start || abs(ii-jj)>sepa_end)continue;
		}
		//count med and long
		int pattern=0;
		if(abs(ii-jj)>=24)long_range++,pattern=1;
		if(abs(ii-jj)>=12 && abs(ii-jj)<=23)med_range++,pattern=2;
		if(abs(ii-jj)<12)short_range++,pattern=3;
		//check topk
		if(topk_ratio>=0)
		{
			if(count>=topk)break;
		}
		//push_back
		if(pattern==1 && long_range<1.0*topk_ratio*length/2) //15L
		{
			posi_out.push_back(posi[index]);
			value_out.push_back(value[index]);
			count++;
		}
		if(pattern==2 && med_range<1.0*topk_ratio*length/3)  //10L
		{
			posi_out.push_back(posi[index]);
			value_out.push_back(value[index]);
			count++;
		}
		if(pattern==3 && short_range<1.0*topk_ratio*length/5) //6L
		{
			posi_out.push_back(posi[index]);
			value_out.push_back(value[index]);
			count++;
		}
	}
	//return
	return count;
}

//------------ example of contact.tbl format ---------//
/*
assign (resid  40 and name cb) (resid  60 and name cb) 3.60 0.10 4.40
assign (resid  24 and name cb) (resid  91 and name cb) 3.60 0.10 4.40
assign (resid  27 and name cb) (resid  93 and name ca) 3.60 0.10 4.40
assign (resid   4 and name cb) (resid  32 and name cb) 3.60 0.10 4.40
assign (resid  41 and name cb) (resid  56 and name cb) 3.60 0.10 4.40
assign (resid  39 and name ca) (resid  80 and name cb) 3.60 0.10 4.40
assign (resid  18 and name cb) (resid  86 and name cb) 3.60 0.10 4.40
assign (resid  41 and name cb) (resid  55 and name cb) 3.60 0.10 4.40
assign (resid  34 and name cb) (resid  62 and name cb) 3.60 0.10 4.40
....
*/
void Output_To_ContactTBL(FILE *fp,string &seqres,
	vector<pair<int, int> > &posi_in, vector <double> &value_in,int start_pos,
	double mean,double vari, double Zsco_cut, double lower_L, double upper_L,
	int protein_length,int limit=31000)
{
	int i;
	int size=(int)value_in.size();
	//-> determine lower_bound and upper_bound
	int lower_bound=protein_length*lower_L;
	int upper_bound=protein_length*upper_L;
	//-> otput
	for(i=0;i<size;i++)
	{
		if(i>=limit)break;
		if(i>upper_bound)break;
		int ii=posi_in[i].first;
		int jj=posi_in[i].second;
		//check Z-score cutoff
		double z_sco=1.0*(value_in[i]-mean)/vari;
		if(i>lower_bound && z_sco<Zsco_cut)break;
		//check Glycine
		string cb1="cb";
		if(seqres[ii]=='G')cb1="ca";
		string cb2="cb";
		if(seqres[jj]=='G')cb2="ca";
		//output to 'contact.tbl'
		fprintf(fp,"assign (resid %3d and name %s) (resid %3d and name %s) 3.60 0.10 4.40 !-> %c %c %lf\n",
			ii+start_pos,cb1.c_str(),jj+start_pos,cb2.c_str(),seqres[ii],seqres[jj],z_sco);
	}
}



//------------ main -------------//
int main(int argc, char** argv)
{
	//------- MAT_To_CaspRR -----//
	{
		if(argc<7)
		{
			fprintf(stderr,"Version 1.00\n");
			fprintf(stderr,"GCNNzsco_To_TBL <input_fasta> <input_symm_mat> <output_tbl> \n");
			fprintf(stderr,"        <Zsco_cutoff> <lower_L> <upper_L> \n");
			fprintf(stderr,"[note]: Zsco_cutoff should be set to 3.0, by default. \n");
			fprintf(stderr,"        lower_L and upper_L should be set to 1 and 3, respectively. \n");
			exit(-1);
		}
		
		string input_fasta=argv[1];
		string input_file=argv[2];
		string output_file=argv[3];
		double Zsco_cut=atof(argv[4]);  //-> Zsco > 3.0
		double lower_L=atof(argv[5]);   //-> lowerL = 1.0
		double upper_L=atof(argv[6]);   //-> upperL = 3.0
		//separation
		int sepa_start=6;
		int sepa_end=99999;
		double topk_ratio=30.0;
		int start_pos=1;
		//input fasta
		string seqres;
		int seq_len=Read_FASTA_SEQRES(input_fasta,seqres);
		if(seq_len<=0)exit(-1);
		//input matrix
		vector < vector < double > > symm_mat;
		int length=Input_Symmetric_Matrix(input_file, symm_mat);
		if(length<=0)exit(-1);
		//check length
		if(seq_len!=length)
		{
			fprintf(stderr,"seq_len %d not equal to matrix_len %d \n",seq_len,length);
			exit(-1);
		}
		//calculate zscore
		double mean,vari;
		Calc_Mat_Zscore(symm_mat, mean, vari);
		//get target name
		string targ_name;
		getBaseName(input_fasta,targ_name,'/','.');
		//process sort
		vector<pair<int, int> > posi;
		vector <double> value;
		SymmMat_To_Column(symm_mat,posi,value);
		vector <int> order;
		Vector_Sort(value, order);
		//process shrink
		vector<pair<int, int> > posi_out;
		vector <double> value_out;
		Shrink_Column_List(length,posi,value,order,posi_out,value_out,
			sepa_start,sepa_end,topk_ratio);
		//output to TBL_file
		FILE *fp=fopen(output_file.c_str(),"wb");
		Output_To_ContactTBL(fp,seqres,posi_out,value_out,start_pos,mean,vari,
			Zsco_cut,lower_L,upper_L,seq_len);
		fclose(fp);
		//exit
		exit(0);
	}
}
