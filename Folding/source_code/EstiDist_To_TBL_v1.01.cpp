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

//============== estimate distance using predicted probability ===========//
void Estimate_Distance_Bound(int length,
	vector <vector <double> > & estimated_mat, vector <vector <double> > & estimated_var,
	vector <vector <double> > & lower_bound, vector <vector <double> > & upper_bound,
	vector <vector <double> > & lower_dist, vector <vector <double> > & upper_dist,
	double Sigma_Val)
{
	//---- init output matrix ----//start
	{
		lower_dist.resize(length);
		upper_dist.resize(length);
		for(int i=0;i<length;i++)
		{
			lower_dist[i].resize(length);
			upper_dist[i].resize(length);
		}
		for(int i=0;i<length;i++)
		{
			for(int j=0;j<length;j++)
			{
				lower_dist[i][j]=-1;
				upper_dist[i][j]=-1;
			}
		}
	}
	//---- init output matrix ----//end
	
	int i,j;
	for(i=0;i<length;i++)
	{
		for(j=i+6;j<length;j++)
		{
			//---- check for invalid position ----//
			if(estimated_mat[i][j]<0)continue;

			//---- assign lower_dist and upper_dist ---//
			double estimated_dist=estimated_mat[i][j];
			double estimated_std1=lower_bound[i][j];
			double estimated_std2=upper_bound[i][j];

			//----- assign final lower and upper bound ----//
			double lower_value,upper_value;
			if(Sigma_Val<0)
			{
				if(estimated_dist<6)
				{
					lower_value=estimated_dist-3.0*estimated_std1;
					upper_value=estimated_dist+1.0*estimated_std2;
				}
				else if(estimated_dist<8)
				{
					lower_value=estimated_dist-2.5*estimated_std1;
					upper_value=estimated_dist+1.0*estimated_std2;
				}
				else if(estimated_dist<10)
				{
					lower_value=estimated_dist-2.0*estimated_std1;
					upper_value=estimated_dist+1.0*estimated_std2;
				}
				else if(estimated_dist<12)
				{
					lower_value=estimated_dist-2.0*estimated_std1;
					upper_value=estimated_dist+1.5*estimated_std2;
				}
				else
				{
					lower_value=estimated_dist-2.0*estimated_std1;
					upper_value=estimated_dist+2.0*estimated_std2;
				}
			}
			else
			{
				lower_value=estimated_dist-Sigma_Val*estimated_std1;
				upper_value=estimated_dist+Sigma_Val*estimated_std2;
			}

			//----- final bound limit assign ----//
			if(lower_value<3.5)lower_value=3.5;
			if(upper_value>20.0)upper_value=20.0;
			lower_dist[i][j]=lower_value;
			upper_dist[i][j]=upper_value;
			lower_dist[j][i]=lower_value;
			upper_dist[j][i]=upper_value;
		}
	}
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
void Output_To_DistanceTBL(FILE *fp,string &seqres,vector <vector <double> > & estimated_mat,
	vector <vector <double> > & lower_dist, vector <vector <double> > & upper_dist)
{
	int i,j;
	int length=(int)seqres.size();
	int start_pos=1;
	//-> otput
	for(i=0;i<length;i++)
	{
		for(j=i;j<length;j++)
		{
			//get distance bound
			double estimated_dist=estimated_mat[i][j];
			double lower_value=lower_dist[i][j];
			double upper_value=upper_dist[i][j];
			//check distance
			if(estimated_dist<0)continue;
			//check Glycine
			string cb1="cb";
			if(seqres[i]=='G')cb1="ca";
			string cb2="cb";
			if(seqres[j]=='G')cb2="ca";
			//output to 'contact.tbl'
			fprintf(fp,"assign (resid %3d and name %s) (resid %3d and name %s) %f %f %f \n",
				i+start_pos,cb1.c_str(),j+start_pos,cb2.c_str(),
				estimated_dist,estimated_dist-lower_value,upper_value-estimated_dist);
		}
	}
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//------- EstiDist_To_TBL -----//
	{
		if(argc<7)
		{
			fprintf(stderr,"Version 1.01\n");
			fprintf(stderr,"EstiDist_To_TBL <input_fasta> <estimate_dist> <estimate_var> \n");
			fprintf(stderr,"        <lower_std_mat> <upper_std_mat> <output_TBL> <sigma> \n");
			fprintf(stderr,"[note]: if sigma is set to -1, then use automatic estimation, \n");
			fprintf(stderr,"        otherwise, use dist +/- sigma*std to estimate distance. \n");
			exit(-1);
		}
		string input_fasta=argv[1];
		string estimate_dist=argv[2];
		string estimate_vari=argv[3];
		string lower_std_mat=argv[4];
		string upper_std_mat=argv[5];
		string output_TBL=argv[6];
		double sigma=atof(argv[7]);

		//input fasta
		string seqres;
		int seq_len=Read_FASTA_SEQRES(input_fasta,seqres);
		if(seq_len<=0)exit(-1);

		//input matrix
		//-> estimated distance matrix
		vector < vector < double > > estimated_mat;
		int length1=Input_Symmetric_Matrix(estimate_dist, estimated_mat);
		if(length1<=0)exit(-1);
		//check length
		if(seq_len!=length1)
		{
			fprintf(stderr,"seq_len %d not equal to estimated_mat %d \n",seq_len,length1);
			exit(-1);
		}
		//-> estimated variance matrix
		vector < vector < double > > estimated_var;
		int length2=Input_Symmetric_Matrix(estimate_vari, estimated_var);
		if(length2<=0)exit(-1);
		//check length
		if(seq_len!=length2)
		{
			fprintf(stderr,"seq_len %d not equal to estimated_var %d \n",seq_len,length2);
			exit(-1);
		}
		//-> estimated lower variance matrix
		vector < vector < double > > lower_bound;
		int length3=Input_Symmetric_Matrix(lower_std_mat, lower_bound);
		if(length3<=0)exit(-1);
		//check length
		if(seq_len!=length3)
		{
			fprintf(stderr,"seq_len %d not equal to lower_bound %d \n",seq_len,length3);
			exit(-1);
		}
		//-> estimated upper variance matrix
		vector < vector < double > > upper_bound;
		int length4=Input_Symmetric_Matrix(upper_std_mat, upper_bound);
		if(length4<=0)exit(-1);
		//check length
		if(seq_len!=length4)
		{
			fprintf(stderr,"seq_len %d not equal to upper_bound %d \n",seq_len,length4);
			exit(-1);
		}

		//generate lower_dist and upper_dist
		vector <vector <double> > lower_dist;
		vector <vector <double> > upper_dist;
		Estimate_Distance_Bound(seq_len,estimated_mat,estimated_var,lower_bound,upper_bound,
			lower_dist,upper_dist,sigma);

		//output to TBL_file
		FILE *fp=fopen(output_TBL.c_str(),"wb");
		Output_To_DistanceTBL(fp,seqres,estimated_mat,lower_dist,upper_dist);
		fclose(fp);
		//exit
		exit(0);
	}
}
