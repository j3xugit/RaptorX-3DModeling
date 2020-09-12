#include <string>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <sys/types.h> 
#include <unistd.h>

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

//------ Random Number Generator ------//
int Random_Number_Generator(void)
{

	/*
	srand( time( NULL ) );
	return rand();
	*/

	timeval t1;
	gettimeofday(&t1, NULL);
	//srand(t1.tv_usec * t1.tv_sec + getpid() );
	srand(t1.tv_usec * getpid() + t1.tv_sec );
	return rand();
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

//==================== Note that we only need to change the below sections ===============//
//-> 1. model number
/*
{* if only regularizing coordinate files (no DG) then 
   enter the number of coordinate files to be regularized (DGSA) *}
{===>} pdb.dg.count=10;
{* number of trial or accepted structures *}
{===>} pdb.end.count=10;
*/

//-> 2. random seed (should be generated automatically)
/*
{* seed for random number generator *}
{* change to get different initial velocities *}
{===>} md.seed=1668066805;
*/

//-> 3. output file
/*
{* base name for output coordinate files *}
{===>} pdb.out.name="1_1pazA";
*/

//-> 4. minimization step (= L*15, where L is the protein length)
/*
{* number of minimization steps *}
{===>} md.pow.step=1845;
*/

//-> 5. parameters for final minimization step
/*
{* scale factor for NOE energy term *}
{===>} md.pow.noe=10;
{* scale factor for dihedral angle energy term *}
{===>} md.pow.cdih=50;     //--> ws_note: 50 = 10*5
*/

//-> 6. parameters for the distance geometry stage
/*
{* scale factor for distance geometry restraint term *}
{===>} md.dg.scale=100.;  //--? not sure to modify or not
{* high temperature scale factor for dihedral angle energy term *}
{===>} md.hot.cdih=5;     //--? not sure to modify or not
{* scale factor for NOE energy term *}
{===>} md.cool.noe=10;
{* slow-cooling scale factor for dihedral angle energy term *}
{===>} md.cool.cdih=50;    //--> ws_note: 50 = 10*5
*/

//-> 7. scaling factor for sse
/*
@CNS_CUSTOMMODULE:scalehotedited ( md=&md;
                          nmr=&nmr;
                          wscale=5;
                          input.noe.scale=&md.cool.noe;
                          input.cdih.scale=&md.hot.cdih; )
      @CNS_CUSTOMMODULE:scalecoolsetupedited ( kdani=$kdani;
                                      ksani=$ksani;
                                      nmr=&nmr;
                                      wscale=5;
                                      input.noe.scale=&md.cool.noe;
                                      input.cdih.scale=&md.cool.cdih;
                                      input.ncycle=$ncycle; )
*/

//============== Main Process ==========//
//-> [note]: NOE_scale shall be set to 10 to 100,
//           Dihedral_scale shall be set to 1 to 10, and less than NOE_scale
//-> [default]: NOE_scale = 10, Dihedral_scale = 5,
void DSGA_File_Mod(string &input_file, string &output_file,
	int model_num, int rand_seed, string &out_name, int protein_length, 
	double NOE_scale, double Dihedral_scale, double initNOE_scale)
{
	ifstream fin;
	string buf,temp;
	FILE *fp=fopen(output_file.c_str(),"wb");
	//read
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",input_file.c_str());
		exit(-1);
	}
	//modify parameter
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		//-> model number
		{
			temp=buf.substr(0,20);
			if(temp=="{===>} pdb.dg.count=")
			{
				fprintf(fp,"{===>} pdb.dg.count=%d;\n",model_num);
				continue;
			}
			temp=buf.substr(0,21);
			if(temp=="{===>} pdb.end.count=")
			{
				fprintf(fp,"{===>} pdb.end.count=%d;\n",model_num);
				continue;
			}
		}
		//-> random seed
		{
			temp=buf.substr(0,15);
			if(temp=="{===>} md.seed=")
			{
				int wsrand;
				if(rand_seed<=0)wsrand=Random_Number_Generator();
				else wsrand=rand_seed;
				fprintf(fp,"{===>} md.seed=%d;\n",wsrand);
				continue;
			}
		}
		//-> output file
		{
			temp=buf.substr(0,20);
			if(temp=="{===>} pdb.out.name=")
			{
				fprintf(fp,"{===>} pdb.out.name=\"%s\";\n",out_name.c_str());
				continue;
			}
		}
		//-> minimization step
		{
			temp=buf.substr(0,19);
			if(temp=="{===>} md.pow.step=")
			{
				fprintf(fp,"{===>} md.pow.step=%d;\n",15*protein_length);
				continue;
			}
		}
		//-> parameters for final minimization step
		{
			temp=buf.substr(0,18);
			if(temp=="{===>} md.pow.noe=")
			{
				fprintf(fp,"{===>} md.pow.noe=%f;\n",NOE_scale);
				continue;
			}
			temp=buf.substr(0,19);
			if(temp=="{===>} md.pow.cdih=")
			{
				//modified by Jinbo Xu
				//fprintf(fp,"{===>} md.pow.cdih=%f;\n",NOE_scale*Dihedral_scale);
				fprintf(fp,"{===>} md.pow.cdih=%f;\n",initNOE_scale*Dihedral_scale);
				continue;
			}
		}
		//-> parameters for the distance geometry stage
		{
			temp=buf.substr(0,19);
			if(temp=="{===>} md.cool.noe=")
			{
				fprintf(fp,"{===>} md.cool.noe=%f;\n",NOE_scale);
				continue;
			}
			temp=buf.substr(0,20);
			if(temp=="{===>} md.cool.cdih=")
			{
				//modified by Jinbo Xu
				//fprintf(fp,"{===>} md.cool.cdih=%f;\n",NOE_scale*Dihedral_scale);
				fprintf(fp,"{===>} md.cool.cdih=%f;\n",initNOE_scale*Dihedral_scale);
				continue;
			}
		}
		//-> scaling factor for sse
		{
			temp=buf.substr(0,33);
			if(temp=="                          wscale=")
			{
				fprintf(fp,"                          wscale=%f;\n",Dihedral_scale);
				continue;
			}
			temp=buf.substr(0,45);
			if(temp=="                                      wscale=")
			{
				fprintf(fp,"                                      wscale=%f;\n",Dihedral_scale);
				continue;
			}
		}
		//-> normal output
		fprintf(fp,"%s\n",buf.c_str());
	}
}


//----------- main -------------//
int main(int argc,char **argv)
{
	//---- DGSA_File_Mod ----//__170105__//
	{
		if(argc<10)
		{
			fprintf(stderr,"DGSA_File_Mod <fasta> <raw_dgsa> <out_dgsa> <model_num> <rand_seed> <out_name> \n");
			fprintf(stderr,"              <NOE_scale> <Dihedral_scale> <initNOE_scale>\n");
			fprintf(stderr,"[note1]: <fasta> and <raw_dgsa> are two input files, \n");
			fprintf(stderr,"         for protein length and taw 'dgsa.inp', respectively. \n");
			fprintf(stderr,"         <out_dgsa> are the output file for modified 'dgsa.inp' \n");
			fprintf(stderr,"[note2]: <model_num>, <rand_seed>, <out_name> are three basic parameters, \n");
			fprintf(stderr,"         if rand_seed is set to -1, then automatically generate one. \n");
			fprintf(stderr,"[note3]: <NOE_scale> <Dihedral_scale> are two advanced parameters, \n");
			fprintf(stderr,"         for distance restraints and dihedral restraints, respectively. \n");
			fprintf(stderr,"         by default, NOE_scale shall set to 10, and Dihedral_scale to 5, \n");
			fprintf(stderr,"         which means that the weight of dihedral restraints are 10*5=50 \n");
			exit(-1);
		}

		//----- argument -----//
		//-> input/output
		string fasta=argv[1];
		string raw_dgsa=argv[2];
		string out_dgsa=argv[3];
		//-> basic parameter
		int model_num=atoi(argv[4]);
		int rand_seed=atoi(argv[5]);
		string out_name=argv[6];
		//-> advanced parameter
		double NOE_scale=atof(argv[7]);
		double Dihedral_scale=atof(argv[8]);
		double initNOE_scale=atof(argv[9]);

		//----- parameter check ----//
		if(model_num<1)
		{
			fprintf(stderr,"model_num %d is not valid \n",model_num);
			exit(-1);
		}
		if(initNOE_scale<0 || initNOE_scale>100)
		{
			fprintf(stderr,"initNOE_scale %lf must be restricted to [0,100] \n",initNOE_scale);
			exit(-1);
		}
		if(Dihedral_scale<0.1 || Dihedral_scale>10)
		{
			fprintf(stderr,"Dihedral_scale %lf must be restricted to [0.1 , 10] \n",Dihedral_scale);
			exit(-1);
		}
		if(NOE_scale<0 || NOE_scale>100)
		{
			fprintf(stderr,"NOE_scale %lf must be restricted to [0,100] \n",NOE_scale);
			exit(-1);
		}

		//----- load fasta sequence for protein length ----//
		string seqres;
		int length=Read_FASTA_SEQRES(fasta,seqres);

		//----- main process to modify the parameters in 'dgsa.inp' ---//
		DSGA_File_Mod(raw_dgsa,out_dgsa,model_num,rand_seed,out_name,length,
			NOE_scale,Dihedral_scale, initNOE_scale);

		//exit
		exit(0);
	}
}

