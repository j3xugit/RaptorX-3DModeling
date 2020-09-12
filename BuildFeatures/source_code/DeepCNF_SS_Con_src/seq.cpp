#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include "seq.h"
using namespace std;

//------- constructor & destructor -------//
SEQUENCE::SEQUENCE(string seqname,string tgt_root,int LoadOrig,int LoadBlast)
{
	seq_name = seqname;
	length=0;        //inner length=0;
	WantOriginal=LoadOrig;
	WantBlast=LoadBlast;
	ReadFeatures_TGT(tgt_root + "/" + seqname+".tgt");
	pair_dis=0;
}
SEQUENCE::~SEQUENCE(void)
{
	if(length!=0)Sequence_Delete_Matrix(length);
	if(pair_dis!=0)DeleteArray3D_(&pair_dis,length,length);
}

//--------- distance related --------//
void SEQUENCE::Read_distance_file(string filename)
{
	//create
	if(pair_dis==0)NewArray3D_(&pair_dis,length,length,12);
	for(int i=0;i<length;i++)for(int j=0;j<length;j++)for(int k=0;k<12;k++)pair_dis[i][j][k] = 0;
	//init
	int pair_count = 0;
	char buf[40960];
	ifstream fea_in(filename.c_str());
	//read
	if(!fea_in.is_open())
	{
		cerr << "Epad file missing " << filename << endl;
		failure=0;
		return;
	}
	//process
	int index1,index2;
	double dis_prob[12],temp_dis=-1,sum_dis=0,max_dis=0;
	max_energy=0,min_energy=0;
	while(fea_in.getline(buf,40960))
	{
		sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&index1,&index2,&dis_prob[0],&dis_prob[1],&dis_prob[2],&dis_prob[3],&dis_prob[4],&dis_prob[5],
			&dis_prob[6],&dis_prob[7],&dis_prob[8],&dis_prob[9],&dis_prob[10],&dis_prob[11]);
		pair_count++;
		sum_dis = 0; temp_dis = -1; max_dis =0;
		for(int i=0;i<12;i++)
		{
			if(dis_prob[i]==0) dis_prob[i] = 0.01;
			pair_dis[index1][index2][i] = 0.01*dis_prob[i];
			pair_dis[index2][index1][i] = pair_dis[index1][index2][i];
		}
	}
	fea_in.close();
}

//======================= Binary epad_prob file process ===================//__140901__//
void SEQUENCE::Read_distance_file_Binary(double *input,int size,int line_)
{
	//-> 0. init
	//create
	if(pair_dis==0)NewArray3D_(&pair_dis,length,length,12);
	for(int i=0;i<length;i++)for(int j=0;j<length;j++)for(int k=0;k<12;k++)pair_dis[i][j][k] = 0;
	//init
	int pair_count = 0;
	int index1,index2;
	double dis_prob[12],temp_dis=-1,sum_dis=0,max_dis=0;
	max_energy=0,min_energy=0;
	//judge
	if(size<=0 || line_<=0)return;

	//-> 1. calculate protein length
	double epsilu=0.1;
	int len=(int)(1.0*sqrt(0.25+2.0*line_)-0.5+epsilu)+6;
	//-> 2. output
	int first=1;
	int line=-1;
	int cur=0;
	int separate=6;
	int increase=0;
	string cur_rec;
	char command[30000];
	for(int i=0;i<size;i++)
	{
		//increase line
		if(i%12==0)
		{
			line++;
			if(first==1)
			{
				first=0;
				sprintf(command,"%d %d ",cur,cur+separate+increase);
				cur_rec+=command;
				increase++;
			}
			else
			{
				//judge increase
				if(cur+separate+increase>=len)
				{
					cur++;
					increase=0;
				}
				//output
				{
					sscanf(cur_rec.c_str(),"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						&index1,&index2,&dis_prob[0],&dis_prob[1],&dis_prob[2],&dis_prob[3],&dis_prob[4],&dis_prob[5],
						&dis_prob[6],&dis_prob[7],&dis_prob[8],&dis_prob[9],&dis_prob[10],&dis_prob[11]);
					pair_count++;
					sum_dis = 0; temp_dis = -1; max_dis =0;
					for(int k=0;k<12;k++)
					{
						if(dis_prob[k]==0) dis_prob[k] = 0.01;
						pair_dis[index1][index2][k] = 0.01*dis_prob[k];
						pair_dis[index2][index1][k] = pair_dis[index1][index2][k];
					}
				}
				cur_rec="";
				sprintf(command,"%d %d ",cur,cur+separate+increase);
				cur_rec+=command;
				increase++;
			}
		}
		//output value
		sprintf(command,"%.4f ",input[i]);
		cur_rec+=command;
	}
	//output
	{
		sscanf(cur_rec.c_str(),"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&index1,&index2,&dis_prob[0],&dis_prob[1],&dis_prob[2],&dis_prob[3],&dis_prob[4],&dis_prob[5],
			&dis_prob[6],&dis_prob[7],&dis_prob[8],&dis_prob[9],&dis_prob[10],&dis_prob[11]);
		pair_count++;
		sum_dis = 0; temp_dis = -1; max_dis =0;
		for(int k=0;k<12;k++)
		{
			if(dis_prob[k]==0) dis_prob[k] = 0.01;
			pair_dis[index1][index2][k] = 0.01*dis_prob[k];
			pair_dis[index2][index1][k] = pair_dis[index1][index2][k];
		}
	}
}


//---------- create & delete --------//
void SEQUENCE::Sequence_Create_Matrix(int length)
{
	isMultiHIS=new short[length];
	ACC=new short[length];
	SS=new short[length];
	DISO=new float[length];
	NewArray2D_(&SS2,length,3);
	NewArray2D_(&SS8,length,8);
	NewArray2D_(&acc_our_10_42,length,3);
	NewArray2D_(&SS2_,length,3);
	NewArray2D_(&SS8_,length,8);
	NewArray2D_(&acc_our_10_42_,length,3);
}
void SEQUENCE::Sequence_Delete_Matrix(int length)
{
	delete [] isMultiHIS;
	delete [] ACC;
	delete [] SS;
	delete [] DISO;
	DeleteArray2D_(&SS2,length);
	DeleteArray2D_(&SS8,length);
	DeleteArray2D_(&acc_our_10_42,length);
	DeleteArray2D_(&SS2_,length);
	DeleteArray2D_(&SS8_,length);
	DeleteArray2D_(&acc_our_10_42_,length);
}

//-------- main part ------//
void SEQUENCE::ReadFeatures_TGT(string filename)
{
	ifstream fin;
	string wbuf,temp;
	char buf[4096];
	fin.open(filename.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",filename.c_str());
		failure=0;
		return;
	}
	//read detail
	for(;;)
	{
		if(!getline(fin,wbuf,'\n'))break;
		if(wbuf.length()>4095)
		{
			temp=wbuf.substr(0,4095);
			strcpy(buf,temp.c_str());
		}
		else strcpy(buf,wbuf.c_str());
		//check sequence
		if(strncmp(buf,"Sequence =",10)==0)
		{
			sequence = (string)(buf+11);
			length = (int)sequence.size();

			//create here !!!! //__120520__//
			Profile_Create_Matrix(length);
			Sequence_Create_Matrix(length);
			//create over !!!!//
	
			//init diso
			for(int i=0;i<length;i++) DISO[i]=0;

			//process isMultiHIS//__120320__//
			for(int i=0;i<length;i++) isMultiHIS[i] = 0;
			//head mark
			int pos1=sequence.find("HHH");
			if(pos1<10 && pos1!=-1)
			{
				int i;
				for(i=0;i<pos1;i++)isMultiHIS[i] = 1;
				for(i=pos1;i<length-1;i++)
				{
					if(sequence[i]=='H')
						isMultiHIS[i]=1;
					else break;
				}
				isMultiHIS[i]=1;
				if(i+1<length)isMultiHIS[i+1]=1;
			}
			//tail mark
			int pos2=sequence.rfind("HHH");
			if(pos2>length-10 && pos2!=-1)
			{
				int i;
				for(i=length-1;i>=pos2;i--)isMultiHIS[i] = 1;
				for(i=pos2-1;i>=1;i--)
				{
					if(sequence[i]=='H')
						isMultiHIS[i]=1;
					else break;
				}
				isMultiHIS[i]=1;
				if(i-1>=0)isMultiHIS[i-1]=1;
			}
			//process isMultiHIS//__120320__//over
			for(int i=0;i<length;i++)
			{
				if(sequence[i]<'A'||sequence[i]>'Z')residue[i]=20;
				else residue[i] = (short)AA1Coding[AA2SUB[sequence[i]-'A']];
			}
		}//end of Sequence
		if(strncmp(buf,"SSEseq   =",10)==0)
		{
			sse_seq = (string)(buf+11);
		}
		if(strncmp(buf,"SSEconf   =",10)==0)
		{
			sse_conf = (string)(buf+11);
		}
		if(strncmp(buf,"ACCseq   =",10)==0)
		{
			acc_seq = (string)(buf+11);
		}
		if(strncmp(buf,"ACCconf   =",10)==0)
		{
			acc_conf = (string)(buf+11);
		}
		//NEFF
		if(strncmp(buf,"NEFF",4)==0)
		{
			sscanf(buf+7,"%f",&NEFF);
		}
		//EVD
		if(strncmp(buf,"EVD",3)==0)
		{
			int curlen=(int)wbuf.size();
			if(curlen==6)
			{
				evd[0]=0;
				evd[1]=1;
			}
			else
			{
				string wwtemp=wbuf.substr(6,curlen);
				istringstream www(wwtemp);
				for(int j=0;j<2;j++) //get remaining
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [EVD col %d] \n",filename.c_str(),j);
						failure=0;
						return;
					}
					evd[j]=atof(temp.c_str());
				}
			}
		}
		//structure
		if(strncmp(buf,"//////////// Original SS3+SS8+ACC file",38)==0)
		{
			//skip first
			if(!getline(fin,wbuf,'\n'))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [structure skip first] \n",filename.c_str());
				failure=0;
				return;
			}
			//skip first
			if(!getline(fin,wbuf,'\n'))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [structure skip first] \n",filename.c_str());
				failure=0;
				return;
			}
			//get remaining
			for(int i=0;i<length;i++)
			{
				//get first
				if(!getline(fin,wbuf,'\n'))
				{
					fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [structure line %d] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				//process string
				istringstream www(wbuf);
				//[ss2]
				for(int j=0;j<3;j++) //get remaining
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS2 line %d col %d] \n",filename.c_str(),i,j);
						failure=0;
						return;
					}
					float wsval=atof(temp.c_str());
					SS2[i][j]= wsval;
				}
				//assign SSE
				if(SS2[i][0]>SS2[i][1] && SS2[i][0] > SS2[i][2])SS[i] = HELIX;
				else if(SS2[i][1]>SS2[i][0] && SS2[i][1]>SS2[i][2])SS[i] = SHEET;
				else SS[i] = LOOP;
				//[ss8]
				for(int j=0;j<8;j++) //get remaining
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS8 line %d col %d] \n",filename.c_str(),i,j);
						failure=0;
						return;
					}
					float wsval=atof(temp.c_str());
					SS8[i][j]= wsval;
				}
				//[acc]
				for(int j=0;j<3;j++) //get remaining
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [ACC line %d col %d] \n",filename.c_str(),i,j);
						failure=0;
						return;
					}
					float wsval=atof(temp.c_str());
					acc_our_10_42[i][j]= wsval;
					if(acc_our_10_42[i][j]<=0.01)
						acc_our_10_42[i][j]=0;
				}
				for(int j=0;j<3;j++)
				{
					if(acc_our_10_42[i][j]>=0.95)
					{
						for(int k=0;k<3;k++)
						{
							if(k==j)
								acc_our_10_42[i][k]=1;
							else
								acc_our_10_42[i][k]=0;
						}
						break;
					}
				}
				//assign ACC
				if(acc_our_10_42[i][0]>acc_our_10_42[i][1] && acc_our_10_42[i][0] > acc_our_10_42[i][2])ACC[i] = 0;
				else if(acc_our_10_42[i][1]>acc_our_10_42[i][0] && acc_our_10_42[i][1]>acc_our_10_42[i][2])ACC[i] = 1;
				else ACC[i] = 2;
			}
		}//structure
		//PSM
		if(strncmp(buf,"//////////// Original PSM",25)==0)
		{
			Process_PSM(fin,filename);
		}//end of PSM
		//PSP
		if(strncmp(buf,"//////////// Original PSP",25)==0)
		{
			Process_PSP(fin,filename);
		}//end of PSP
		//HMM
		if(strncmp(buf,"//////////// Original HHM",25)==0)
		{
			Process_HMM(fin,filename);
		}//end of HHM
		//DISO
		if(strncmp(buf,"//////////// Original DIS",25)==0)
		{
			//skip first
			if(!getline(fin,wbuf,'\n'))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [diso skip first] \n",filename.c_str());
				failure=0;
				return;
			}
			for(int i=0;i<length;i++)
			{
				//get first
				if(!getline(fin,wbuf,'\n'))
				{
					fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [diso line %d] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				//process string
				istringstream www(wbuf);
				//->skip
				for(int j=0;j<5;j++) //get remaining
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [diso line %d col %d] \n",filename.c_str(),i,j);
						failure=0;
						return;
					}
				}
				//->record
				float wsval=atof(temp.c_str());
				DISO[i]=wsval;
			}
		}
	}//end of for(;;)
}

