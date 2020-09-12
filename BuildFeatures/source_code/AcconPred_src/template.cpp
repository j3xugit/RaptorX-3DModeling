#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include "template.h"
using namespace std;

//------- constructor & destructor -------//
TEMPLATE::TEMPLATE(string tempName,string root,int LoadOrig,int LoadBlast)
{
	StartPos = -1;
	tLength = -1;
	temp_name = tempName;
	length=0;        //inner length=0;
	WantOriginal=LoadOrig;
	WantBlast=LoadBlast;
	ReadFeatures_TPL(root + "/" + tempName + ".tpl");
	dis_matrix=0;
}
TEMPLATE::~TEMPLATE(void)
{
	if(length!=0)Template_Delete_Matrix(length);
	if(dis_matrix!=0)DeleteArray2D_(&dis_matrix,length);
}

//-> type: 1 for CA and the other for CB
float TEMPLATE::DistOf2AAs(int i, int j, int type)
{
	if(isMissing[i]==0 && isMissing[j]==0)
	{
		if(type==1)return sqrt((CA[i][0]-CA[j][0])*(CA[i][0]-CA[j][0]) + (CA[i][1]-CA[j][1])*(CA[i][1]-CA[j][1]) + (CA[i][2]-CA[j][2])*(CA[i][2]-CA[j][2]));
		else return sqrt((CB[i][0]-CB[j][0])*(CB[i][0]-CB[j][0]) + (CB[i][1]-CB[j][1])*(CB[i][1]-CB[j][1]) + (CB[i][2]-CB[j][2])*(CB[i][2]-CB[j][2]));
	}
	return -1;
}

//----- epad related ---//__140810__//
void TEMPLATE::Compute_All_Distance(void)
{
	//create
	if(dis_matrix==0)NewArray2D_(&dis_matrix,length,length);
	for(int i=0;i<length;i++)for(int j=0;j<length;j++)dis_matrix[i][j] = 0;
	//calculate
	for(int i=0;i<length;i++)
	{
		for(int j=i+1;j<length;j++)
		{
			if(isMissing[i]==0 && isMissing[j]==0)
			{
				dis_matrix[i][j] = DistOf2AAs(i,j,2);
				dis_matrix[j][i] = dis_matrix[i][j];
			}
		}
	}
}


//---------- create & delete --------//
void TEMPLATE::Template_Create_Matrix(int length)
{
	isMissing=new short[length];
	isMultiHIS=new short[length];
	NewArray2D_(&CA,length,3);
	NewArray2D_(&CB,length,3);
	ACCp=new short[length];
	ACC=new short[length];
	CLE=new short[length];
	SS=new short[length];
	NewArray2D_(&SS2,length,3);
	NewArray2D_(&SS8,length,8);
	NewArray2D_(&acc,length,3);
	NewArray2D_(&SS2_,length,3);
	NewArray2D_(&SS8_,length,8);
	NewArray2D_(&acc_,length,3);
	//contact number
	contactNum=new short[length];
	contactNum_B=new short[length];
}
void TEMPLATE::Template_Delete_Matrix(int length)
{
	delete [] isMissing;
	delete [] isMultiHIS;
	DeleteArray2D_(&CA,length);
	DeleteArray2D_(&CB,length);
	delete [] ACCp;
	delete [] ACC;
	delete [] CLE;
	delete [] SS;
	DeleteArray2D_(&SS2,length);
	DeleteArray2D_(&SS8,length);
	DeleteArray2D_(&acc,length);
	DeleteArray2D_(&SS2_,length);
	DeleteArray2D_(&SS8_,length);
	DeleteArray2D_(&acc_,length);
	//contact number
	delete [] contactNum;
	delete [] contactNum_B;
}

//-------- main part ------//
void TEMPLATE::ReadFeatures_TPL(string filename)
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
		//SEQRES
		if(strncmp(buf,"SEQRES",6)==0)
		{
			sequence = (string)(buf+18);
			length = (int)sequence.size();
			if(StartPos < 0)
				StartPos = 0, tLength = length;
			if(tLength <=0 || StartPos + tLength > length)
				tLength = length - StartPos;	
			sequence = sequence.substr(StartPos,tLength);

			//create here !!!! //__120520__//
			Profile_Create_Matrix(tLength);
			Template_Create_Matrix(tLength);
			//create over !!!!//

			//process isMultiHIS//__120320__//
			for(int i=0;i<tLength;i++) isMultiHIS[i] = 0;
			//head mark
			int pos1=sequence.find("HHH");
			if(pos1<10 && pos1!=-1)
			{
				int i;
				for(i=0;i<pos1;i++)isMultiHIS[i] = 1;
				for(i=pos1;i<tLength-1;i++)
				{
					if(sequence[i]=='H')
						isMultiHIS[i]=1;
					else break;
				}
				isMultiHIS[i]=1;
				if(i+1<tLength)isMultiHIS[i+1]=1;
			}
			//tail mark
			int pos2=sequence.rfind("HHH");
			if(pos2>tLength-10 && pos2!=-1)
			{
				int i;
				for(i=tLength-1;i>=pos2;i--)isMultiHIS[i] = 1;
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
				if(i<StartPos) continue;
				if(i>StartPos+tLength-1) continue;
				if(sequence[i-StartPos]<'A'||sequence[i-StartPos]>'Z')residue[i-StartPos]=20;
				else residue[i-StartPos] = (short)AA1Coding[AA2SUB[sequence[i-StartPos]-'A']];
			}
		}//end of SEQRES
		//DSSP
		if(strncmp(buf,"DSSP",4)==0)
		{
			dssp_sequence = (string)(buf+18);
			dssp_sequence = dssp_sequence.substr(StartPos,tLength);
		}
		//NEFF
		if(strncmp(buf,"NEFF",4)==0)
		{
			sscanf(buf+7,"%f",&NEFF);
		}
		//FEAT
		if(strncmp(buf,"//////////// Features",20)==0)
		{
			//skip first
			if(!getline(fin,wbuf,'\n'))
			{
				fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT skip first] \n",filename.c_str());
				failure=0;
				return;
			}
			//get remaining
			for(int i=0;i<length;i++)
			{
				//get first
				if(!getline(fin,wbuf,'\n'))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				if(wbuf.length()<52)
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				if(i<StartPos) continue;
				if(i>StartPos+tLength-1) continue;
				//get features
				stringstream www(wbuf);
				for(int j=0;j<2;j++)  //skip header
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d skip %d] \n",filename.c_str(),i,j);
						failure=0;
						return;
					}
				}
				int next_read=0;
				//[1] missing
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d missing] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				if(temp=="0")
				{
					isMissing[i-StartPos]=0;
					next_read=1;
				}
				else if(temp=="1")
				{
					isMissing[i-StartPos]=1;
					next_read=0;
					CA[i-StartPos][0]=CA[i-StartPos][1]=CA[i-StartPos][2]=0;
					CB[i-StartPos][0]=CB[i-StartPos][1]=CB[i-StartPos][2]=0;
				}
				else
				{
					fprintf(stderr,"CONTENT BAD AT TEMPLATE FILE %s -> [FEAT line %d missing] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				//[2] SSE
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d SSE] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				//ss8
				int aa;
				if(temp=="H")aa=0;
				else if(temp=="G")aa=1;
				else if(temp=="I")aa=2;
				else if(temp=="E")aa=3;
				else if(temp=="B")aa=4;
				else if(temp=="T")aa=5;
				else if(temp=="S")aa=6;
				else if(temp=="L")aa=7;
				else
				{
					fprintf(stderr,"CONTENT BAD AT TEMPLATE FILE %s -> [FEAT line %d SSE] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				for(int nn=0;nn<8;nn++)
				{
					if(aa==nn)SS8[i-StartPos][nn]=1;
					else SS8[i-StartPos][nn]=0;
				}
				//ss3
				if(temp=="H")
				{
					SS[i-StartPos]=HELIX;
					SS2[i-StartPos][0]=0;
					SS2[i-StartPos][1]=1;
					SS2[i-StartPos][2]=0;
				}
				else if(temp=="E")
				{
					SS[i-StartPos]=SHEET;
					SS2[i-StartPos][0]=0;
					SS2[i-StartPos][1]=0;
					SS2[i-StartPos][2]=1;
				}
				else
				{
					SS[i-StartPos]=LOOP;
					SS2[i-StartPos][0]=1;
					SS2[i-StartPos][1]=0;
					SS2[i-StartPos][2]=0;
				}
				//[3] CLE
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d CLE] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				CLE[i-i-StartPos]=temp[0]-'A';
				//[4] ACC
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d ACC] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				if(temp=="0"){
					ACC[i-StartPos] = temp[0]-'0';
					acc[i-StartPos][0]=1;
					acc[i-StartPos][1]=0;
					acc[i-StartPos][2]=0;
				}
				else if(temp=="1"){
					ACC[i-StartPos] = temp[0]-'0';
					acc[i-StartPos][0]=0;
					acc[i-StartPos][1]=1;
					acc[i-StartPos][2]=0;
				}
				else if(temp=="2"){
					ACC[i-StartPos] = temp[0]-'0';
					acc[i-StartPos][0]=0;
					acc[i-StartPos][1]=0;
					acc[i-StartPos][2]=1;
				}
				else
				{
					fprintf(stderr,"CONTENT BAD AT TEMPLATE FILE %s -> [FEAT line %d ACC] \n",filename.c_str(),i);
					failure=0;
					return;
				}

				//[5] pACC
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d pACC] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				int pacc=atoi(temp.c_str());
				if(pacc<0 || pacc>100)
				{
					fprintf(stderr,"CONTENT BAD AT TEMPLATE FILE %s -> [FEAT line %d pACC %d] \n",filename.c_str(),i,pacc);
					failure=0;
					return;
				}
				ACCp[i-StartPos]=pacc;
				//[6] CA_contact
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d CA_contact] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				int ca_contact=atoi(temp.c_str());
				if(ca_contact<0 || ca_contact>20)
				{
					fprintf(stderr,"CONTENT BAD AT TEMPLATE FILE %s -> [FEAT line %d CA_contact %d] \n",filename.c_str(),i,ca_contact);
					failure=0;
					return;
				}
				contactNum[i-StartPos]=ca_contact;
				//[7] CB_contact
				if(!(www>>temp))
				{
					fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d CB_contact] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				int cb_contact=atoi(temp.c_str());
				if(cb_contact<0 || cb_contact>20)
				{
					fprintf(stderr,"CONTENT BAD AT TEMPLATE FILE %s -> [FEAT line %d CB_contact %d] \n",filename.c_str(),i,cb_contact);
					failure=0;
					return;
				}
				contactNum_B[i-StartPos]=cb_contact;
				//[8] CA,CB coordinate
				if(next_read==1)
				{
					//CA
					for(int j=0;j<3;j++)
					{
						if(!(www>>temp))
						{
							fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d CA coordinate] \n",filename.c_str(),i);
							failure=0;
							return;
						}
						float aa=atof(temp.c_str());
						CA[i-StartPos][j]=aa;
					}
					//CB
					for(int j=0;j<3;j++)
					{
						if(!(www>>temp))
						{
							fprintf(stderr,"FORMAT BAD AT TEMPLATE FILE %s -> [FEAT line %d CA coordinate] \n",filename.c_str(),i);
							failure=0;
							return;
						}
						float aa=atof(temp.c_str());
						CB[i-StartPos][j]=aa;
					}
					//check
					int check=0;
					if(check==1)
					{
						float x1=CA[i-StartPos][0]-CB[i-StartPos][0];
						float x2=CA[i-StartPos][1]-CB[i-StartPos][1];
						float x3=CA[i-StartPos][2]-CB[i-StartPos][2];
						float dist=sqrt(x1*x1+x2*x2+x3*x3);
						if(dist>2.5)
						{
							contactNum_B[i-StartPos]=contactNum[i-StartPos]; //modified at //__2012_05_18__// adviced by Jinbo
						}
					}
				}

			}
		}//end of FEAT
		//SS2
		if(strncmp(buf,"//////////// Original SS2",25)==0)
		{
			//skip first
			if(!getline(fin,wbuf,'\n'))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS2 skip first] \n",filename.c_str());
				failure=0;
				return;
			}
			//get remaining
			for(int i=0;i<length;i++)
			{
				//get first
				if(!getline(fin,wbuf,'\n'))
				{
					fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS2 line %d] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				if(wbuf.length()<6)
				{
					fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS2 line %d] \n",filename.c_str(),i);
					failure=0;
					return;
				}
				//process string
				istringstream www(wbuf);
				for(int j=0;j<3;j++) //skip header
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS2 line %d header] \n",filename.c_str(),i);
						failure=0;
						return;
					}
				}
				for(int j=0;j<3;j++) //get remaining
				{
					if(!(www>>temp))
					{
						fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [SS2 line %d col %d] \n",filename.c_str(),i,j);
						failure=0;
						return;
					}
				}
			}
		}//end of SS2
		//PSM
		if(strncmp(buf,"//////////// Original PSM",25)==0)
		{
			Process_PSM(fin,filename,StartPos);
		}//end of PSM
		//PSP
		if(strncmp(buf,"//////////// Original PSP",25)==0)
		{
			Process_PSP(fin,filename,StartPos);
		}//end of PSP
		//HMM
		if(strncmp(buf,"//////////// Original HHM",25)==0)
		{
			Process_HMM(fin,filename,StartPos);
		}//end of HHM
	}//end of for(;;)
	length = tLength;
}
