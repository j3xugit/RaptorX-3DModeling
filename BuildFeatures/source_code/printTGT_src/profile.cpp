#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include "profile.h"
using namespace std;


//------- required matrix ------//
const float gonnet[20][20]= 
{ { 1.7378,  0.870964,0.933254,0.933254, 1.12202,  0.954993, 1,        1.12202,  0.831764, 0.831764,  0.758578, 0.912011, 0.851138, 0.588844, 1.07152,  1.28825,  1.14815,   0.436516,  0.60256,  1.02329},
  { 0.870964,2.95121, 1.07152, 0.933254, 0.60256,  1.41254,  1.09648,  0.794328, 1.14815,  0.57544,   0.60256,  1.86209,  0.676083, 0.47863,  0.812831, 0.954993, 0.954993,  0.691831,  0.660693, 0.630957},
  { 0.933254,1.07152, 2.39883, 1.65959,  0.660693, 1.1749,   1.23027,  1.09648,  1.31826,  0.524807,  0.501187, 1.20226,  0.60256,  0.489779, 0.812831, 1.23027,  1.12202,   0.436516,  0.724436, 0.60256},
  { 0.933254,0.933254,1.65959, 2.95121,  0.47863,  1.23027,  1.86209,  1.02329,  1.09648,  0.416869,  0.398107, 1.12202,  0.501187, 0.354813, 0.851138, 1.12202,  1,         0.301995,  0.524807, 0.512861},
  { 1.12202, 0.60256, 0.660693,0.47863, 14.1254,   0.57544,  0.501187, 0.630957, 0.74131,  0.776247,  0.707946, 0.524807, 0.812831, 0.831764, 0.489779, 1.02329,  0.891251,  0.794328,  0.891251, 1},
  { 0.954993,1.41254, 1.1749,  1.23027,  0.57544,  1.86209,  1.47911,  0.794328, 1.31826,  0.645654,  0.691831, 1.41254,  0.794328, 0.549541, 0.954993, 1.04713,  1,         0.537032,  0.676083, 0.707946},
  { 1,       1.09648, 1.23027, 1.86209,  0.501187, 1.47911,  2.29087,  0.831764, 1.09648,  0.537032,  0.524807, 1.31826,  0.630957, 0.40738,  0.891251, 1.04713,  0.977237,  0.371535,  0.537032, 0.645654},
  { 1.12202, 0.794328,1.09648, 1.02329,  0.630957, 0.794328, 0.831764, 4.57088,  0.724436, 0.354813,  0.363078, 0.776247, 0.446684, 0.301995, 0.691831, 1.09648,  0.776247,  0.398107,  0.398107, 0.467735},
  { 0.831764,1.14815, 1.31826, 1.09648,  0.74131,  1.31826,  1.09648,  0.724436, 3.98107,  0.60256,   0.645654, 1.14815,  0.74131,  0.977237, 0.776247, 0.954993, 0.933254,  0.831764,  1.65959,  0.630957},
  { 0.831764,0.57544, 0.524807,0.416869, 0.776247, 0.645654, 0.537032, 0.354813, 0.60256,  2.51189,   1.90546,  0.616595, 1.77828,  1.25893,  0.549541, 0.660693, 0.870964,  0.660693,  0.851138, 2.04174},
  { 0.758578,0.60256, 0.501187,0.398107, 0.707946, 0.691831, 0.524807, 0.363078, 0.645654, 1.90546,   2.51189,  0.616595, 1.90546,  1.58489,  0.588844, 0.616595, 0.74131,   0.851138,  1,        1.51356},
  { 0.912011,1.86209, 1.20226, 1.12202,  0.524807, 1.41254,  1.31826,  0.776247, 1.14815,  0.616595,  0.616595, 2.0893,   0.724436, 0.467735, 0.870964, 1.02329,  1.02329,   0.446684,  0.616595, 0.676083},
  { 0.851138,0.676083,0.60256, 0.501187, 0.812831, 0.794328, 0.630957, 0.446684, 0.74131,  1.77828,   1.90546,  0.724436, 2.69153,  1.44544,  0.57544,  0.724436, 0.870964,  0.794328,  0.954993, 1.44544},
  { 0.588844,0.47863, 0.489779,0.354813, 0.831764, 0.549541, 0.40738,  0.301995, 0.977237, 1.25893,   1.58489,  0.467735, 1.44544,  5.01187,  0.416869, 0.524807, 0.60256,   2.29087,   3.23594,  1.02329},
  { 1.07152, 0.812831,0.812831,0.851138, 0.489779, 0.954993, 0.891251, 0.691831, 0.776247, 0.549541,  0.588844, 0.870964, 0.57544,  0.416869, 5.7544,   1.09648,  1.02329,   0.316228,  0.489779, 0.660693},
  { 1.28825, 0.954993,1.23027, 1.12202,  1.02329,  1.04713,  1.04713,  1.09648,  0.954993, 0.660693,  0.616595, 1.02329,  0.724436, 0.524807, 1.09648,  1.65959,  1.41254,   0.467735,  0.645654, 0.794328},
  { 1.14815, 0.954993,1.12202, 1,        0.891251, 1,        0.977237, 0.776247, 0.933254, 0.870964,  0.74131,  1.02329,  0.870964, 0.60256,  1.02329,  1.41254,  1.77828,   0.446684,  0.645654, 1},
  { 0.436516,0.691831,0.436516,0.301995, 0.794328, 0.537032, 0.371535, 0.398107, 0.831764, 0.660693,  0.851138, 0.446684, 0.794328, 2.29087,  0.316228, 0.467735, 0.446684, 26.3027,    2.5704,   0.549541},
  { 0.60256, 0.660693,0.724436,0.524807, 0.891251, 0.676083, 0.537032, 0.398107, 1.65959,  0.851138,  1,        0.616595, 0.954993, 3.23594,  0.489779, 0.645654, 0.645654,  2.5704,    6.0256,   0.776247},
  { 1.02329, 0.630957,0.60256, 0.512861, 1,        0.707946, 0.645654, 0.467735, 0.630957, 2.04174,   1.51356,  0.676083, 1.44544,  1.02329,  0.660693, 0.794328, 1,         0.549541,  0.776247, 2.18776} };


//------- constructor & destructor -------//
PROFILE::PROFILE(void)
{
	failure=1;       //default: success
	WantBlast=0;     //load PSP and PSM
	WantOriginal=0;  //load original
	length=0;        //inner length=0;
}
PROFILE::~PROFILE(void)
{
	if(length!=0)Profile_Delete_Matrix(length);
}

//---------- create & delete --------//
void PROFILE::Profile_Create_Matrix(int length)
{
	residue=new short[length];
	if(WantBlast==1)
	{
		NewArray2D_(&PSP,length,20);
		NewArray2D_(&PSM,length,20);
	}
	NewArray2D_(&ProfHMM,length,10);
	if(WantOriginal==1)NewArray2D_(&ProfHMM_original,length,10);
	NewArray2D_(&EmissionScore,length,20);
	if(WantOriginal==1)NewArray2D_(&EmissionScore_original,length,20);
	NewArray2D_(&EmissionProb,length,20);
	if(WantOriginal==1)NewArray2D_(&EmissionProb_original,length,20);
}
void PROFILE::Profile_Delete_Matrix(int length)
{
	delete [] residue;
	if(WantBlast==1)
	{
		DeleteArray2D_(&PSP,length);
		DeleteArray2D_(&PSM,length);
	}
	DeleteArray2D_(&ProfHMM,length);
	if(WantOriginal==1)DeleteArray2D_(&ProfHMM_original,length);
	DeleteArray2D_(&EmissionScore,length);
	if(WantOriginal==1)DeleteArray2D_(&EmissionScore_original,length);
	DeleteArray2D_(&EmissionProb,length);
	if(WantOriginal==1)DeleteArray2D_(&EmissionProb_original,length);
}

//----- profile related process -------//
void PROFILE::Process_PSM(ifstream &fin,string &filename,int StartPos)
{
	//--- init check ----//start
	if(WantBlast==0)return;
	//--- init check ----//over

	string wbuf;
	string temp;
	//skip first
	if(!getline(fin,wbuf,'\n'))
	{
		fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSM skip first] \n",filename.c_str());
		failure=0;
		return;
	}
	//get remaining
	for(int i=0;i<length;i++)
	{
		//get first
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSM line %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(wbuf.length()<20)
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSM line %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		//process string
		istringstream www(wbuf);
		for(int j=0;j<20;j++)
		{
			if(!(www>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSM line %d col %d] \n",filename.c_str(),i,j);
				failure=0;
				return;
			}
			float wsval=atof(temp.c_str());
			PSM[i-StartPos][j]= 0.01*wsval;
		}
	}
}
void PROFILE::Process_PSP(ifstream &fin,string &filename,int StartPos)
{
	//--- init check ----//start
	if(WantBlast==0)return;
	//--- init check ----//over

	string wbuf;
	string temp;
	//skip first
	if(!getline(fin,wbuf,'\n'))
	{
		fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSP skip first] \n",filename.c_str());
		failure=0;
		return;
	}
	//get remaining
	for(int i=0;i<length;i++)
	{
		//get string
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSP line %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(wbuf.length()<40)
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSP line %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		//process string
		istringstream www(wbuf);
		for(int j=0;j<22;j++)  //skip header
		{
			if(!(www>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSP line %d skip %d] \n",filename.c_str(),i,j);
				failure=0;
				return;
			}
		}
		for(int j=0;j<20;j++)  //process PSP
		{
			if(!(www>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [PSP line %d col %d] \n",filename.c_str(),i,j);
				failure=0;
				return;
			}
			float wsval=atof(temp.c_str());
			PSP[i-StartPos][j]= wsval;
		}
		//final check
		float wssum=0;
		for(int j=0;j<20;j++)wssum+=PSP[i-StartPos][j];
		if(wssum<80) //this is mainly due to NEFF=1
		{
			if(wssum!=0)fprintf(stderr,"WARNING %s -> PSP line %d low normalized %f \r",filename.c_str(),i,wssum);
			for(int j=0;j<20;j++)PSP[i-StartPos][j]=0; //this bug is discovered by Jinbo
			PSP[i-StartPos][residue[i-StartPos]]=100;
		}
		if(wssum>110) //this is wired
		{
			fprintf(stderr,"WARNING %s -> PSP line %d high normalized %f \n",filename.c_str(),i,wssum);
			for(int j=0;j<20;j++)PSP[i-StartPos][j]=0; //this bug is discovered by Jinbo
			PSP[i-StartPos][residue[i-StartPos]]=100;
		}
	}
}

//----------- load HHpred generated features -----------//
void PROFILE::Process_HMM(ifstream &fin,string &filename,int StartPos)
{
	string wbuf;
	string temp;
	//skip first
	for(int i=0;i<5;i++)
	{
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM skip first] \n",filename.c_str());
			failure=0;
			return;
		}
	}
	//process string
	for(int i=0;i<length;i++)
	{
		//get string 1
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(wbuf.length()<20)
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		istringstream _em(wbuf);
		for(int j=0;j<2;j++) //skip header
		{
			if(!(_em>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d header %d] \n",filename.c_str(),i,j);
				failure=0;
				return;
			}
		}
		for(int j=0;j<20;j++) //get remaining
		{
			if(!(_em>>temp))
			{
				fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line1 %d col %d] \n",filename.c_str(),i,j);
				failure=0;
				return;
			}
			if(temp[0]=='*')
			{
				EmissionScore[i-StartPos][j]=-99999;
				if(WantOriginal==1)EmissionScore_original[i-StartPos][j]=-99999;
			}
			else
			{
				EmissionScore[i-StartPos][j]=-atof(temp.c_str());
				if(WantOriginal==1)EmissionScore_original[i-StartPos][j] = EmissionScore[i-StartPos][j];
			}
			EmissionProb[i-StartPos][j]=pow(2.0,1.0*EmissionScore[i-StartPos][j]/1000.);
			if(WantOriginal==1)EmissionProb_original[i-StartPos][j] = EmissionProb[i-StartPos][j];
		}
		//check
		float wssum=0;
		for(int j=0;j<20;j++)wssum+=EmissionProb[i-StartPos][j];
		if(wssum<0.9) //this is mainly due to NEFF=1
		{
			for(int j=0;j<20;j++)
			{
				EmissionProb[i-StartPos][j]=0; //this bug is discovered by Jinbo
				if(WantOriginal==1)EmissionProb_original[i-StartPos][j]=0;
			}
			EmissionProb[i-StartPos][residue[i-StartPos]]=1;
			if(WantOriginal==1)EmissionProb_original[i-StartPos][residue[i-StartPos]]=1;
		}
		if(wssum>1.1) //this is wired
		{
			for(int j=0;j<20;j++)
			{
				EmissionProb[i-StartPos][j]=0; //this bug is discovered by Jinbo
				if(WantOriginal==1)EmissionProb_original[i-StartPos][j]=0;
			}
			EmissionProb[i-StartPos][residue[i-StartPos]]=1;
			if(WantOriginal==1)EmissionProb_original[i-StartPos][residue[i-StartPos]]=1;
		}

		//get string 2
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(wbuf.length()<7)
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HMM line2 %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		stringstream _sin(wbuf);
		//[1] M_M
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_M] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][M_M]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][M_M] = 0;
		}
		else
		{
			ProfHMM[i-StartPos][M_M]=exp(-atoi(temp.c_str())/1000.0*0.6931);  //transformed probability
			if(WantOriginal==1)ProfHMM_original[i-StartPos][M_M] = atoi(temp.c_str())/1000.0;    //original value of HHM file
		}
		//[2] M_I
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_I] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][M_I]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][M_I]=0;
		}
		else
		{
			ProfHMM[i-StartPos][M_I]=exp(-atof(temp.c_str())/1000.0*0.6931);
			if(WantOriginal==1)ProfHMM_original[i-StartPos][M_I] = atof(temp.c_str())/1000.0;
		}
		//[3] M_D
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_D] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][M_D]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][M_D] = ProfHMM[i-StartPos][M_D];
		}
		else
		{
			ProfHMM[i-StartPos][M_D]=exp(-atof(temp.c_str())/1000.0*0.6931);
			if(WantOriginal==1)ProfHMM_original[i-StartPos][M_D] = atof(temp.c_str())/1000.0;
		}
		//[4] I_M
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans I_M] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][I_M]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][I_M] = atof(temp.c_str())/1000.0;
		}
		else
		{
			ProfHMM[i-StartPos][I_M]=exp(-atoi(temp.c_str())/1000.0*0.6931);
			if(WantOriginal==1)ProfHMM_original[i-StartPos][I_M] = atof(temp.c_str())/1000.0;
		}
		//[5] I_I
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans I_I] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][I_I]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][I_I] = ProfHMM[i-StartPos][I_I];
		}
		else
		{
			ProfHMM[i-StartPos][I_I]=exp(-atof(temp.c_str())/1000.0*0.6931);
			if(WantOriginal==1)ProfHMM_original[i-StartPos][I_I] = atof(temp.c_str())/1000.0;
		}
		//[6] D_M
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans D_M] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][D_M]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][D_M]=0;
		}
		else
		{
			ProfHMM[i-StartPos][D_M]=exp(-atof(temp.c_str())/1000.0*0.6931);
			if(WantOriginal==1)ProfHMM_original[i-StartPos][D_M] = atof(temp.c_str())/1000.0;
		}
		//[7] D_D
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans D_D] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		if(temp[0]=='*')
		{
			ProfHMM[i-StartPos][D_D]=0;
			if(WantOriginal==1)ProfHMM_original[i-StartPos][D_D]=0;
		}
		else
		{
			ProfHMM[i-StartPos][D_D]=exp(-atof(temp.c_str())/1000.0*0.6931);
			if(WantOriginal==1)ProfHMM_original[i-StartPos][D_D] = atof(temp.c_str())/1000.0;
		}
		//final process
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans M_NEFF] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		ProfHMM[i-StartPos][_NEFF]=atof(temp.c_str());
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans I_NEFF] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		ProfHMM[i-StartPos][I_NEFF]=atof(temp.c_str());
		if(!(_sin>>temp))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line2 %d trans D_NEFF] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		ProfHMM[i-StartPos][D_NEFF]=atof(temp.c_str());
		//final
		ProfHMM[i-StartPos][_NEFF]/=1000;
		ProfHMM[i-StartPos][I_NEFF]/=1000;
		ProfHMM[i-StartPos][D_NEFF]/=1000;
		double rm; //0.6*ProfHMM[i][M_M] + 0.2*ProfHMM[i][M_I] + 0.2*ProfHMM[i][M_D];
		rm = 0.1;
		ProfHMM[i-StartPos][M_M] = (ProfHMM[i-StartPos][_NEFF]*ProfHMM[i-StartPos][M_M] + rm*0.6)/(rm+ProfHMM[i-StartPos][_NEFF]);
		ProfHMM[i-StartPos][M_I] = (ProfHMM[i-StartPos][_NEFF]*ProfHMM[i-StartPos][M_I] + rm*0.2)/(rm+ProfHMM[i-StartPos][_NEFF]);
		ProfHMM[i-StartPos][M_D] = (ProfHMM[i-StartPos][_NEFF]*ProfHMM[i-StartPos][M_D] + rm*0.2)/(rm+ProfHMM[i-StartPos][_NEFF]);
		double ri;//0.75*ProfHMM[i][I_I] + 0.25*ProfHMM[i][I_M];
		ri = 0.1;
		ProfHMM[i-StartPos][I_I] = (ProfHMM[i-StartPos][I_NEFF]*ProfHMM[i-StartPos][I_I] + ri*0.75)/(ri+ProfHMM[i-StartPos][I_NEFF]);
		ProfHMM[i-StartPos][I_M] = (ProfHMM[i-StartPos][I_NEFF]*ProfHMM[i-StartPos][I_M] + ri*0.25)/(ri+ProfHMM[i-StartPos][I_NEFF]);
		double rd;// 0.75*ProfHMM[i][D_D] + 0.25*ProfHMM[i][D_M];
		rd = 0.1;
		ProfHMM[i-StartPos][D_D] = (ProfHMM[i-StartPos][D_NEFF]*ProfHMM[i-StartPos][D_D] + rd*0.75)/(rd+ProfHMM[i-StartPos][D_NEFF]);
		ProfHMM[i-StartPos][D_M] = (ProfHMM[i-StartPos][D_NEFF]*ProfHMM[i-StartPos][D_M] + rd*0.25)/(rd+ProfHMM[i-StartPos][D_NEFF]);
		double ws_tmp_neff=ProfHMM[i-StartPos][_NEFF]-1;

		//get string 3
		if(!getline(fin,wbuf,'\n'))
		{
			fprintf(stderr,"FORMAT BAD AT FEATURE FILE %s -> [HHM line3 %d] \n",filename.c_str(),i);
			failure=0;
			return;
		}
		for(int j=0;j<20;j++)
		{
			double g=0;
			int ind = AA1Coding[j];
			for(int k=0;k<20;k++)g+=EmissionProb[i-StartPos][k]*gonnet[AA1Coding[k]][ind]*pow(2.0,-1.0*HMMNull[j]/1000.0);
			EmissionScore[i-StartPos][j] = (ws_tmp_neff*EmissionProb[i-StartPos][j]+g*10)/(ws_tmp_neff+10);
		}
		for(int j=0;j<20;j++)
		{
			EmissionProb[i-StartPos][j] = EmissionScore[i-StartPos][j];
			EmissionScore[i-StartPos][j] = log(EmissionProb[i-StartPos][j])/log(2.);
		}
	}
}
