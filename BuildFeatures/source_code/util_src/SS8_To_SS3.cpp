#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;



//------------- SS8 -> SS3 ----------//
void WS_SS8_To_SS3(string &infile)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)return;
	//skip 3
	for(int i=0;i<3;i++)if(!getline(fin,buf,'\n'))return;
	printf("#RaptorX-SS3: three-class secondary structure prediction results \n");
	printf("#probabilities are in the order of H E C, the 3 secondary structure types used in DSSP \n");
	printf("\n");
	//process
	string pos;
	string ami;
	double tmp;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//skip 3
		www>>pos;
		www>>ami;
		www>>temp;
		//helix
		double helix=0;
		for(int i=0;i<3;i++)
		{
			www>>tmp;
			helix+=tmp;
		}
		//sheet
		double sheet=0;
		for(int i=0;i<2;i++)
		{
			www>>tmp;
			sheet+=tmp;
		}
		//coil
		double coil=0;
		for(int i=0;i<3;i++)
		{
			www>>tmp;
			coil+=tmp;
		}
		//get max
		char state='C';
		if(helix > sheet && helix > coil)state='H';
		if(sheet > helix && sheet > coil)state='E';
		if(coil > helix && coil > sheet) state='C';
		//printf
		printf(" %s %s %c %5.3f %5.3f %5.3f \n",pos.c_str(),ami.c_str(),state,helix,sheet,coil);
	}
}

//-------- main --------//
int main(int argc,char **argv)
{
	//---- SS8 to SS3 -------//
	{
		if(argc<2)
		{
			fprintf(stderr,"SS8_To_SS3 <ss8_file> \n");
			exit(-1);
		}
		string ss8_file=argv[1];
		WS_SS8_To_SS3(ss8_file);
		exit(0);
	}
}

