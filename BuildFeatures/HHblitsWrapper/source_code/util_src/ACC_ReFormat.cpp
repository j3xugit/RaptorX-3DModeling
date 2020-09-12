#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;


//------------- fasta read ------------//
void FASTA_Read(string &infile,string &outfile)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)return;
	//process
	if(!getline(fin,buf,'\n'))return;
	if(!getline(fin,buf,'\n'))return;
	outfile=buf;
}

//------------- ACC reformat ----------//
void ACC_ReFormat(string &infile,string &sequence)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)return;
	//process
	int count=0;
	printf("#RaptorX-ACC: three-class solvent accessibility prediction results \n");
	printf("#probabilities are in the order of B (Bury, pACC: 0-10), M (Medium, pACC: 10-42) and E (Exposed, pACC: 42-100), \n");
	printf("#where pACC is the normalized solvent accessibility value calculated by DSSP \n");
	printf("\n");
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//Bury
		double bury=0;
		www>>bury;
		//Medium
		double medium=0;
		www>>medium;
		//Exposed
		double exposed=0;
		www>>exposed;
		//get max
		char state='M';
		if(bury > medium && bury > exposed)state='B';
		if(medium > bury && medium > exposed)state='M';
		if(exposed > bury && exposed > medium) state='E';
		//printf
		count++;
		printf(" %d %c %c %5.3f %5.3f %5.3f \n",count,sequence[count-1],state,bury,medium,exposed);
	}
}

//-------- main --------//
int main(int argc,char **argv)
{
	//---- ACC_ReFormat -------//
	{
		if(argc<3)
		{
			fprintf(stderr,"ACC_ReFormat <acc_file> <seq_file> \n");
			exit(-1);
		}
		string acc_file=argv[1];
		string seq_file=argv[2];
		string sequence;
		FASTA_Read(seq_file,sequence);
		ACC_ReFormat(acc_file,sequence);
		exit(0);
	}
}

