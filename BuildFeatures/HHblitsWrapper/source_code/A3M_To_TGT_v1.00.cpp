#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
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

//----------------------//
//-- Getline_Ending ---//
//--------------------//
void getline_end(string &input,char kill)
{
	int len=(int)input.length();
	if(input[len-1]==kill)input=input.substr(0,len-1);
}


// amino acid related data structure
const int maxAcc[21]={101,242,171,161,127,190,195,73,193,174,182,200,192,214,127,122,144,259,235,148,-1};
const int AA2SUB[26]={0,20,1,2,3,4,5,6,7,20,8,9,10,11,20,12,13,14,15,16,20,17,18,20,19,20};
const char AA1Coding[21]={0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18,20};
string AA3Coding[26]={"ALA","XXX","CYS","ASP","GLU","PHE","GLY","HIS","ILE","XXX","LYS","LEU","MET","ASN","XXX","PRO","GLN","ARG","SER","THR","XXX","VAL","TRP","XXX","TYR","XXX"};

// acc definition
int BURIED = 0;
int INTERM = 1;
int EXPOSE = 2;
float BuryCut = 10.25;
float ExpCut = 42.9;

// maxlength
const int MaxLen = 10000;
const int BufLen = 100000;
int fullLen;

// data structure related
float PSP[MaxLen][40];
float SSE[MaxLen][3];
float ACC[MaxLen][3];
char SSEconf[MaxLen];
char ACCconf[MaxLen];
char fullSSE[MaxLen];
char fullACC[MaxLen];
int HHM[MaxLen][40];

// hhm related definition
const int M_M = 0;
const int M_I = 1;
const int M_D = 2;
const int I_M = 3;
const int I_I = 4;
const int D_M = 5;
const int D_D = 6;
const int _NEFF = 7;
const int I_NEFF = 8;
const int D_NEFF = 9;
float ProfHMM[MaxLen][20], EmissionScore[MaxLen][21];
float NEFF;



//---- kill space ----//
void Kill_Invalid(string &in,string &out)
{
	int i;
	int len=(int)in.length();
	out="";
	for(i=0;i<len;i++)
	{
		if(in[i]!=' ' && in[i]!='-' )out=out+in[i];
	}
}


//---- readme ----//
void Usage(char *arg) 
{
	printf("Version: 1.00 \n");
	printf("Usage: %s -i input_fasta_file -I input_a3m_file [ -o output_tgt_file ] [ -t temporary_root ] \n",arg);
}
//---- parameter editor ----//
static option long_options[] =
{
	{"input",   required_argument, NULL, 'i'},
	{"Input",   required_argument, NULL, 'I'},
	{"output",  no_argument,       NULL, 'o'},
	{"temp",    no_argument,       NULL, 't'},
	{0, 0, 0, 0}
};

//------------ main -----------//
int main(int argc, char** argv)
{

	if(argc<2)
	{
		Usage(argv[0]);
		exit(-1);
	}
	string input_fasta="";
	string input_a3m="";
	string output_tgt="";
	string temporary_root="";

	extern char* optarg;
	char c = 0;
	int option_index=0;
	while ((c = getopt_long(argc, argv, "i:I:o:t:",long_options,&option_index)) != EOF) 
	{
		switch (c) 
		{
			case 'i':
				input_fasta = optarg;
				break;
			case 'I':
				input_a3m = optarg;
				break;
			case 'o':
				output_tgt = optarg;
				break;
			case 't':
				temporary_root = optarg;
				break;
			default:
				Usage(argv[0]);
				exit(-1);
		}
	}
	if(input_fasta=="" || input_a3m=="" )
	{
		Usage(argv[0]);
		exit(-1);
	}

//---- mkdir tmp dir ---//__121130__//
	if(temporary_root=="")temporary_root="tmp/";
	string mkdir_cmd = "mkdir -p "+temporary_root;
	system(mkdir_cmd.c_str());
//--- over ---//


	string user_command=argv[0];
	for(int i=1;i<argc;i++)user_command = user_command+" " + argv[i];
	//get base name
	string targetName;
	getBaseName(input_fasta,targetName,'/','.');
	if(output_tgt=="")output_tgt=targetName+".tgt";
	string sequence_ = "";
	string sequence = "";
	char buf[BufLen];

// generate SEQ file
	string seqfile0=input_fasta;
	ifstream seq0_in(seqfile0.c_str());
	if(!seq0_in.is_open())
	{
		cerr << "Cannot open input_fasta_file " << seqfile0 << endl;
		exit(1);
	}
	while(seq0_in.getline(buf,BufLen))
	{
		string wstmp=buf;
		getline_end(wstmp,0x0D);
		if(buf[0]!='>') sequence_ = sequence_ + wstmp;
	}
	Kill_Invalid(sequence_,sequence);
	fullLen = sequence.length();
	seq0_in.close();

// fix unknown amino acids other than ARNDCQEGHILKMFPSTWYVX
	for(int i=0;i<fullLen;i++)
	{
		char c=sequence[i];
		if(c<'A' || c>'Z')
		{
			sequence[i]='X';
			continue;
		}
		if(AA2SUB[c-'A']==20)
		{
			sequence[i]='X';
			continue;
		}
	}

//===== copy input_fasta to tmp ======//__121130__//
//	string copy_cmd = "cp "+input_fasta+" tmp/"+targetName+".seq";
//	system(copy_cmd.c_str());
	string seqfile1=temporary_root+"/"+targetName+".seq";
	FILE *fpseq=fopen(seqfile1.c_str(),"wb");
	fprintf(fpseq,">%s\n",targetName.c_str());
	fprintf(fpseq,"%s\n",sequence.c_str());
	fclose(fpseq);

//===== input MTX file =======//
	string mtxfile1=temporary_root+"/"+targetName+".mtx";

// command to generate features 
	string genMAIN_cmd = "util/genMAIN.sh " + input_a3m + " "+temporary_root+"/";
	string genSS3_cmd = "util/genSS3.sh " + mtxfile1 + " "+temporary_root+"/";
	string genDIS_cmd = "util/genDIS.sh " + mtxfile1 + " "+temporary_root+"/";
	string genSS8_cmd = "util/genSS8.sh " + targetName + " "+temporary_root+"/";
	string genACC_cmd = "util/genACC.sh " + targetName + " "+temporary_root+"/";
// generate MAIN file
	{
		system(genMAIN_cmd.c_str());
	}
// generate SS3 file
	{
		system(genSS3_cmd.c_str());
	}
// generate DIS file
	{
		system(genDIS_cmd.c_str());
	}
// generate SS8 file
	{
		system(genSS8_cmd.c_str());
	}
// generate ACC file
	{
		system(genACC_cmd.c_str());
	}
/////////////////////////////////// all feature generated


	int length = sequence.length();

// get current time
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);

// read-in feature files
	// read SS2 file
	string ss2file = temporary_root+"/" + targetName + ".ss2";
	ifstream ss2_in(ss2file.c_str());	
	if(!ss2_in.is_open())
	{
		cerr << "Cannot open ss2file " << ss2file << endl;
		exit(1);
	}
	ss2_in.getline(buf,BufLen);			
	ss2_in.getline(buf,BufLen);			
	for(int index=0;index<fullLen;index++)
	{
		if(!(ss2_in.getline(buf,BufLen)))
		{
			cerr << "ss2file format bad " << endl;
			exit(1);
		}
		sscanf(buf+9, "%f%f%f", &SSE[index][0], &SSE[index][1], &SSE[index][2]);
	}
	ss2_in.close();
	//assign sse
	for(int i=0;i<fullLen;i++)
	{
		if(SSE[i][0]>SSE[i][1] && SSE[i][0] > SSE[i][2])
		{
			fullSSE[i] = 'L';
			int c = (int)(10*SSE[i][0]);
			if(c<0)c=0;
			else if(c>9)c=9;
			SSEconf[i]=c+'0';
		}
		else if(SSE[i][1]>SSE[i][0] && SSE[i][1]>SSE[i][2])
		{
			fullSSE[i] = 'H';
			int c = (int)(10*SSE[i][1]);
			if(c<0)c=0;
			else if(c>9)c=9;
			SSEconf[i]=c+'0';
		}
		else
		{
			fullSSE[i] = 'E';
			int c = (int)(10*SSE[i][2]);
			if(c<0)c=0;
			else if(c>9)c=9;
			SSEconf[i]=c+'0';
		}
	}
	fullSSE[fullLen]='\0';
	SSEconf[fullLen]='\0';
	string sse_sequence = fullSSE;
	string sse_confidence = SSEconf;


	// read ACC file
	string accfile = temporary_root+"/" + targetName + ".acc";
	ifstream acc_in(accfile.c_str());
	if(!acc_in.is_open())
	{
		cerr << "Cannot open accfile " << accfile << endl;
		exit(1);
	}
	for(int index=0;index<fullLen;index++)
	{
		if(!(acc_in.getline(buf,BufLen)))
		{
			cerr << "accfile format bad " << endl;
			exit(1);
		}
		sscanf(buf, "%f%f%f", &ACC[index][0], &ACC[index][1], &ACC[index][2]);
	}
	acc_in.close();
	//assign sse
	for(int i=0;i<fullLen;i++)
	{
		if(ACC[i][0]>ACC[i][1] && ACC[i][0] > ACC[i][2])
		{
			fullACC[i] = 'B';
			int c = (int)(10*ACC[i][0]);
			if(c<0)c=0;
			else if(c>9)c=9;
			ACCconf[i]=c+'0';
		}
		else if(ACC[i][1]>ACC[i][0] && ACC[i][1]>ACC[i][2])
		{
			fullACC[i] = 'M';
			int c = (int)(10*ACC[i][1]);
			if(c<0)c=0;
			else if(c>9)c=9;
			ACCconf[i]=c+'0';
		}
		else
		{
			fullACC[i] = 'E';
			int c = (int)(10*ACC[i][2]);
			if(c<0)c=0;
			else if(c>9)c=9;
			ACCconf[i]=c+'0';
		}
	}
	fullACC[fullLen]='\0';
	ACCconf[fullLen]='\0';
	string acc_sequence = fullACC;
	string acc_confidence = ACCconf;


// read HHM / NEFF
	string hhmfile = temporary_root+"/" + targetName + ".hhm";
	ifstream hhm_in(hhmfile.c_str());
	if(!hhm_in.is_open())
	{
		cerr << "Cannot open hhmfile " << hhmfile << endl;
		exit(1);
	}
	while(hhm_in.getline(buf,BufLen) && buf[0]!='#')
	{
		if(buf[0]=='N'&&buf[1]=='E'&&buf[2]=='F'&&buf[3]=='F')
		{
			sscanf(buf+4,"%f",&NEFF);
		}
	}

// combine all features together

// TEMP file
//	string out_filename = "TARGET/" + targetName + ".tgt";
	string out_filename=output_tgt;
	ofstream fout(out_filename.c_str());
// output all features
// index AA isMissing SS coreIndex ACC Xca Yca Zca
	fout << "Version 1.1 " << "-> user command = " << user_command << endl;
	fout << "Sequence Name  = " << targetName << endl;
	fout << "Length   = " << fullLen << endl;
	fout << "Sequence = " << sequence << endl;
	fout << "SSEseq   = " << sse_sequence << endl;
	fout << "SSEconf  = " << sse_confidence << endl;
	fout << "ACCseq   = " << acc_sequence << endl;
	fout << "ACCconf  = " << acc_confidence << endl;
	fout << "NEFF = " << NEFF << endl;
	fout << "EVD = " << endl;
	fout << "Date = " << tm.tm_year + 1900 << "-" << tm.tm_mon+1 << "-" << tm.tm_mday << "  " << tm.tm_hour << ":" << tm.tm_hour << ":" << tm.tm_sec << endl;	


// PSM info
	fout << "//////////// Original PSM file (generated by BLAST makemat -P targetname )" << endl;
	fout << "//////////// " << tm.tm_year + 1900 << "-" << tm.tm_mon+1 << "-" << tm.tm_mday << "  " << tm.tm_hour << ":" << tm.tm_hour << ":" << tm.tm_sec << endl;
	// read MTX
	string mtxfile = temporary_root+"/" + targetName + ".mtx";
	ifstream mtx_in(mtxfile.c_str());
	if(!mtx_in.is_open())
	{
		cerr << "Cannot open mtxfile " << mtxfile << endl;
		exit(1);
	}
	for(int i=0;i<14;i++)
		mtx_in.getline(buf,BufLen);
	for(int i=0;i<fullLen;i++)
	{
		char aa_scores[23][20];
		for(int j=0;j<23;j++)
		{
			mtx_in >> aa_scores[j];
		}
		sprintf(buf, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
			aa_scores[1], aa_scores[16], aa_scores[13], aa_scores[4],
			aa_scores[3], aa_scores[15], aa_scores[5], aa_scores[7],
			aa_scores[8], aa_scores[9], aa_scores[11], aa_scores[10],
			aa_scores[12], aa_scores[6], aa_scores[14], aa_scores[17],
			aa_scores[18], aa_scores[20], aa_scores[22], aa_scores[19]);
		fout << buf << endl;
		for(int j=0;j<5;j++)
			mtx_in >> buf;
	}
	mtx_in.close();
	fout << endl << endl << endl;


// PSP info
	fout << "//////////// Original PSP file (generated by BLAST -b 0 -j 5 -h 0.001 -d nr90)" << endl;
	fout << "//////////// " << tm.tm_year + 1900 << "-" << tm.tm_mon+1 << "-" << tm.tm_mday << "  " << tm.tm_hour << ":" << tm.tm_hour << ":" << tm.tm_sec << endl;
	// read PSP
	string pspfile = temporary_root+"/" + targetName + ".psp";
	ifstream psp_in(pspfile.c_str());	
	if(!psp_in.is_open())
	{
		cerr << "Cannot open pspfile " << pspfile << endl;
		exit(1);
	}		
	psp_in.getline(buf,BufLen);
	psp_in.getline(buf,BufLen);
	psp_in.getline(buf,BufLen);
	for(int index=0;index<fullLen;index++)
	{
		if(!(psp_in.getline(buf,BufLen)))
		{
			cerr << "pspfile format bad " << endl;
			exit(1);
		}
		fout << buf << endl;
	}	
	psp_in.close();
	fout << endl << endl << endl;


// disorder info
	fout << "//////////// Original DIS file (generated by DISOPRED version 2)" << endl;
	fout << "//////////// " << tm.tm_year + 1900 << "-" << tm.tm_mon+1 << "-" << tm.tm_mday << "  " << tm.tm_hour << ":" << tm.tm_hour << ":" << tm.tm_sec << endl;
	// read PSP
	string disofile = temporary_root+"/" + targetName + ".diso";
	ifstream diso_in(disofile.c_str());
	if(!diso_in.is_open())
	{
		cerr << "Cannot open disofile " << disofile << endl;
		exit(1);
	}
	diso_in.getline(buf,BufLen);
	diso_in.getline(buf,BufLen);
	diso_in.getline(buf,BufLen);
	diso_in.getline(buf,BufLen);
	diso_in.getline(buf,BufLen);
	for(int index=0;index<fullLen;index++)
	{
		if(!(diso_in.getline(buf,BufLen)))
		{
			cerr << "disofile format bad " << endl;
			exit(1);
		}
		fout << buf << endl;
	}
	diso_in.close();
	fout << endl << endl << endl;



// HMM info
	fout << "//////////// Original HHM file (generated by buildali2.pl + hhmake)" << endl;
	fout << "//////////// " << tm.tm_year + 1900 << "-" << tm.tm_mon+1 << "-" << tm.tm_mday << "  " << tm.tm_hour << ":" << tm.tm_hour << ":" << tm.tm_sec << endl;
	int hmm_count=0;
	for(int index=0;index<3*fullLen+4;index++)
	{
		if(!(hhm_in.getline(buf,BufLen)))
		{
			cerr << "hhmfile format bad " << endl;
			exit(1);
		}
		fout << buf << endl;
	}
	hhm_in.close();
	fout << endl << endl << endl;

// structural info
	fout << "//////////// Original SS3+SS8+ACC file (generated by PSIPRED + RaptorX-SS8 + RaptorX-ACC )" << endl;
	fout << "//////////// " << tm.tm_year + 1900 << "-" << tm.tm_mon+1 << "-" << tm.tm_mday << "  " << tm.tm_hour << ":" << tm.tm_hour << ":" << tm.tm_sec << endl;
	fout << "SS3:   H     E     L  | SS8:   H     G     I     E     B     T     S     L  | ACC:  Bury   Medium  Exposed" << endl;
	fout.close();	
	string fin_combine="";
	fin_combine = fin_combine + "util/combine_file.sh "+temporary_root+"/"+targetName+".ss2 "+temporary_root+"/"+targetName+".ss8 "+temporary_root+"/"+targetName+".acc "+output_tgt;
	system(fin_combine.c_str());


// delete junk
	string rm_command;
	//rm_command = "";
	//rm_command = rm_command + "rm -f tmp/"+targetName+"_sse*";
	//system(rm_command.c_str());
	//rm_command = "";
	//rm_command = rm_command + "rm -f tmp/"+targetName+".mn";
	//system(rm_command.c_str());
	rm_command = "";
	rm_command = rm_command + "rm -f error.log";
	system(rm_command.c_str());

// ---------- build ss3 and acc using SS8_To_SS3 and ACC_ReFormat ------------- //__2014_10_30__//
	string nnn_command;
	nnn_command="";
	nnn_command = nnn_command + "util/SS8_To_SS3 "+temporary_root+"/"+targetName+".ss8 > "+temporary_root+"/"+targetName+".ss3";
	system(nnn_command.c_str());
	nnn_command="";
	nnn_command = nnn_command + "util/ACC_ReFormat "+temporary_root+"/"+targetName+".acc "+temporary_root+"/"+targetName+".seq > "+temporary_root+"/"+targetName+".acc_";
	system(nnn_command.c_str());
	nnn_command="";
	nnn_command = nnn_command + "mv "+temporary_root+"/"+targetName+".acc_ "+temporary_root+"/"+targetName+".acc";
	system(nnn_command.c_str());
}

