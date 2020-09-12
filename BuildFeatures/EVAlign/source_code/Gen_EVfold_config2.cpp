#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
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


//-------- create evfold list -------//
//-> template
/*
pipeline: protein_monomer

stages:
    - align

global:
    prefix: output/1pazA
    sequence_id: 1pazA
    sequence_file: 1pazA.fasta

    region:
    theta: 0.8
    cpu: 1

align:
    protocol: standard
    first_index: 1
    use_bitscores: True
    domain_threshold: 0.5
    sequence_threshold: 0.5
    iterations: 5
    database: uniref90
    compute_num_effective_seqs: False
    seqid_filter:
    minimum_sequence_coverage: 50
    minimum_column_coverage: 70
    extract_annotation: True
    cpu: 1
    nobias: False
    reuse_alignment: True
    checkpoints_hmm: False
    checkpoints_ali: False

databases:
     uniref90: /home/RaptorX/wstest/JackHmmer_Test/uniref90.fasta

tools:
     jackhmmer: /home/RaptorX/wstest/JackHmmer_Test/jackhmmer
     hhfilter: NULL
*/

//----------- generate config -------------//
void Generate_EVfold_config(string &tmpout, string &seqid, string &file, 
	int cpu, int bitsco, double thres, int iter, string &uniref, string &root)
{
	//--- determine using bitsco or not ----//
	string bitstr="";
	if(bitsco==1)bitstr="True";
	else bitstr="False";

	//--- print config file ----//
	printf("pipeline: protein_monomer\n");
	printf("stages:\n");
	printf("    - align\n");
	printf("global:\n");
	printf("    prefix: %s\n",tmpout.c_str());
	printf("    sequence_id: %s\n",seqid.c_str());
	printf("    sequence_file: %s\n",file.c_str());
	printf("    region:\n");
	printf("    theta: 0.8\n");
	printf("    cpu: %d\n",cpu);
	printf("align:\n");
	printf("    protocol: standard\n");
	printf("    first_index: 1\n");
	printf("    use_bitscores: %s\n",bitstr.c_str());
	printf("    domain_threshold: %.3e\n",thres);    //-> default: 0.5 for bitsco, 0.001 for evalue
	printf("    sequence_threshold: %.3e\n",thres);  //-> default: 0.5 for bitsco, 0.001 for evalue
	printf("    iterations: %d\n",iter);
	printf("    database: uniref90\n");
	printf("    compute_num_effective_seqs: False\n");
	printf("    seqid_filter:\n");
	printf("    minimum_sequence_coverage: 50\n");
	printf("    minimum_column_coverage: 70\n");
	printf("    extract_annotation: True\n");
	printf("    cpu: %d\n",cpu);
	printf("    nobias: False\n");
	printf("    reuse_alignment: True\n");
	printf("    checkpoints_hmm: False\n");
	printf("    checkpoints_ali: False\n");
	printf("databases:\n");
	//printf("     uniref90: %s/databases/%s.fasta\n",root.c_str(),uniref.c_str());
	printf("     uniref90: %s\n",uniref.c_str());
	printf("tools:\n");
	printf("     jackhmmer: %s/util/jackhmmer\n",root.c_str());
	printf("     hhfilter: NULL\n");
}


//------ main ------//
int main(int argc,char **argv)
{
	//----- Generate EVfold config ---//
	{
		if(argc<9)
		{
			fprintf(stderr,"./Gen_EVfold_config <tmpout> <seqfile> <cpu> <bitsco> <thres> \n");
			fprintf(stderr,"     <iter> <uniref> <root>\n");
			exit(-1);
		}
		//--- arguments ----//
		string tmpout=argv[1];
		string seqfile=argv[2];
		int cpu=atoi(argv[3]);
		int bitsco=atoi(argv[4]);
		double thres=atof(argv[5]);
		int iter=atoi(argv[6]);
		string uniref=argv[7];
		string root=argv[8];
		
		//--- process ----//
		string seqid;
		getBaseName(seqfile,seqid,'/','.');
		Generate_EVfold_config(tmpout,seqid,seqfile,cpu,bitsco,thres,iter,uniref,root);

		//--- exit ----//
		exit(0);
	}
}

