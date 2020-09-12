/* MetaPSICOV Alignment Stats Program - by David T. Jones June 2014 */

/* Copyright (C) 2014 University College London */

/* V1.01 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>
#include <vector>
using namespace std;

#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))

#define MAXSEQLEN 5000



const char *rescodes = "ARNDCQEGHILKMFPSTWYVXXX";

/* Betancourt and Thirumalai Contact Potential */
float contmat[20][20] = {
    {-0.20, 0.27, 0.24, 0.30,-0.26, 0.21, 0.43,-0.03, 0.21,-0.35,-0.37, 0.20,-0.23,-0.33, 0.07, 0.15, 0.00,-0.40,-0.15,-0.38 },
    { 0.27, 0.13, 0.02,-0.71, 0.32,-0.12,-0.75, 0.14, 0.04, 0.18, 0.09, 0.50, 0.17, 0.08,-0.02, 0.12, 0.00,-0.41,-0.37, 0.17 },
    { 0.24, 0.02,-0.04,-0.12, 0.28,-0.05,-0.01, 0.10, 0.10, 0.55, 0.36,-0.14, 0.32, 0.29, 0.13, 0.14, 0.00,-0.09, 0.01, 0.39 },
    { 0.30,-0.71,-0.12, 0.27, 0.38, 0.12, 0.40, 0.17,-0.22, 0.54, 0.62,-0.69, 0.62, 0.48, 0.25, 0.01, 0.00, 0.06,-0.07, 0.66 },
    {-0.26, 0.32, 0.28, 0.38,-1.34, 0.04, 0.46,-0.09,-0.19,-0.48,-0.50, 0.35,-0.49,-0.53,-0.18, 0.09, 0.00,-0.74,-0.16,-0.51 },
    { 0.21,-0.12,-0.05, 0.12, 0.04, 0.14, 0.10, 0.20, 0.22, 0.14, 0.08,-0.20,-0.01,-0.04,-0.05, 0.25, 0.00,-0.11,-0.18, 0.17 },
    { 0.43,-0.75,-0.01, 0.40, 0.46, 0.10, 0.45, 0.48,-0.11, 0.38, 0.37,-0.87, 0.24, 0.34, 0.26, 0.10, 0.00,-0.15,-0.16, 0.41 },
    {-0.03, 0.14, 0.10, 0.17,-0.09, 0.20, 0.48,-0.20, 0.23, 0.21, 0.14, 0.12, 0.08, 0.11,-0.01, 0.10, 0.00,-0.24,-0.04, 0.04 },
    { 0.21, 0.04, 0.10,-0.22,-0.19, 0.22,-0.11, 0.23,-0.33, 0.19, 0.10, 0.26,-0.17,-0.19,-0.05, 0.15, 0.00,-0.46,-0.21, 0.18 },
    {-0.35, 0.18, 0.55, 0.54,-0.48, 0.14, 0.38, 0.21, 0.19,-0.60,-0.79, 0.21,-0.60,-0.65, 0.05, 0.35, 0.00,-0.65,-0.33,-0.68 },
    {-0.37, 0.09, 0.36, 0.62,-0.50, 0.08, 0.37, 0.14, 0.10,-0.79,-0.81, 0.16,-0.68,-0.78,-0.08, 0.26, 0.00,-0.70,-0.44,-0.80 },
    { 0.20, 0.50,-0.14,-0.69, 0.35,-0.20,-0.87, 0.12, 0.26, 0.21, 0.16, 0.38, 0.22, 0.11, 0.12, 0.10, 0.00,-0.28,-0.40, 0.16 },
    {-0.23, 0.17, 0.32, 0.62,-0.49,-0.01, 0.24, 0.08,-0.17,-0.60,-0.68, 0.22,-0.56,-0.89,-0.16, 0.32, 0.00,-0.94,-0.51,-0.47 },
    {-0.33, 0.08, 0.29, 0.48,-0.53,-0.04, 0.34, 0.11,-0.19,-0.65,-0.78, 0.11,-0.89,-0.82,-0.19, 0.10, 0.00,-0.78,-0.49,-0.67 },
    { 0.07,-0.02, 0.13, 0.25,-0.18,-0.05, 0.26,-0.01,-0.05, 0.05,-0.08, 0.12,-0.16,-0.19,-0.07, 0.17, 0.00,-0.73,-0.40,-0.08 },
    { 0.15, 0.12, 0.14, 0.01, 0.09, 0.25, 0.10, 0.10, 0.15, 0.35, 0.26, 0.10, 0.32, 0.10, 0.17, 0.13, 0.00, 0.07, 0.07, 0.25 },
    { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
    {-0.40,-0.41,-0.09, 0.06,-0.74,-0.11,-0.15,-0.24,-0.46,-0.65,-0.70,-0.28,-0.94,-0.78,-0.73, 0.07, 0.00,-0.74,-0.55,-0.62 },
    {-0.15,-0.37, 0.01,-0.07,-0.16,-0.18,-0.16,-0.04,-0.21,-0.33,-0.44,-0.40,-0.51,-0.49,-0.40, 0.07, 0.00,-0.55,-0.27,-0.27 },
    {-0.38, 0.17, 0.39, 0.66,-0.51, 0.17, 0.41, 0.04, 0.18,-0.68,-0.80, 0.16,-0.47,-0.67,-0.08, 0.25, 0.00,-0.62,-0.27,-0.72 }
};


/* Dump a rude message to standard error and exit */
void fail(char *errstr)
{
	fprintf(stderr, "\n*** %s\n\n", errstr);
	exit(-1);
}

/* Convert AA letter to numeric code (0-21) */
int aanum(int ch)
{
	const static int aacvs[] =
	{
		999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
		21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
	};
	
	return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Calculate sequence weight */
int compute_similarity(int i,int j,char  **aln, double idthresh, int seqlen)
{
	int nthresh = idthresh * seqlen;
	for (int kk=0; nthresh > 0 && kk<seqlen; kk++)
		if (aln[i][kk] != aln[j][kk])
			nthresh--;
	return nthresh;
}

/* Calculate pal_local and mi */
double calculate_localMI(int i,int j,char  **aln, int nseqs, double wtsum, 
	vector <vector <double> > &pa, vector <float> &weight)
{
	double pab_local[21][21];

	/* Calculate pseudocount = 1 */
	for (int aa=0; aa<21; aa++)
		for (int bb=0; bb<21; bb++)
			pab_local[aa][bb] = 1.0 / 21.0;

	/* singleton */
	for (int kk=0; kk<nseqs; kk++)
	{
		int aa = aln[kk][i];
		int bb = aln[kk][j];
		if (aa < 21 && bb < 21)
			pab_local[aa][bb] += weight[kk];
	}

	/* pariwise */
	for (int aa=0; aa<21; aa++)
		for (int bb=0; bb<21; bb++)
		{
			pab_local[aa][bb] /= 21.0 + wtsum;
			pab_local[bb][aa] = pab_local[aa][bb];
		}

	/* final assignment */
	double mi=0;
	for (int aa=0; aa<20; aa++)
		for (int bb=0; bb<20; bb++)
			mi += pab_local[aa][bb] * log(pab_local[aa][bb]/pa[i][aa]/pa[j][bb]);

	//-- return ---//
	return mi;
}

/* Calculate potsum */
double calculate_Potsum(int i,int j,char  **aln, int nseqs, double wtsum,
	vector <float> &weight)
{
	//-- calculate --//
	double potsum=0;
	for (int kk=0; kk<nseqs; kk++)
	{
		int aa = aln[kk][i];
		int bb = aln[kk][j];
		if (aa < 20 && bb < 20 && aa>=0 && bb>=0)
			potsum += weight[kk] * contmat[aa][bb];
	}
	potsum /= wtsum;

	//-- return ---//
	return potsum;
}


/* ------- main ---------*/
int main(int argc, char **argv)
{
	int a, b, c, d, i, j, k, l, n, nseqs, seqlen;
	double aacomp[21], idthresh=0.38;
	char buf[MAXSEQLEN], seq[MAXSEQLEN];
	char **aln;
	FILE *ifp, *singofp, *pairofp;

	/* ----- parse arguments -----*/
	if (argc != 4)
		fail("Usage: alnstats alnfile singoutfile pairoutfile");

	ifp = fopen(argv[1], "r");
	if (!ifp)
		fail("Unable to open alignment file!");


	/* ----- PART I: initialize data structure ------*/ /*start*/
	//-> get nseqs
	for (nseqs=0;; nseqs++)
		if (!fgets(seq, MAXSEQLEN, ifp))
			break;

	//-> get seqlen
	rewind(ifp);
	if (!fgets(seq, MAXSEQLEN, ifp))
		fail("Bad alignment file!");
	seqlen = strlen(seq)-1;

	//-> init aln
	aln = new char*[nseqs];
	for(i=0;i<nseqs;i++)
		aln[i]=new char[seqlen+1];

	//-> first sequence
	for (j=0; j<seqlen; j++)
		aln[0][j] = aanum(seq[j]);

	//-> remaining sequences
	for (i=1; i<nseqs; i++)
	{
		if (!fgets(seq, MAXSEQLEN, ifp))
			break;

		if (seqlen != strlen(seq)-1)
			fail("Length mismatch in alignment file!");

		for (j=0; j<seqlen; j++)
			aln[i][j] = aanum(seq[j]);
	}
	fclose(ifp);

	//-> print nseqs and seqlen
	printf("neseq = %d , seqlen = %d \n", nseqs, seqlen);
	/* ----- initialize data structure ------*/ /*end*/



	/* ----- PART II: Calculate sequence weights ----- */ /*start*/
	vector<float> weight(nseqs,0);

	#pragma omp parallel
	{
		vector<int> sim_local(nseqs, 0);
		#pragma omp for schedule(dynamic)
		for (int ii=0; ii<nseqs; ii++)
			for (int jj=ii+1; jj<nseqs; jj++)
			{
				int nthresh = compute_similarity(ii,jj,aln, idthresh, seqlen);
				if (nthresh > 0)
				{
					sim_local[ii] ++;
					sim_local[jj] ++;
				}
			}
		#pragma omp critical
		{
			for(int ii=0; ii<nseqs; ++ii)weight[ii]+=sim_local[ii];
		}
	}

	for (i=0; i<nseqs; i++)
		weight[i] += 1.0;

	double wtsum=0;
	for (i=0; i<nseqs; i++)
		wtsum += (weight[i] = 1.0 / weight[i]);
	/* ----- PART II: Calculate sequence weights ----- */ /*end*/



	/* ----- PART III: Calculate singlet frequencies with pseudocount = 1 */ /*start*/
	vector <vector <double> > pa;
	vector <vector <double> > aaa;
	pa.resize(seqlen);
	aaa.resize(seqlen);
	for(i=0;i<seqlen;i++)pa[i].resize(21,0);
	for(i=0;i<seqlen;i++)aaa[i].resize(21,0);

	#pragma omp parallel
	{
		vector <vector<double> > pa_local;
		vector <vector<double> > aa_local;
		pa_local.resize(seqlen);
		aa_local.resize(seqlen);
		for (int ii=0; ii<seqlen; ii++)pa_local[ii].resize(21,0);
		for (int ii=0; ii<seqlen; ii++)aa_local[ii].resize(21,0);
		#pragma omp for schedule(dynamic)
		for (int ii=0; ii<seqlen; ii++)
			for (int kk=0; kk<nseqs; kk++)
			{
				int aa = aln[kk][ii];
				if (aa < 21 && aa>=0)
				{
					pa_local[ii][aa] += weight[kk];
					aa_local[ii][aa] += weight[kk];
				}
			}
		#pragma omp critical
		{
			for (int ii=0; ii<seqlen; ii++)
				for (int aa=0; aa<21; aa++)
				{
					pa[ii][aa]+=pa_local[ii][aa];
					aaa[ii][aa]+=aa_local[ii][aa];
				}
		}
	}

	for (int ii=0; ii<seqlen; ii++)
		for (int aa=0; aa<21; aa++)
		{
			pa[ii][aa] += 1.0;
			pa[ii][aa] /= 21.0 + wtsum;
		}

	for (int aa=0; aa<21; aa++) 
	{
		aacomp[aa]=0;
		for (int ii=0; ii<seqlen; ii++)
			aacomp[aa]+=aaa[ii][aa];
	}
	/* ----- PART III: Calculate singlet frequencies with pseudocount = 1 */ /*end*/



	/* ----- PART IV: Calculate pariwise frequencies with pseudocount = 1 */ /*start*/
	double **psmat=new double*[seqlen];
	for (i=0; i<seqlen; i++) psmat[i]=new double[seqlen];
	double **mimat=new double*[seqlen];
	for (i=0; i<seqlen; i++) mimat[i]=new double[seqlen];

	#pragma omp parallel for schedule(dynamic)
	for (int ii=0; ii<seqlen; ii++)
		for (int jj=ii+1; jj<seqlen; jj++)
		{
			//--| calculate Contact_Potential
			double potsum=calculate_Potsum(ii,jj,aln, nseqs, wtsum, weight);
			psmat[ii][jj] = potsum;
			psmat[jj][ii] = potsum;
			//--| calculate Mutual_Information
			double mi=calculate_localMI(ii,jj,aln, nseqs, wtsum, pa, weight);
			mimat[ii][jj] = mi;
			mimat[jj][ii] = mi;
		}

	float *misum=new float[seqlen];
	for (int ii=0; ii<seqlen; ++ii)
	{
		misum[ii]=0.0;
		mimat[ii][ii]=0.0;
		psmat[ii][ii]=0.0;
		for (int jj=0; jj<seqlen; jj++)
			misum[ii]+=mimat[ii][jj];
	}

	float mimean=0.0;
	for (int ii=0; ii<seqlen; ii++)
		for (int jj=ii+1; jj<seqlen; jj++)
			mimean += mimat[ii][jj];
	mimean /= seqlen * (seqlen-1) * 0.5;
	/* ----- PART IV: Calculate pariwise frequencies with pseudocount = 1 */ /*end*/



	/* ----- PART V: final output -------- */ /*start*/
	if (!(singofp = fopen(argv[2], "w")))
		fail("Cannot open sing output file!");

	if (!(pairofp = fopen(argv[3], "w")))
		fail("Cannot open pair output file!");

	fprintf(singofp, "%d\n", seqlen);
	fprintf(singofp, "%d\n", nseqs);
	fprintf(singofp, "%f\n", wtsum);

	//-> output singleton header
	for (a=0; a<21; a++)
		fprintf(singofp, "%6.4f%c", aacomp[a] / wtsum / seqlen, (a == 20) ? '\n' : ' ');

	//-> output pairwise
	for (i=0; i<seqlen; i++)
	{
		//--| calculate singleton content
		float entropy=0;
		for (j=0; j<21; j++)
		{
			fprintf(singofp, "%5.3f ", pa[i][j]);
			entropy += pa[i][j] * log(pa[i][j]);
		}
		entropy = -entropy / log(21.0);
		fprintf(singofp, " %6.4f\n", entropy);

		//--| calculate pairwise content
		for (j=i+1; j<seqlen; j++)
		{
			double mip = (float)(mimat[i][j]) - misum[i]*misum[j]/SQR(seqlen-1.0)/mimean;
			fprintf(pairofp, "%d %d %f %f %f\n", i+1, j+1, psmat[i][j], (float)(mimat[i][j]), mip);
		}
	}

	fclose(singofp);
	fclose(pairofp);
	/* ----- PART V: final output -------- */ /*end*/


	/* ------ delete data structure ----- */
	for(i=0;i<nseqs;i++)
		delete [] aln[i];
	delete [] aln;

	for(i=0;i<seqlen;i++)
		delete [] psmat[i];
	delete [] psmat;

	for(i=0;i<seqlen;i++)
		delete [] mimat[i];
	delete [] mimat;

	delete [] misum;

	//---- exit ----//
	return 0;
}
