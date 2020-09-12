/* RaptorX Alignment Stats Program - by Jinbo Xu */
/* adapted from DeepCov */

/* V1.00 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAXSEQLEN 2000


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
		999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2, 21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    	};

    	return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Allocate matrix */
void   *allocmat(int rows, int columns, int size)
{
	int i;
    	void  **p;

    	p = malloc(rows * sizeof(void *));

    	if (p == NULL)
		fail("allocmat: malloc [] failed!");
    	for (i = 0; i < rows; i++)
		if ((p[i] = calloc(columns, size)) == NULL)
	    		fail("allocmat: malloc [][] failed!");

    	return p;
}

double CalcSeqWeight(int** aln, int nseqs, int seqlen, float* weight){
	int i, j, k;
	float idthresh = 0.38;
	float wtsum;
	
    	for (i=0; i<nseqs; i++)
		weight[i] = 1.0;

    	for (i=0; i<nseqs; i++){
		for (j=i+1; j<nseqs; j++){
	    		int nthresh = idthresh * seqlen;
	
	    		for (k=0; nthresh >= 0 && k<seqlen; k++)
				if (aln[i][k] != aln[j][k])
		    			nthresh--;
	    
	    		if (nthresh > 0){
				weight[i]++;
				weight[j]++;
	    		}
		}
    	}

	wtsum = 0;
	for (i=0; i<nseqs; i++){
		float tmp = 1.0/weight[i];
		weight[i] = tmp;
		wtsum += tmp;
	}
	return wtsum;
}

/* aln, weight, nseqs, seqlen, bounds are input; fullmatrix and simplematrix are results */
/* aln is the MSA represented by numbers (gap shall be represented by 20), weight is the sequence weight */
/* if weight is set to 0, then the function will calculate the sequence weight */
void CalcPairMatrixFromMSA_internal(int** aln, float* weight, int nseqs, int seqlen, int bounds[4], int flag, float* fullmatrix, float* simplematrix)
{
	int a, b, i, j, k;
	int i2, j2, start;
    	float *** pab, pa[MAXSEQLEN][21];
	int blockSize;
	int top, left, bottom, right;
	int nRows, nCols;
	float wtsum;

	wtsum = 0.0;
	for (i=0; i<nseqs; i++)
		wtsum += weight[i];
	/* in principle, wtsum shall be at least 1; if weight is not set, we calculate it here */
	if (wtsum < 0.5){
		wtsum = CalcSeqWeight(aln, nseqs, seqlen, weight);
	}

	/* set the boundary. Sometimes we just need to calculate a submatrix to save time */
	top = bounds[0];
	left = bounds[1];
	bottom = bounds[2];
	right = bounds[3];

	if (top<0)
		top =0;
	if (left<0)
		left =0;
	if (bottom<0)
		bottom = seqlen;
	if (right<0)
		right = seqlen;
	
	nRows = bottom - top;
	nCols = right - left;
	blockSize = 21 * 21;

	/*allocate memory for pab and map pab to fullmatrix */
	pab = malloc( nRows * sizeof(void*) );
	for (i=0; i<nRows; i++){
		float* tmp = fullmatrix + i * nRows * blockSize;
		pab[i] = malloc( nCols * sizeof(float*) );
		for (j=0; j<nCols; j++)
			pab[i][j] = tmp + j*blockSize;
	}

    	/* Calculate singlet frequencies with pseudocount = 1 */
    	for (i=0; i<seqlen; i++){
		for (a=0; a<21; a++)
	    		pa[i][a] = 1.0;
	
		for (k=0; k<nseqs; k++){
	    		a = aln[k][i];
	    		if (a < 21)
				pa[i][a] += weight[k];
		}
	
		for (a=0; a<21; a++)
	    		pa[i][a] /= 21.0 + wtsum;
    	}

	/* Calculate pair frequencies with pseudocount = 1 */
    	for (i=top; i<bottom; i++){
		i2 = i-top;
		for (j=left; j<right; j++){
			j2 = j-left;
	    		for (a=0; a<21; a++){
				start = a*21;
				for (b=0; b<21; b++)
		    			pab[i2][j2][start+b] = 1.0 / 21.0;
			}
		    
	    		for (k=0; k<nseqs; k++){
				a = aln[k][i];
				b = aln[k][j];
				if (a < 21 && b < 21)
		    			pab[i2][j2][a*21+b] += weight[k];
	    		}
	    
	    		for (a=0; a<21; a++){
				start = a*21;
				for (b=0; b<21; b++){
		    			pab[i2][j2][start+b] /= 21.0 + wtsum;
				}
			}
		}
    	}

	/* calculate the final results */
	if (flag==1){
   		/* calc full covariance matrix */ 
		float val;
		int start;

		for (i=top; i<bottom; i++){
			i2 = i-top;
			for (j=left; j<right; j++){
				j2 = j-left;
    				for (a=0; a<21; a++){
					start = a*21;
					for (b=0; b<21; b++){
		    				val = pab[i2][j2][start+b] - pa[i][a] * pa[j][b];
						pab[i2][j2][start+b] = val;
					}
				}
			}
		}

	}else{
		/* calc full MI matrix */
		float val;
		int start;

		/* calculate log of pa */
		float lgpa[MAXSEQLEN][21];
		for (i=0; i<seqlen; i++)
			for (a=0; a<21; a++)
				lgpa[i][a] = log(pa[i][a]);

		for (i=top; i<bottom; i++){
			i2 = i-top;
			for (j=left; j<right; j++){
				j2 = j-left;
    				for (a=0; a<21; a++){
					start = a*21;
					for (b=0; b<21; b++){
		    				val = log(pab[i2][j2][start+b]) - lgpa[i][a] -lgpa[j][b];
						pab[i2][j2][start+b] = val;
					}
				}
			}
		}
	}

	/* add code here to generate results for simplematrix */
    

	/* free pab */
	for (i=0; i<nRows; i++)
		free(pab[i]);
	free(pab);
    
}

/* convert an MSA in string representation to numeric representation*/
/* the caller shall allocate and free memory for result*/
void Seq2Numbers(char** aln, int nseqs, int seqlen, int** result){
	int i, j;
	for (i=0; i<nseqs; i++)
		for(j=0; j<seqlen; j++)
			result[i][j] = aanum(aln[i][j]);
}

void CalcPairMatrixFromMSA(char** aln, float* weight, int nseqs, int seqlen, int bounds[4], int flag, float* fullmatrix, float* simplematrix)
{

	int** MSA;
	int i;

	/* allocate memory for MSA */
	MSA = malloc(nseqs * sizeof(int*));
	for (i=0; i<nseqs; i++)
		MSA[i] = malloc(sizeof(int) * seqlen);

	Seq2Numbers(aln, nseqs, seqlen, MSA);
	CalcPairMatrixFromMSA_internal(MSA, weight, nseqs, seqlen, bounds, flag, fullmatrix, simplematrix);

	/* free memory */
	for (i=0; i<nseqs; i++)
		free(MSA[i]);
	free(MSA);
}


/* read an MSA and return its seqlen */
/* ifp is the input, aln stores the resultant MSA in numeric form */
/* the caller shall allocate memory for aln and weight and then free them */
/* currently weight is not used, but in future we may read seq weight from file */
int ReadMSA(char* alnfile, int nseqs, int seqlen, int** aln, float* weight)
{
	char seq[MAXSEQLEN];
	int i, j;
	FILE* ifp;

	ifp = fopen(alnfile, "r");
        if (!ifp){
                fail("Unable to open alignment file!");
                fail(alnfile);
        }

	/* read protein sequences and convert them into numeric representation */
    	for (i=0; i<nseqs; i++){
		if (!fgets(seq, MAXSEQLEN, ifp)){
			fail("the file does not have enough sequences");
			fail(alnfile);
		}
	
		if (seqlen != strlen(seq)-1)
	    		fail("Length mismatch in alignment file!");
	
    		for (j=0; j<seqlen; j++)
			aln[i][j] = aanum(seq[j]);
    	}

	fclose(ifp);

	/* for debug 
	for (i=0; i<nseqs; i++){
		printf("\n");
		for (j=0; j<seqlen; j++)
			printf("%d ", aln[i][j]);
	}
	*/
}

/* determine the number of sequences and seqlen in an MSA */
/* the results are stored in nseqs and seqlen */
int DetermineNumNLen(char* alnfile, int* nseqs_ptr, int* seqlen_ptr){
	char seq[MAXSEQLEN];
	int nseqs, seqlen;

	FILE* ifp;
    	ifp = fopen(alnfile, "r");
    	if (!ifp){
		fail("Unable to open alignment file!");
		fail(alnfile);
	}
	
	fgets(seq, MAXSEQLEN, ifp);
	seqlen = strlen(seq)-1;

	/* count the number of sequences */
    	for (nseqs=1;; nseqs++)
		if (!fgets(seq, MAXSEQLEN, ifp))
	    		break;

	*nseqs_ptr = nseqs;
	*seqlen_ptr = seqlen;
	fclose(ifp);
}

/* alnfile, seqLen, bounds=[top, left, bottom, right], flag are input; fullmatrix and simplematrix are results */
/* flag is used to indicate if covariance matrix or mutual information matrix shall be calculated */
/* flag =1 for covriance matrix and 0 for MI matrix */
void CalcPairMatrixFromFile(char* alnfile, int bounds[4], int flag, float* fullmatrix, float* simplematrix)
{
	int i;
    	float *weight;
    	char  seq[MAXSEQLEN];
	int **aln;
	FILE* ifp;
	int nseqs, seqlen;

	DetermineNumNLen(alnfile, &nseqs, &seqlen);

	/* allocate memory */
	aln = malloc( sizeof(int*) * nseqs);
	for (i=0; i<nseqs; i++)
		aln[i] = malloc( sizeof(int) * seqlen);

    	weight = malloc(nseqs * sizeof(float));
    	if (!weight)
		perror("Cannot allocate memory to weight!");
	memset(weight, 0, nseqs * sizeof(float) );
	/*
	for (i=0; i<nseqs; i++)
		weight[i] = 1.0;
	*/

	ReadMSA(alnfile, nseqs, seqlen, aln, weight);

	CalcPairMatrixFromMSA_internal(aln, weight, nseqs, seqlen, bounds, flag, fullmatrix, simplematrix);

	/* free memory */ 
	free(weight);

	for (i=0; i<nseqs; i++)
		free(aln[i]);
	free(aln);
	
}

int main(int argc, char** argv){
	char *alnfile, *outfile;
	int nseqs, seqLen, flag, bounds[4];
	float *fullmatrix, *simplematrix;
	long fmSize;
	int i, j, k;

	if (argc<2)
		fail("Usage: CLib4MSA alnfile [mode]");
	alnfile = argv[1];
	flag = 1;
	if (argc>=3){
		flag = atoi(argv[2]);
		if (flag != 1)
			flag=0;
	}
		
	DetermineNumNLen(alnfile, &nseqs, &seqLen);

	printf("#seqs=%d, seqLen=%d\n", nseqs, seqLen);
	
	bounds[0] = 0;
	bounds[1] = 0;
	bounds[2] = -1;
	bounds[3] = -1;
	
	fmSize = seqLen * seqLen * 21 * 21;

	fullmatrix = malloc(fmSize * sizeof(float) );
	if (!fullmatrix)
		perror("Cannot allocate memory to fullmatrix!");

	simplematrix = malloc( seqLen*seqLen * sizeof(float) );
	if (!simplematrix)
		perror("Cannot allocate memory to simplematrix!");

	CalcPairMatrixFromFile(alnfile, bounds, flag, fullmatrix, simplematrix);

	/* print some results */
	for (i=0; i<3; i++)
		for(j=i+1; j<5; j++){
			printf("%d %d: ", i, j);
			for (k=0; k<441; k++)
				printf("%.6f ", fullmatrix[ i*seqLen*441 + j*441 + k]);
			printf("\n");
		} 

	/* free memory */
	free(fullmatrix);
	free(simplematrix);

	return 0;	
}
