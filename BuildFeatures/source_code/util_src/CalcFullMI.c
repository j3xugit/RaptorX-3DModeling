/* DeepCOV Alignment Stats Program - by David T. Jones June 2017 */

/* Copyright (C) 2017 University College London */

/* V1.00 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))

#define MAXSEQLEN 2000


/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Convert AA letter to numeric code (0-21) */
int
                aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
	21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p;

    p = malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

int             main(int argc, char **argv)
{
    int             a, b, i, j, k, l, nseqs, llen, posn, seqlen, ar=1, pairflg=0;
    double sum, score, di, pdir, wtsum, idthresh=0.38, (** pab)[21][21];
    float pa[MAXSEQLEN][21];
    float logpa[MAXSEQLEN][21];
    float *weight, val;
    char  buf[4096], name[512], seq[MAXSEQLEN], **aln;
    int flag;

    FILE *ifp, *pairofp;

    if (argc < 3 )
	fail("Usage: CalcFullMI alnfile outputfile [mode]");

    ifp = fopen(argv[1], "r");
    if (!ifp)
	fail("Unable to open alignment file!");

    /* when flag=0, calculate full MI matrix. otherwise, full covariance matrix */
    flag = 0; 
    if (argc>=4)
	flag = atoi(argv[3]);
    if (flag != 1)
	flag = 0;

    for (nseqs=0;; nseqs++)
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;

    if (!(aln = malloc(nseqs * sizeof(char *))) || !(weight = malloc(nseqs * sizeof(float))))
	fail("Out of memory!");

    rewind(ifp);
    
    if (!fgets(seq, MAXSEQLEN, ifp))
	fail("Bad alignment file!");
    
    seqlen = strlen(seq)-1;

    if (!(aln[0] = malloc(seqlen)))
	fail("Out of memory!");

    for (j=0; j<seqlen; j++)
	aln[0][j] = aanum(seq[j]);
    
    for (i=1; i<nseqs; i++)
    {
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;
	
	if (seqlen != strlen(seq)-1)
	    fail("Length mismatch in alignment file!");
	
	if (!(aln[i] = malloc(seqlen)))
	    fail("Out of memory!");
	
	for (j=0; j<seqlen; j++)
	    aln[i][j] = aanum(seq[j]);
    }

    fclose(ifp);

    pab = allocmat(seqlen, seqlen, sizeof(double) * 21 * 21);

    /* Calculate sequence weights */
    for (i=0; i<nseqs; i++)
	weight[i] = 1.0;
    
    for (i=0; i<nseqs; i++)
    {
	for (j=i+1; j<nseqs; j++)
	{
	    int nthresh = idthresh * seqlen;
	
	    for (k=0; nthresh >= 0 && k<seqlen; k++)
		if (aln[i][k] != aln[j][k])
		    nthresh--;
	    
	    if (nthresh > 0)
	    {
		weight[i]++;
		weight[j]++;
	    }
	}
    }

    for (wtsum=i=0; i<nseqs; i++)
	wtsum += (weight[i] = 1.0 / weight[i]);

    /* Calculate singlet frequencies with pseudocount = 1 */
    for (i=0; i<seqlen; i++)
    {
	for (a=0; a<21; a++)
	    pa[i][a] = 1.0;
	
	for (k=0; k<nseqs; k++)
	{
	    a = aln[k][i];
	    if (a < 21)
		pa[i][a] += weight[k];
	}
	
	for (a=0; a<21; a++)
	    pa[i][a] /= 21.0 + wtsum;
    }

    for (i=0; i<seqlen; i++)
	for (a=0; a<21; a++)
		logpa[i][a] = log(pa[i][a]);

    for (i=0; i<seqlen; i++)
    {
	for (j=0; j<seqlen; j++)
	{
	    /* Calculate pair frequencies with pseudocount = 1 */
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		    pab[i][j][a][b] = 1.0 / 21.0;
		    
	    for (k=0; k<nseqs; k++)
	    {
		a = aln[k][i];
		b = aln[k][j];
		if (a < 21 && b < 21)
		    pab[i][j][a][b] += weight[k];
	    }
	    
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		{
		    pab[i][j][a][b] /= 21.0 + wtsum;
		}
	}
    }
    
    if (!(pairofp = fopen(argv[2], "wb")))
	fail("Cannot open output file!");

    /* output seq length */
    fwrite(&seqlen, sizeof(int), 1, pairofp);

    /* output position-specific frequency matrix, which has 21 times seqlen values */
    if (fwrite(pa, sizeof(float), 21*seqlen, pairofp) != 21*seqlen)
	fail("Write failed!");

    for (i=0; i<seqlen; i++)
	for (j=0; j<seqlen; j++)
    		for (a=0; a<21; a++)
			for (b=0; b<21; b++)
			{
				if (flag==1)
		    			val = pab[i][j][a][b] - pa[i][a] * pa[j][b];
				else
		    			val = log(pab[i][j][a][b]) - logpa[i][a] - logpa[j][b];
		    		if (fwrite(&val, sizeof(float), 1, pairofp) != 1)
					fail("Write failed!");
			}
    
    fclose(pairofp);
    
    return 0;
}
