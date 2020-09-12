/* MetaPSICOV Alignment Stats Program - by David T. Jones June 2014 */

/* Copyright (C) 2014 University College London */

/* V1.01 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))

#define MAXSEQLEN 5000

char  **aln;
int nseqs;

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVXXX";

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
    void          **p, *rp;

    rp = malloc(rows * sizeof(void *) + sizeof(int));

    if (rp == NULL)
	fail("allocmat: malloc [] failed!");

    *((int *)rp) = rows;

    p = rp + sizeof(int);

    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

int             main(int argc, char **argv)
{
    int             a, b, c, d, i, j, ii, jj, k, l, n, llen, maxnps=3, posn, qstart, seqlen, endflg = FALSE, tots, maxtots, lstart, nids, s;
    double aacomp[21], sumj[21], sum, score, pab[21][21], pa[MAXSEQLEN][21], pav[MAXSEQLEN][21], di, pdir, hahb, hx, hy, hxy, mi, mip, z, oldvec[21], change, wtsum, idthresh=0.38, potsum;
    float *weight;
    char            buf[4096], name[512], seq[MAXSEQLEN], qseq[160], *cp;
    FILE *ifp, *singofp, *pairofp;

    if (argc != 4)
	fail("Usage: alnstats alnfile singoutfile pairoutfile");

    ifp = fopen(argv[1], "r");
    if (!ifp)
	fail("Unable to open alignment file!");

    for (nseqs=0;; nseqs++)
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;

    aln = malloc(nseqs * sizeof(char *));
    weight = malloc(nseqs * sizeof(float));

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

	
   printf("neseq = %d \n", nseqs );    
/* Calculate sequence weights */
    for (i=0; i<nseqs; i++)
	weight[i] = 1.0;
    
    for (i=0; i<nseqs; i++)
    {
	for (j=i+1; j<nseqs; j++)
	{
	    int nthresh = idthresh * seqlen;
	
	    for (k=0; nthresh > 0 && k<seqlen; k++)
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
	    aacomp[a] += weight[k];
	}
	
	for (a=0; a<21; a++)
	    pa[i][a] /= 21.0 + wtsum;
    }

    float misum[MAXSEQLEN], mimean=0.0, **mimat;

    mimat = allocmat(seqlen, seqlen, sizeof(float));
    
    for (i=0; i<seqlen; i++)
	misum[i] = 0.0;
	    
    for (i=0; i<seqlen; i++)
    {
	for (j=i+1; j<seqlen; j++)
	{
	    /* Calculate pair frequencies with pseudocount = 1 */
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		    pab[a][b] = 1.0 / 21.0;
		    
	    for (k=0; k<nseqs; k++)
	    {
		a = aln[k][i];
		b = aln[k][j];
		if (a < 21 && b < 21)
		    pab[a][b] += weight[k];
	    }
	    
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		{
		    pab[a][b] /= 21.0 + wtsum;
		    pab[b][a] = pab[a][b];
		    
//		    printf("%d/%d %d/%d %f %f %f %f\n", i+1, a, j+1, b, pab[i][j][a][b], pa[i][a] , pa[j][b], pab[i][j][a][b] - pa[i][a] * pa[j][b]);
		}

	    for (mi=a=0; a<20; a++)
		for (b=0; b<20; b++)
		    mi += pab[a][b] * log(pab[a][b]/pa[i][a]/pa[j][b]);

	    mimat[i][j] = mimat[j][i] = mi;
	    misum[i] += mi;
	    misum[j] += mi;
	    mimean += mi;
	}
    }
    
    mimean /= seqlen * (seqlen-1) * 0.5;

    if (!(singofp = fopen(argv[2], "w")))
	fail("Cannot open sing output file!");

    if (!(pairofp = fopen(argv[3], "w")))
	fail("Cannot open pair output file!");

    fprintf(singofp, "%d\n", seqlen);
    fprintf(singofp, "%d\n", nseqs);
    fprintf(singofp, "%f\n", wtsum);

    for (a=0; a<21; a++)
	fprintf(singofp, "%6.4f%c", aacomp[a] / wtsum / seqlen, (a == 20) ? '\n' : ' ');

    for (i=0; i<seqlen; i++)
    {
	float entropy;
	
	for (entropy=j=0; j<21; j++)
	{
	    fprintf(singofp, "%5.3f ", pa[i][j]);
	    entropy += pa[i][j] * log(pa[i][j]);
	}

	entropy = -entropy / log(21.0);

	fprintf(singofp, " %6.4f\n", entropy);
	
	for (j=i+1; j<seqlen; j++)
	{
	    for (potsum=k=0; k<nseqs; k++)
	    {
		a = aln[k][i];
		b = aln[k][j];
		if (a < 20 && b < 20)
		    potsum += weight[k] * contmat[a][b];
	    }

	    potsum /= wtsum;

	    mip = mimat[i][j] - misum[i]*misum[j]/SQR(seqlen-1.0)/mimean;
	    
	    fprintf(pairofp, "%d %d %f %f %f\n", i+1, j+1, potsum, mimat[i][j], mip);
	}
    }
    
    fclose(singofp);
    fclose(pairofp);
    
    return 0;
}
