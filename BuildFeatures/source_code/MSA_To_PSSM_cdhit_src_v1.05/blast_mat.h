#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

//----- BLAST AA <-> ARND mapping ---------//
extern int BLAST_AA_To_ARND_Mapping(char amino);
extern char BLAST_ARND_To_AA_Mapping(int num);

//----- BLAST ARND <-> ACDE mapping -------//
extern int BLAST_ACDE_To_ARND_Mapping(int num);
extern int BLAST_ARND_To_ACDE_Mapping(int num);

//----- BLAST 28-digit <-> ARND mapping ---//
extern int BLAST_Digit_To_ARND_Mapping(int amino);
extern int BLAST_ARND_To_Digit_Mapping(int amino);

/** Retrieve the background letter probabilities implicitly used in
 * constructing the score matrix matrix_name. */
extern int Blast_FrequencyDataIsAvailable(const char *matrix_name);

/** Return true if substitution data is available for the given matrix name. */
extern int Blast_GetSubstituteForMatrix(int *sub_sco, const char *matrix_name);

/** Get joint probabilities for the named matrix. */
extern int Blast_GetJointProbsForMatrix(double *probs, const char *matrixName);

/** Return true if frequency data is available for the given matrix name. */
extern int Blast_GetMatrixBackgroundFreq(double *back, const char *matrix_name);

/** Return true if frequency ratio is available for the given matrix name. */
extern int Blast_GetFreqRatioForMatrix(double *freq_ratio, const char *matrix_name);

//----- get matrix-related gap values -----//
extern int Blast_GetMatrixRelatedGapValues(double &lambda,double &H,double &K,
	int gap_open,int gap_extend,const char *matrix_name);

/** Get standard background frequency */
extern void Blast_GetStandardBackgroundFreq(double *std_back);
