#include "blast_pssm.h"
#include <math.h>


/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
BLAST_PSSM::BLAST_PSSM(void)
{
	//-> assign macros
	BLOCK_LEN = 5;                 //for block length less than BLOCK_LEN, apply blastpgp
	// normal macros
	MAX_IND_OBSERVATIONS=400;
	PSEUDO_MAX=1000000;
	kZeroObsPseudo = 30.0;         //arbitrary constant to use for columns with zero observations
	kPSIScaleFactor=200;           //used in Scaling_PSSM_Matrix, scaling factor
	kPositScalingNumIterations=10; //used in Scaling_PSSM_Matrix, iteration numbers
	kPositScalingPercent=0.05;     //used in Scaling_PSSM_Matrix, percent variance
	// constant values, were #defines in posit.c // -> used in ColumnSpecific_Pseudocounts
	kPseudoMult = 500.0;       /*Was PSEUDO_MULTIPLIER */
	kPseudoNumerator = 0.0457; /*numerator of entropy-based method, was PSEUDO_NUMERATOR */
	kPseudoExponent = 0.8;     /*exponent of denominator, was PSEUDO_EXPONENT */
	kPseudoSmallInitial = 5.5; /*small number of pseudocounts to avoid 0 probabilities */
	kPosEpsilon = 0.0001;      /*minimum return value of s_computeRelativeEntropy */

	//-> init data structure
	// query length
	Query_Length=0;
	Cur_Lambda=-1;
	Cur_H=-1;
	Cur_K=-1;
	// aligned block
	Aligned_Block_Left=NULL;
	Aligned_Block_Right=NULL;
	Temp_Block_Left=NULL;
	Temp_Block_Right=NULL;
	// pseudo count
	PosDistinctDistrib=NULL;
	PosDistinctNum=NULL;
	PosNumParticipating=NULL;
	Pseudo_Prob=NULL;
	// sequence weight
	PosObservations=NULL;
	Match_Weights=NULL;
	Colum_Weights=NULL;
	Colum_Weights_P=NULL;
	// PSSM
	Freq_Ratio=NULL;
	Freq_Entro=NULL;
	PSSM_Ini=NULL;
	PSSM_Fin=NULL;

	//-> temporary data structure
	Pseudo_Expno=new double[MAX_IND_OBSERVATIONS];
	Inter_Prob=new double[BLAST_AASIZE];
	Column_Residue_Num=new int[BLAST_AASIZE];
}

/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
BLAST_PSSM::~BLAST_PSSM(void)
{
	//aligned block
	if(Aligned_Block_Left!=NULL)delete [] Aligned_Block_Left;
	if(Aligned_Block_Right!=NULL)delete [] Aligned_Block_Right;
	if(Temp_Block_Left!=NULL)delete [] Temp_Block_Left;
	if(Temp_Block_Right!=NULL)delete [] Temp_Block_Right;
	//pseudo count
	if(PosDistinctDistrib!=NULL)delete [] PosDistinctDistrib;
	if(PosDistinctNum!=NULL)delete [] PosDistinctNum;
	if(PosNumParticipating!=NULL)delete [] PosNumParticipating;
	if(Pseudo_Prob!=NULL)delete [] Pseudo_Prob;
	//sequence weight
	if(PosObservations!=NULL)delete [] PosObservations;
	if(Match_Weights!=NULL)delete [] Match_Weights;
	if(Colum_Weights!=NULL)delete [] Colum_Weights;
	if(Colum_Weights_P!=NULL)delete [] Colum_Weights_P;
	//PSSM
	if(Freq_Ratio!=NULL)delete [] Freq_Ratio;
	if(Freq_Entro!=NULL)delete [] Freq_Entro;
	if(PSSM_Ini!=NULL)delete [] PSSM_Ini;
	if(PSSM_Fin!=NULL)delete [] PSSM_Fin;
	//temporary
	delete [] Pseudo_Expno;
	delete [] Inter_Prob;
	delete [] Column_Residue_Num;
}

//================= data create/delete operations ===================//
void BLAST_PSSM::Init_Pseudo_Expno(void)
{
	Init_ExpNum_Observations(Freq_Back,Pseudo_Expno);    
}
void BLAST_PSSM::Create_DataStructure(int query_len)
{
	Query_Length=query_len;
	//aligned block
	Aligned_Block_Left=new int[Query_Length];
	Aligned_Block_Right=new int[Query_Length];
	Temp_Block_Left=new int[Query_Length];
	Temp_Block_Right=new int[Query_Length];
	//pseudo count
	PosDistinctDistrib=new int[Query_Length*(BLAST_AASIZE+1)];
	PosDistinctNum=new int[Query_Length];
	PosNumParticipating=new int[Query_Length];
	Pseudo_Prob=new double[Query_Length];
	//sequence weight
	PosObservations=new double[Query_Length];
	Match_Weights=new double[Query_Length*BLAST_AASIZE];
	Colum_Weights=new double[Query_Length];
	Colum_Weights_P=new double[Query_Length];
	//PSSM
	Freq_Ratio=new double[Query_Length*BLAST_AASIZE];
	Freq_Entro=new double[Query_Length];
	PSSM_Ini=new int[Query_Length*BLAST_AASIZE];
	PSSM_Fin=new int[Query_Length*BLAST_AASIZE];
}
void BLAST_PSSM::Delete_DataStructure(void)
{
	//aligned block
	if(Aligned_Block_Left!=NULL)
	{
		delete [] Aligned_Block_Left;
		Aligned_Block_Left=NULL;
	}
	if(Aligned_Block_Right!=NULL)
	{
		delete [] Aligned_Block_Right;
		Aligned_Block_Right=NULL;
	}
	if(Temp_Block_Left!=NULL)
	{
		delete [] Temp_Block_Left;
		Temp_Block_Left=NULL;
	}
	if(Temp_Block_Right!=NULL)
	{
		delete [] Temp_Block_Right;
		Temp_Block_Right=NULL;
	}
	//pseudo count
	if(PosDistinctDistrib!=NULL)
	{
		delete [] PosDistinctDistrib;
		PosDistinctDistrib=NULL;
	}
	if(PosDistinctNum!=NULL)
	{
		delete [] PosDistinctNum;
		PosDistinctNum=NULL;
	}
	if(PosNumParticipating!=NULL)
	{
		delete [] PosNumParticipating;
		PosNumParticipating=NULL;
	}
	if(Pseudo_Prob!=NULL)
	{
		delete [] Pseudo_Prob;
		Pseudo_Prob=NULL;
	}
	//sequence weight
	if(PosObservations!=NULL)
	{
		delete [] PosObservations;
		PosObservations=NULL;
	}
	if(Match_Weights!=NULL)
	{
		delete [] Match_Weights;
		Match_Weights=NULL;
	}
	if(Colum_Weights!=NULL)
	{
		delete [] Colum_Weights;
		Colum_Weights=NULL;
	}
	if(Colum_Weights_P!=NULL)
	{
		delete [] Colum_Weights_P;
		Colum_Weights_P=NULL;
	}
	//PSSM
	if(Freq_Ratio!=NULL)
	{
		delete [] Freq_Ratio;
		Freq_Ratio=NULL;
	}
	if(Freq_Entro!=NULL)
	{
		delete [] Freq_Entro;
		Freq_Entro=NULL;
	}
	if(PSSM_Ini!=NULL)
	{
		delete [] PSSM_Ini;
		PSSM_Ini=NULL;
	}
	if(PSSM_Fin!=NULL)
	{
		delete [] PSSM_Fin;
		PSSM_Fin=NULL;
	}
}


//=================  pseudocount related operations ====================//

/** initialize the expected number of observations
  use background probabilities for this matrix
  Calculate exp. # of distinct aa's as a function of independent trials
  copy of posit.c:initializeExpNumObservations

  @param expno table of expectations [out]
  @param backgroundProbabilities residue background probs [in]
*/
//[note]: since *expno only depends on back_prob, we could calculate them at head
void BLAST_PSSM::Init_ExpNum_Observations(
	double *back_prob,                                       //-> input, background probability
	double *expno)                                           //-> output, expno table of expectations
{
	int     j,k ; /*loop indices*/
	double  weighted_sum; /*20 - this is how many distinct amino acids are expected*/
	expno[0] = 0;
	for (j=1;j<MAX_IND_OBSERVATIONS;++j)
	{
		weighted_sum = 0;
		for (k=0;k<BLAST_AASIZE;++k)weighted_sum += exp(j*log(1.0-back_prob[k]));
		expno[j] = BLAST_AASIZE-weighted_sum;
	}
}


/** A method to estimate the effetive number of observations
  in the interval for the specified columnNumber 
  copy of posit.c:effectiveObservations
  @param align_blk data structure describing the aligned blocks [in]
  @param seq_weights data structure of sequence weights [in]
  @param columnNumber column in the PSSM [in]
  @param queryLength length of the query sequence
  @param expno table of expectations [in]
*/
//[note]: this function is block specific, so we should only calculate once for each block
double BLAST_PSSM::Effective_Observations(                 //-> return, effective observations
	int *align_blk_left,int *align_blk_right,                //-> input, align block left/right boundary
	int position, double *expno,                             //-> input, table of expectations [fixed]
	int *posDistrib, int *posNum)                            //-> input, number of sequences at each position
{
	int     i,k; /*loop indices*/
	double  indep;           /*number of independent observations to return*/
	int halfNumColumns;      /*half the number of columns in the interval, rounded down*/
	int totalDistinctCounts; /*total number of distinct letters in columns used*/
	double aveDistinctAA;    /*average number of distinct letters in columns used*/
	int columnsAccountedFor; /*how many of the columns had their
	                            distinct count totaled so far*/

/*  Calculate the average number of distinct amino acids in the half of the
    columns within the block in question with the most distinct amino acids;
    +2 in the parentheses is for rounding up.*/    
	int halfNumColumns_ = (align_blk_right[position]-align_blk_left[position]+2)/2;
	halfNumColumns = halfNumColumns_>1?halfNumColumns_:1;

	//calc AverageDistinctAA
	k = BLAST_AASIZE;
	columnsAccountedFor = 0;
	totalDistinctCounts = 0;
	while (columnsAccountedFor < halfNumColumns) 
	{
		totalDistinctCounts += (posDistrib[position*(BLAST_AASIZE+1)+k] *k);
		columnsAccountedFor += posDistrib[position*(BLAST_AASIZE+1)+k];
		if (columnsAccountedFor > halfNumColumns) 
		{
			totalDistinctCounts -= ((columnsAccountedFor - halfNumColumns) * k);
			columnsAccountedFor = halfNumColumns;
		}
		k--;
	}
	aveDistinctAA = ((double) totalDistinctCounts)/((double) columnsAccountedFor);

/*  Then use the following code to calculate the number of
        independent observations corresponding to
        aveDistinctAA. */
	for (i=1;i<MAX_IND_OBSERVATIONS && expno[i]<=aveDistinctAA;++i);
	indep = (i==MAX_IND_OBSERVATIONS) ? i : (i-(expno[i]-aveDistinctAA)/(expno[i]-expno[i-1]));
	indep = (indep < posNum[position]) ? indep : posNum[position];
	indep = ((indep - 1) > 0) ? (indep - 1) : 0;
	return indep;
}


/** copy of posit.c:columnSpecificPseudocounts 
 @param posSearch data structure of sequence weights [in]
 @param columnNumber column in the PSSM [in]
 @param backgroundProbabilities residue background probs [in]
 @param observations for each column an estimate of observed residues [in]
*/
double BLAST_PSSM::ColumnSpecific_Pseudocounts(            //-> return, eistimated pseudo count
	double *input_prob,                                      //-> input, column specific prob, size is 20
	double *inter_prob,                                      //-> input, just an intermediate term, size is 20
	double *back_prob,                                       //-> input, matrix-specific back prob, size is 20
	double observations)                                     //-> input, an estimate of observed residues
{
	int c;
  /* Constant values, were #defines in posit.c */
//	double kPseudoMult = 500.0;       /*Was PSEUDO_MULTIPLIER */
//	double kPseudoNumerator = 0.0457; /*numerator of entropy-based method, was PSEUDO_NUMERATOR */
//	double kPseudoExponent = 0.8;     /*exponent of denominator, was PSEUDO_EXPONENT */
//	double kPseudoSmallInitial = 5.5; /*small number of pseudocounts to avoid 0 probabilities */
//	double kPosEpsilon = 0.0001;      /*minimum return value of s_computeRelativeEntropy */
 
	//-> 1. adjust the probabilities by assigning observations weight
//	to initialProbabilities and standardWeight to standardProbabilities
	double overallSum = 0.0;
	for(c = 0; c < BLAST_AASIZE; c++) 
	{
		inter_prob[c] = input_prob[c]*observations + back_prob[c]*kPseudoSmallInitial;
		overallSum += inter_prob[c];
	}
	for(c = 0; c < BLAST_AASIZE; c++) inter_prob[c] = inter_prob[c]/overallSum;

	//-> 2. compute relative entropy of first distribution to second distribution
	double relativeEntropy = 0.0;    /*relative entropy of this column to background probs.*/
	for(c = 0; c < BLAST_AASIZE; c++)
	{
		if (inter_prob[c] > kPosEpsilon)
			relativeEntropy += (inter_prob[c] * log (inter_prob[c]/back_prob[c]));
	}
	if (relativeEntropy < kPosEpsilon) relativeEntropy = kPosEpsilon;

	//-> 3. calculate entropy-based pseudo count
	double pseudoDenominator = pow(relativeEntropy, kPseudoExponent); /*intermediate term*/
	double alpha = kPseudoNumerator/pseudoDenominator;                /*intermediate term*/
	double returnValue;                                        /* final returned pseudo count*/
	if (alpha < (1.0 - kPosEpsilon))returnValue = kPseudoMult * alpha/ (1- alpha);
	else returnValue = PSEUDO_MAX;

	//-> 4. final return
	return returnValue;
}

//====================== vice calculation steps ========================//

/** Calculates the position based weights using a modified version of the
 * Henikoff's algorithm presented in "Position-based sequence weights". 
 * Skipped optimization about identical previous sets.
 * @param msa multiple sequence alignment data structure [in]
 * @param aligned_blocks aligned regions' extents [in]
 * @param position position of the query to calculate the sequence weights for
 * [in]
 * @param aligned_seqs array containing the indices of the sequences 
 * participating in the multiple sequence alignment at the requested 
 * position [in]
 * @param seq_weights sequence weights data structure [out]
 */
void BLAST_PSSM::Compute_SeqPos_Weights(
	vector <string> &MSA_in, int *column_residue_num,        //-> input, MSA must be valid and uppere case !!
	int *aligned_blocks_left, int *aligned_blocks_right,     //-> input, aligned blocks
	int position, vector <int> &aligned_seqs,                //-> input, aligned sequences for the given position
	int *posDistrib, int *posDistinct,                       //-> output for positional distribute and distinct
	vector <double> &norm_seq_weights)                       //-> output for positional normalized sequence weight
{
	/* Index into aligned block for requested position */
	int i = 0;
	int k = 0;

	/* Index into array of aligned sequences */
	int asi = 0;
	int sigma = 0;

	/* This flag will be true if more than one different type of residue is
	* found in a column in the extent of the alignment that corresponds to 
	* the position being examined. (replaces Sigma in old code) */
	int distinct_residues_found = 0;

	/* Number of different characters occurring in matches within a 
	* multi-alignment block including identical columns (replaces
	* intervalSigma in old code) 
	* FIXME: alternate description
	* number of distinct residues in all columns in the extent of the 
	* alignment corresponding to a position
	*/
	int seq_num=(int)aligned_seqs.size();
	vector <double> row_sigma (seq_num,0.0);
	norm_seq_weights.resize(seq_num);

	//---- init posDistrib ---//
	for(k=0;k<=BLAST_AASIZE;k++)posDistrib[position*(BLAST_AASIZE+1)+k]=0;

	//------ iterate against aligned block ------//
	int wsfuck=(aligned_blocks_right[position]-aligned_blocks_left[position]+1);
	int wsfuck_K=BLOCK_LEN;  // default : 5
	double wsfuck_ratio=1.0*(wsfuck-1)/(wsfuck_K-1);
	if(wsfuck_ratio>1)wsfuck_ratio=1;
	for (i = (int)aligned_blocks_left[position]; i <= (int)aligned_blocks_right[position]; i++) 
	{
		/* number of distinct residues found in a column of the alignment
		* extent correspoding to a query position */
		int num_local_std_letters = 0; 
		//---- init column_residue_num ---//
		for(k=0;k<BLAST_AASIZE;k++)column_residue_num[k]=0;
	
		/* Count number of residues in column i of the alignment extent
		* corresponding to position */
		for (asi = 0; asi < seq_num; asi++)
		{
			int kSeqIdx = aligned_seqs[asi];
			int kResidue = BLAST_AA_To_ARND_Mapping(MSA_in[kSeqIdx][i]); //-> must return 0->19
			if(kResidue==-1)
			{
				fprintf(stderr,"bad residue %c at input seq %d \n",MSA_in[kSeqIdx][i],kSeqIdx);
				continue;
			}
			if (column_residue_num[kResidue] == 0) num_local_std_letters++;
			column_residue_num[kResidue]++;
		}
		sigma+=num_local_std_letters;
	
		//-> calculate sigma
		posDistrib[position*(BLAST_AASIZE+1)+num_local_std_letters]++; 
		if (num_local_std_letters > 1) 
		{
			/* num_local_std_letters == 1 means that all residues in
			* that column of the alignment extent are the same residue */
			distinct_residues_found = 1;
		}
		
		/* Calculate row_sigma, an intermediate value to calculate the
		* normalized sequence weights */
		for (asi = 0; asi < seq_num; asi++)
		{
			int kSeqIdx = aligned_seqs[asi];
			int kResidue = BLAST_AA_To_ARND_Mapping(MSA_in[kSeqIdx][i]); //-> must return 0->19
			if(kResidue==-1)
			{
				fprintf(stderr,"bad residue %c at input seq %d \n",MSA_in[kSeqIdx][i],kSeqIdx);
				continue;
			}
			/* This is a modified version of the Henikoff's idea in
			* "Position-based sequence weights" paper. The modification
			* consists in using the alignment extents. */
			row_sigma[asi] += (1.0 / ( pow((double)(column_residue_num[kResidue]),wsfuck_ratio) * num_local_std_letters) );
		}
	}
	posDistinct[position]=sigma;
	
	//------- final process --------//	
	if (distinct_residues_found) 
	{
		//--- calculate normalization factor ----//
		double weight_sum = 0.0;
		for (asi = 0; asi < seq_num; asi++) 
		{
			norm_seq_weights[asi] = row_sigma[asi] / (aligned_blocks_right[position]-aligned_blocks_left[position]+1);
			weight_sum += norm_seq_weights[asi];
		}
		//-- Normalize --//
		for (asi = 0; asi < seq_num; asi++) 
		{
			norm_seq_weights[asi] /= weight_sum;
		}
	} 
	else 
	{
		/* All residues in the extent of this position's alignment are the same
		* residue, therefore we distribute the sequence weight equally among
		* all participating sequences */
		for (asi = 0; asi < seq_num; asi++)
		{
			norm_seq_weights[asi] = (1.0/(double) seq_num);
		}
	}
}


//---- vice functions to compare two aligned_seqs
int BLAST_PSSM::Compare_Aligned_Seqs(vector <int> &in1, vector <int> &in2)
{
	int s1=(int)in1.size();
	int s2=(int)in2.size();
	if(s1!=s2)return 0;  //-> not equal
	int i;
	for(i=0;i<s1;i++)
	{
		if(in1[i]!=in2[i])return 0; //not equal
	}
	return 1;  //equal 
}

//---- vice functions to copy two aligned_seqs, in1 <- in2
int BLAST_PSSM::Copy_Aligned_Seqs(vector <int> &in1, vector <int> &in2)
{
	in1.clear();
	int size=(int)in2.size();
	int i;
	for(i=0;i<size;i++)in1.push_back(in2[i]);
	return size;
}

//---- vice functions to calculate aligned_seqs from MSA
int BLAST_PSSM::Calc_Aligned_Seqs(
	vector <string> &MSA_in,int pos,                         //-> input, canonical MSA
	vector <int> &aligned_seqs)                              //-> output, aligned sequences for a given position
{
	int k;
	int totnum=(int)MSA_in.size();
	int num=0;
	aligned_seqs.clear();
	for(k=0;k<totnum;k++)
	{
		if(MSA_in[k][pos]!='-')
		{
			aligned_seqs.push_back(k);
			num++;
		}
	}
	return num;
}

//----- vice function to check sequence weights ----//
int BLAST_PSSM::CheckSequenceWeights(
	int query_len, int *posNum,                              //-> input
	double *match_weights)                                   //-> input
{
	int i,k;
	for(i=0;i<query_len;i++)
	{
		if(posNum[i]<=1)continue;
		double sum=0.0;
		for(k=0;k<BLAST_AASIZE;k++) sum+=match_weights[i*BLAST_AASIZE+k];
		if(sum<0.99 || sum>1.01)return 0; //failed
	}
	return 1; //success
}


//=====================  main calculation steps ========================//

/** Main function to compute aligned blocks' properties for each position 
 * within multiple alignment (stage 3) 
 * Corresponds to posit.c:posComputeExtents
 * @param msa multiple sequence alignment data structure [in]
 * @param aligned_block data structure describing the aligned blocks'
 * properties for each position of the multiple sequence alignment [out]
 * @return PSIERR_BADPARAM if arguments are NULL
 *         PSI_SUCCESS otherwise
 */
void BLAST_PSSM::Compute_Alignment_Blocks(
	vector <string> &MSA_in,                                 //-> input, canonical MSA
	int *temp_block_left, int *temp_block_right,             //-> input, temporary aligned block
	int *aligned_blocks_left, int *aligned_blocks_right)     //-> output, final aligned_block
{
	int i,k;
	int totnum=(int)MSA_in.size();
	int totlen=(int)MSA_in[0].length();

	//init aligned block
	for(i=0;i<totlen;i++)
	{
		aligned_blocks_left[i]=-1;
		aligned_blocks_right[i]=totlen;
	}

	//calculate aligned block (start from the 1st aligned sequence)
	for(k=1;k<totnum;k++)
	{
		//init temporary aligned block for currrent sequence
		{
			for(i=0;i<totlen;i++)
			{
				temp_block_left[i]=-1;
				temp_block_right[i]=totlen;
			}
		}

		//extends left
		{
			int prev = 0;
			int curr = 0;
			if( MSA_in[k][prev]!='-' ) temp_block_left[prev]=prev;
			for (curr = prev + 1; curr < totlen; curr++, prev++) 
			{
				if ( MSA_in[k][curr]=='-') continue;
				if ( MSA_in[k][prev]!='-') temp_block_left[curr]=temp_block_left[prev];
				else temp_block_left[curr]=curr;
			}
		}

		//extends right
		{
			int last = totlen-1;
			int curr_ =0;
			if( MSA_in[k][last]!='-' ) temp_block_right[last]=last;
			for (curr_ = last - 1; curr_ >=0; curr_--, last--) 
			{
				if ( MSA_in[k][curr_]=='-') continue;
				if ( MSA_in[k][last]!='-') temp_block_right[curr_]=temp_block_right[last];
				else temp_block_right[curr_]=curr_;
			}
		}

		//assign current sequence to aligned block
		{
			for(i=0;i<totlen;i++)
			{
				if( MSA_in[k][i]!='-')
				{
					aligned_blocks_left[i]=aligned_blocks_left[i]>temp_block_left[i]?aligned_blocks_left[i]:temp_block_left[i];
					aligned_blocks_right[i]=aligned_blocks_right[i]<temp_block_right[i]?aligned_blocks_right[i]:temp_block_right[i];
				}
			}
		}
	}

	//final check
	for(i=0;i<totlen;i++)
	{
		if(aligned_blocks_left[i]==-1 || aligned_blocks_right[i]==totlen)
		{
			aligned_blocks_left[i]=i;
			aligned_blocks_right[i]=i;
		}
	}

/*
//wstest
for(i=0;i<totlen;i++)
{
	printf("%4d %4d\n",aligned_blocks_left[i],aligned_blocks_right[i]);
}
*/

}


/** Main function to calculate the sequence weights. Should be called with the
 * return value of PSIComputeAlignmentBlocks (stage 4)
 * Corresponds to posit.c:posComputeSequenceWeights
 * @param msa multiple sequence alignment data structure [in]
 * @param aligned_blocks data structure describing the aligned blocks'
 * properties for each position of the multiple sequence alignment [in]
 * @param nsg_compatibility_mode set to true to emulate the structure group's
 * use of PSSM engine in the cddumper application. By default should be FALSE
 * [in]
 * @param seq_weights data structure containing the data needed to compute the
 * sequence weights [out]
 * @return PSIERR_BADPARAM if arguments are NULL, PSIERR_OUTOFMEM in case of
 * memory allocation failure, PSIERR_BADSEQWEIGHTS if the sequence weights fail
 * to add up to 1.0, PSI_SUCCESS otherwise
 */
int BLAST_PSSM::Compute_Sequence_Weights(
	vector <string> &MSA_in,                                 //-> input, canonical MSA
	int *aligned_blocks_left, int *aligned_blocks_right,     //-> input, aligned block
	double *pseudo_expno, int *column_residue_num,           //-> input, table of expectations
	int *posDistrib, int *posDistinct, int *posNum,          //-> input, temporary for pseudo_count
	double *match_weights, double *PosObserve,               //-> output, match_weight and pseudo_count
	double *column_weights)                                  //-> output, column_weight
{
	int pos = 0;                  /* position index */
	int asi = 0;                  /* sequence index */
	int num;
	int totlen=(int)MSA_in[0].length();
	vector <int> aligned_seqs;
	vector <int> aligned_seqs_prev;
	vector <double> Norm_Seq_Weights;

	//-- init match_weights ---//
	for(pos=0;pos<totlen*BLAST_AASIZE;pos++)match_weights[pos]=0.0;
	for(pos=0;pos<totlen;pos++)column_weights[pos]=0.0;
	for(pos=0;pos<totlen;pos++)PosObserve[pos]=0.0;

	//-- calculate each position ---//
	aligned_seqs.clear();
	aligned_seqs_prev.clear();
	for (pos = 0; pos < totlen; pos++) 
	{
		//copy current to prev
		Copy_Aligned_Seqs(aligned_seqs_prev,aligned_seqs);
		//get current aligned sequences
		num=Calc_Aligned_Seqs(MSA_in,pos,aligned_seqs);
		//ignore positions of no matches //
		posNum[pos]=num;
//		if(num<=1)continue;

		//compare current to prev
		if( !Compare_Aligned_Seqs(aligned_seqs,aligned_seqs_prev)) //-> not equal !!
		{
			//calculate posDistinctDistrib
			Compute_SeqPos_Weights(MSA_in,column_residue_num,
				aligned_blocks_left,aligned_blocks_right,
				pos,aligned_seqs,posDistrib,posDistinct,Norm_Seq_Weights);

			//calculate EffectiveObservations
			double effect_observe=Effective_Observations(
				aligned_blocks_left,aligned_blocks_right,
				pos,pseudo_expno,posDistrib,posNum);
			PosObserve[pos]=effect_observe;
		}
		else                                                       //-> equal !!
		{
			int index;
			for(index=0;index<=BLAST_AASIZE;index++)
				posDistrib[pos*(BLAST_AASIZE+1)+index]=posDistrib[(pos-1)*(BLAST_AASIZE+1)+index];
			posDistinct[pos]=posDistinct[pos-1];
			PosObserve[pos]=PosObserve[pos-1];
			/* seq_weights->norm_seq_weights are unchanged from the previous
			* iteration, leaving them ready to be used in
			* _PSICalculateMatchWeights */
		}

		// Uses seq_weights->norm_seq_weights to populate match_weights //
		for (asi = 0; asi < num; asi++) 
		{
			int kSeqIdx = aligned_seqs[asi];
			int kResidue = BLAST_AA_To_ARND_Mapping(MSA_in[kSeqIdx][pos]); //-> must return 0->19
			if(kResidue==-1)
			{
				fprintf(stderr,"bad residue %c at input seq %d \n",MSA_in[kSeqIdx][pos],kSeqIdx);
				continue;
			}
			match_weights[pos*BLAST_AASIZE+kResidue]+= Norm_Seq_Weights[asi];
			column_weights[pos]+=Norm_Seq_Weights[asi];
		}

/*
//wstest
for(int r=0;r<BLAST_AASIZE;r++)
{
	int www=BLAST_ACDE_To_ARND_Mapping(r);
	printf("%e ",match_weights[pos*BLAST_AASIZE+www]);
}
printf("\n");
*/

	}

	//-- check match weights ---//
	int retv=CheckSequenceWeights(totlen,posNum,match_weights);
	if(retv!=1)
	{
		fprintf(stderr,"sequence weight error here !!! \n");
		return -1;
	}
	return 0;
}


/** Main function to compute the PSSM's frequency ratios (stage 5).
 * Implements formula 2 in Nucleic Acids Research, 2001, Vol 29, No 14.
 * Corresponds to posit.c:posComputePseudoFreqs
 * @param msa multiple sequence alignment data structure [in]
 * @param seq_weights data structure containing the data needed to compute the
 * sequence weights [in]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in]
 * @param aligned_blocks data structure describing the aligned blocks'
 * properties for each position of the multiple sequence alignment [in]
 * @param pseudo_count pseudo count constant [in]
 * @param nsg_compatibility_mode set to true to emulate the structure group's
 * use of PSSM engine in the cddumper application. By default should be FALSE
 * @param internal_pssm PSSM being computed [out]
 * @return PSIERR_BADPARAM if arguments are NULL, PSI_SUCCESS otherwise
 */
//[note1]: according to NAR 1997 Vol.25 No.17, Gapped BLAST and PSI-BLAST, page 3396, eq.5
// Qi = (alpha*fi + beta*gi) / (alpha + beta)
// where, alpha is observe_weight, beta is the psedo_weight,
// fi is the real calculated frequency, gi is the pseudocount frequency.
// In more details, gi = Sigma(j) fj/pj * pij, where pij is the background pair-probability
//[note2]: remember pij = fij*pi*pj, so gi=pi* Sigma(j) fj*fij,
// if we return frequency_ratio, i.e., Qi/pi = ( alpha*fi/pi + beta*Sigma(j)fj*fij ) / (alpha + beta)
// the above function is calculated in the following procedure.
void BLAST_PSSM::Compute_Freq_Ratio(
	int query_len,                                           //-> input, query length
	double *std_prob, double *tmp_prob,                      //-> input, standard frequency
	double *back_prob_single, double *back_freq_ratio,       //-> input, matrix based frequency
	double *match_weights, double *PosObserve,               //-> input, match weights and pseudo count
	double *freq_ratios,double *freq_entro,                  //-> output, for frequency ratio, entro 
	double *pseudocounts)                                    //-> output, for pseudocounts
{
	int p = 0;                    /* index on positions */
	int r = 0;                    /* index on residues */
	//-- calculate each position ---//
	for (p = 0; p < query_len; p++) 
	{
		double columnCounts = 0.0; /*column-specific pseudocounts*/
		double observations = 0.0; /* pre-calculated column observations */
		double pseudoWeight; /*multiplier for pseudocounts term*/

		//-- calculate pseudo count --//
		double *column_match_weights=match_weights+p*BLAST_AASIZE;
		observations=PosObserve[p];
		columnCounts=ColumnSpecific_Pseudocounts(column_match_weights,tmp_prob,back_prob_single,observations);

		//-- check for pseudo count
		if (columnCounts >= PSEUDO_MAX) 
		{
			pseudoWeight = kZeroObsPseudo;
			observations = 0;
		}
		else 
		{
			pseudoWeight = columnCounts;
		}
		pseudocounts[p]=pseudoWeight;

		//-- calculate each residue ---//
		for (r = 0; r < BLAST_AASIZE; r++)
		{
			int i = 0;             /* loop index */
			/* beta( Sum_j(f_j f_ij) ) in formula 2 */          
			/* Renamed to match the formula in the paper */
			double pseudo = 0.0;  
			double kBeta = pseudoWeight;
			double numerator = 0.0;         /* intermediate term */
			double denominator = 0.0;       /* intermediate term */
			double qOverPEstimate = 0.0;    /* intermediate term */

			/* As specified in 2001 paper, underlying matrix frequency ratios are used here */
			for (i = 0; i < BLAST_AASIZE; i++) 
			{
				pseudo += match_weights[p*BLAST_AASIZE+i] * back_freq_ratio[r*BLAST_AASIZE+i];
			}
			pseudo *= kBeta;

			//-> calculate final pseudo count value
			numerator = (observations * match_weights[p*BLAST_AASIZE+r] / std_prob[r]) + pseudo;
			denominator = observations + kBeta;
			qOverPEstimate = numerator/denominator;

			/* Note artificial multiplication by standard probability to normalize */
			//-> this will lead to Qi, i.e., estimated frequency, instead of frequency ratio.
			freq_ratios[p*BLAST_AASIZE+r] = qOverPEstimate*std_prob[r];
		}

		//-- calculate entropy for each position ---//
		double wsr=0.0;
		for(r = 0; r < BLAST_AASIZE; r++)
		{
		        if ( freq_ratios[p*BLAST_AASIZE+r] > kPosEpsilon)
                	wsr += freq_ratios[p*BLAST_AASIZE+r] * log (freq_ratios[p*BLAST_AASIZE+r]/std_prob[r]) / NCBIMATH_LN2;
		}
		freq_entro[p]=wsr;

/*
//wstest
for(r=0;r<BLAST_AASIZE;r++)
{
	int pos=BLAST_ACDE_To_ARND_Mapping(r);
	printf("%e ",freq_ratios[p*BLAST_AASIZE+pos]);
}
printf("\n");
*/

	}
}


//--------- a small step to calculate column specific weights that considers pseudocounts --------------//
void BLAST_PSSM::Calc_Column_Weight_P(
	int query_len, double *column_weights,                //-> input, original column weight without pseudo
	int *aligned_blocks_left, int *aligned_blocks_right,  //-> input, temporary aligned block
	int *posDistinct,double *pseudocounts,                //-> input, pseudo count related
	double *column_weights_p)                             //-> output, column weight considering pseudo
{
	int i;
	for(i=0;i<query_len;i++)
	{
		if(pseudocounts[i]<=0.0)column_weights_p[i]=0.0;
		else
		{
			double gapless_column_weights = 1.0*column_weights[i] / pseudocounts[i];
			int aligned_block_size = aligned_blocks_right[i]-aligned_blocks_left[i]+1;
			if(aligned_block_size==0)column_weights_p[i]=0.0;
			else column_weights_p[i] = gapless_column_weights * ( 1.0 * posDistinct[i] / aligned_block_size - 1);
		}
	}
}


/** Converts the PSSM's frequency ratios obtained in the previous stage to a 
 * PSSM of scores. (stage 6) 
 * @param internal_pssm PSSM being computed [in|out]
 * @param query query sequence in ncbistdaa encoding. The length of this
 * sequence is read from internal_pssm->ncols [in]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in]
 * @param std_probs array containing the standard residue probabilities [in]
 * @return PSIERR_BADPARAM if arguments are NULL, PSI_SUCCESS otherwise
 */
int BLAST_PSSM::Compute_PSSM_Matrix(
	string &query,                                           //-> input, query sequence
	double *std_prob, double *back_freq_ratio,               //-> input, matrix based frequency
	double *freq_ratios, double ideal_lambda,                //-> input, previously calculated frequency 
	int *PSSM_ini)                                           //-> output, for initial PSSM (un-scaled)
{
	int i = 0;
	int j = 0;
	int query_len=(int)query.length();
	
	/* Each column is a position in the query */
	for (i = 0; i < query_len; i++) 
	{
		int kResidue = BLAST_AA_To_ARND_Mapping(query[i]); //-> must return 0->19
		if(kResidue==-1)
		{
			fprintf(stderr,"bad residue %c at query seq \n",query[i]);
			return -1;
		}

		//-- calculate each residue --//
		for (j = 0; j < BLAST_AASIZE; j++) 
		{
			/* Division compensates for multiplication in _PSIComputeFreqRatios */
			double qOverPEstimate = freq_ratios[i*BLAST_AASIZE+j] / std_prob[j];
			PSSM_ini[i*BLAST_AASIZE+j] = BLAST_Nint (kPSIScaleFactor*log(qOverPEstimate)/ideal_lambda);
		}
	}

	//return success
	return 0;
}


/** Scales the PSSM (stage 7)
 * @param query query sequence in ncbistdaa encoding. The length of this
 * sequence is read from internal_pssm->ncols [in]
 * @param std_probs array containing the standard background residue 
 * probabilities [in]
 * @param internal_pssm PSSM being computed [in|out]
 * @param sbp score block structure initialized for the scoring system used
 * with the query sequence [in|out]
 * @return PSIERR_BADPARAM if arguments are NULL, PSIERR_POSITIVEAVGSCORE if
 * the average score of the generated PSSM is positive, PSI_SUCCESS otherwise
 */
int BLAST_PSSM::Scaling_PSSM_Matrix(
	int query_len,                                           //-> input, query length
	double ideal_lambda,double gap_std_K, double ideal_K,    //-> input, KARLIN parameters
	double *std_prob,                                        //-> input, standard background probability
	int *PSSM_ini,                                           //-> input, initial PSSM (un-scaled)
	int *PSSM_fin,                                           //-> output, PSSM with size L*20
	double &lambda,double &H,double &K,double &gap_K)        //-> output, KARLIN parameters
{
	int retv;
	int first_time = 1;  //TRUE
	int too_high = 1;    //TRUE
	int index = 0;     /* loop index */
	double factor;
	double factor_low = 1.0;
	double factor_high = 1.0;
	double new_lambda = 0.0;        /* Karlin-Altschul parameter calculated from scaled matrix*/
	double cur_scale_factor = kPSIScaleFactor;
	
	//-- determine factor_low and factor_high
	factor = 1.0;
	for ( ; ; ) 
	{
		int i = 0;
		int j = 0;

		//calculate PSSM_fin
		for (i = 0; i < query_len; i++) 
		{
			for (j = 0; j < BLAST_AASIZE; j++) 
			{
				PSSM_fin[i*BLAST_AASIZE+j] = 
					BLAST_Nint(factor*PSSM_ini[i*BLAST_AASIZE+j]/cur_scale_factor);
			}
		}

		//update LambdaK
		{
			retv=Calc_PSSM_LambdaK(PSSM_fin,query_len,std_prob,gap_std_K,ideal_K,lambda,H,K,gap_K);
			if(retv!=0)
			{
				fprintf(stderr,"update_lambda error !! \n");
				return -1;
			}
			new_lambda=lambda;
		}

		//judge terminate condition to determine factor low/high
		if (new_lambda > ideal_lambda) 
		{
			if (first_time) 
			{
				factor_high = 1.0 + kPositScalingPercent;
				factor = factor_high;
				factor_low = 1.0;
				too_high = 1;
				first_time = 0;
			} 
			else 
			{
				if (too_high == 0) break;
				factor_high += (factor_high - 1.0);
				factor = factor_high;
			}
		} 
		else if (new_lambda > 0) 
		{
			if (first_time) 
			{
				factor_high = 1.0;
				factor_low = 1.0 - kPositScalingPercent;
				factor = factor_low;
				too_high = 0;
				first_time = 0;
			} 
			else 
			{
				if (too_high == 1) break;
				factor_low += (factor_low - 1.0);
				factor = factor_low;
			}
		} 
		else 
		{
			fprintf(stderr,"new_lambda %f less than 0. \n",new_lambda);
			return -1;
		}
	}

	/* Binary search for kPositScalingNumIterations times */
	for (index = 0; index < kPositScalingNumIterations; index++) 
	{
		int i = 0;
		int j = 0;
		factor = (factor_high + factor_low)/2;

		//calculate PSSM_fin
		for (i = 0; i < query_len; i++) 
		{
			for (j = 0; j < BLAST_AASIZE; j++) 
			{
				PSSM_fin[i*BLAST_AASIZE+j] = 
					BLAST_Nint(factor*PSSM_ini[i*BLAST_AASIZE+j]/cur_scale_factor);
			}
		}

		//update LambdaK
		{
			retv=Calc_PSSM_LambdaK(PSSM_fin,query_len,std_prob,gap_std_K,ideal_K,lambda,H,K,gap_K);
			if(retv!=0)
			{
				fprintf(stderr,"update_lambda error !! \n");
				return -1;
			}
			new_lambda=lambda;
		}

		//update factor
		if (new_lambda > ideal_lambda) factor_low = factor;
		else factor_high = factor;
	}
	
	//return
	return 0;
}

//===================== final process: MSA_To_PSSM =====================//
//-> Given MSA, calculate PSSM
//[note]: the MSA input must be valid and upper case.
int BLAST_PSSM::MSA_To_PSSM(
	vector <string> &MSA_in,                                 //-> input, MSA be valid and upper case, with first be query
	int block_len)                                           //-> input, if block length < block_len, then use blastpgp
{
	//init
	BLOCK_LEN=block_len;
	int query_len=(int)MSA_in[0].length();
	int retv;

	//-> 1. Compute_Alignment_Blocks
	Compute_Alignment_Blocks(MSA_in,Temp_Block_Left,Temp_Block_Right,
		Aligned_Block_Left,Aligned_Block_Right);
	//-> 2. Compute Sequence Weights
	retv=Compute_Sequence_Weights(MSA_in,Aligned_Block_Left,Aligned_Block_Right,
		Pseudo_Expno,Column_Residue_Num,
		PosDistinctDistrib,PosDistinctNum,PosNumParticipating,
		Match_Weights,PosObservations,Colum_Weights);
	if(retv!=0)return -1;
	//-> 3. Compute Freq Ratios
	Compute_Freq_Ratio(query_len,Freq_Std,Inter_Prob,Freq_Back,Freq_Ratios,
		Match_Weights,PosObservations,Freq_Ratio,Freq_Entro,Pseudo_Prob);
	//-> 4.	Compute column specific weights considering pseudocounts
	Calc_Column_Weight_P(query_len,Colum_Weights,Aligned_Block_Left,Aligned_Block_Right,
		PosDistinctNum,Pseudo_Prob,Colum_Weights_P);
	//-> 5. Compute initial PSSM Matrix
	retv=Compute_PSSM_Matrix(MSA_in[0],Freq_Std,Freq_Ratios,
		Freq_Ratio,Ideal_Lambda,PSSM_Ini);
	if(retv!=0)return -1;
	//-> 6. Scaling PSSM Matrix
	retv=Scaling_PSSM_Matrix(query_len,Ideal_Lambda,Gap_K,Ideal_K,Freq_Std,
		PSSM_Ini,PSSM_Fin,Cur_Lambda,Cur_H,Cur_K,Cur_Gap_K);
	if(retv!=0)return -1;

	//final return
	return 0;
}
