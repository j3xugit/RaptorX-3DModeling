#pragma once
#include "blast_stat.h"


//====class: BLAST_PSSM ====//
//=> calculate BLAST PSSM
class BLAST_PSSM : public BLAST_Stat
{
public:
	BLAST_PSSM(void);
	~BLAST_PSSM(void);

//----------- macros -----------------//
public:
	int BLOCK_LEN;             // 5   /* for block length less than BLOCK_LEN, use blasgpgp */
	int MAX_IND_OBSERVATIONS;  // 400 /**< max number of independent observation for pseudocount calculation */
	int PSEUDO_MAX;            // 1000000 /**< effective infinity */
	double kZeroObsPseudo;     // 30.0 /*arbitrary constant to use for columns with zero observations */
	int kPSIScaleFactor;       // 200 /**< PSSM scale in initial stage */
	int kPositScalingNumIterations; // 10;
	double kPositScalingPercent;    // 0.05
	// used in ColumnSpecific_Pseudocounts
	double kPseudoMult;         // 500.0;   /*Was PSEUDO_MULTIPLIER */
	double kPseudoNumerator;    // 0.0457;  /*numerator of entropy-based method, was PSEUDO_NUMERATOR */
	double kPseudoExponent;     // 0.8;     /*exponent of denominator, was PSEUDO_EXPONENT */
	double kPseudoSmallInitial; // 5.5;     /*small number of pseudocounts to avoid 0 probabilities */
	double kPosEpsilon;         // 0.0001;  /*minimum return value of s_computeRelativeEntropy */

//----------- variables --------------//
public:
	//---- query length  ---------------------//
	int Query_Length;            //-> query length, L
	double Cur_Lambda;           //-> current lambda
	double Cur_H;                //-> current H
	double Cur_K;                //-> current K
	double Cur_Gap_K;            //-> current gap lambda
	//---- aligned block related -------------//
	int *Aligned_Block_Left;     //-> aligned block left, size is L
	int *Aligned_Block_Right;    //-> aligned block right, size is L
	int *Temp_Block_Left;        //-> temporary aligned block left, size is L
	int *Temp_Block_Right;       //-> temporary aligned block right, size is L
	//---- pseudo count related --------------//
	int *PosDistinctDistrib;     //-> /*number of distinct letters at each position*/, size is L*21 !!
	int *PosDistinctNum;         //-> column specific number of distinct residues, size is L
	int *PosNumParticipating;    //-> /*number of sequences at each position*/, size is L
	double *Pseudo_Prob;         //-> column specific pseudo probability, size is L
	//---- sequence weight related -----------//
	double *PosObservations;     //-> /*number of expected observations at each position*/, size is L
	double *Match_Weights;       //-> match state sequence weights, size is L*20
	double *Colum_Weights;       //-> column specific weight, sizie is L
	double *Colum_Weights_P;     //-> column specific weight considering pseudo, size is L
	//---- PSSM related ----------------------//
	double *Freq_Ratio;          //-> frequency ratio at each position, size is L*20 (double)
	double *Freq_Entro;          //-> the entropy of frequency ratio at each position, sizie is L (double)
	int *PSSM_Ini;               //-> initial or temporary PSSM, size is L*20 (double)
	int *PSSM_Fin;               //-> final output PSSM, size is L*20 (int)
	//---- temporary data structure ----------//
	double *Pseudo_Expno;        //-> table of expectations, size is 400
	double *Inter_Prob;          //-> intermediate probability, size is 20
	int *Column_Residue_Num;     //-> intermediate column residue number count, size is 20

//----------- functions --------------//
public:
	//---- data create/delete ----//
	void Init_Pseudo_Expno(void);
	void Create_DataStructure(int query_len);
	void Delete_DataStructure(void);

	//---- pseudo count related functions ----//
	void Init_ExpNum_Observations(
		double *back_prob, double *expno);                     //-> input and output
	double Effective_Observations(
		int *align_blk_left,int *align_blk_right,              //-> input
		int position, double *pseudo_expno,                    //-> input
		int *posDistrib, int *posNum);                         //-> input
	double ColumnSpecific_Pseudocounts(                      
		double *input_prob,                                    //-> input, column specific prob, size is 20
		double *inter_prob,                                    //-> input, just an intermediate term, size is 20
		double *back_prob,                                     //-> input, matrix-specific back prob, size is 20
		double observations);                                  //-> input, an estimate of observed residues

	//---- vice calculation functions -----//
	//-> [2.1] Compute Position Sequence Weights
	//[note]: applied in [2]
	void Compute_SeqPos_Weights(
		vector <string> &MSA_in, int *column_residue_num,      //-> input
		int *aligned_blocks_left, int *aligned_blocks_right,   //-> input, aligned block
		int position, vector <int> &aligned_seqs,              //-> input, aligned sequences for the given position
		int *posDistrib, int *posDistinct,                     //-> output for positional distribute and distinct
		vector <double> &norm_seq_weights);                    //-> output for positional match_weight
	//-> [2.2] aligned_seqs operations
	int Compare_Aligned_Seqs(                                      //-> return same or not
		vector <int> &in1, vector <int> &in2);                 //-> input, two aligned_seqs for comparison
	int Copy_Aligned_Seqs(
		vector <int> &in1, vector <int> &in2);                 //-> input, two aligned_seqs for comparison
	int Calc_Aligned_Seqs(
		vector <string> &MSA_in,int pos,                       //-> input, canonical MSA
		vector <int> &aligned_seqs);                           //-> output, aligned sequences for a given position
	//-> [2.3] CheckSequenceWeights
	int CheckSequenceWeights(                                      //-> return valid or not
		int query_len, int *posNum,                            //-> input
		double *match_weights);                                //-> input

	//---- main calculation steps ----//
	//-> [1] Compute Alignment Blocks
	//[note]: aligned_blocks length is L, 
	//        aligned_blocks_size = aligned_blocks_right - aligned_blocks_left + 1
	void Compute_Alignment_Blocks(
		vector <string> &MSA_in,                               //-> input, canonical MSA
		int *temp_block_left, int *temp_block_right,           //-> input, temporary aligned block
		int *aligned_blocks_left, int *aligned_blocks_right);  //-> output for aligned block

	//-> [2] Compute Sequence Weights
	//[note]: 
	//----- for match_weight -----//
	//1. match_weights size is L*20 (major data_structure for Freq_Ratio, 
	//      record AA frequence ratio for each position. )
	//----- for pseudo_count ----//
	//2. posDistinctDistrib size is L*21 (major data_structure for pseudocount, 
	//      record AA appearing number for each position, so the maximal number should be 21)
	//3. posNumParticipating size is L (minor data_structure for pseudocount, record sequence number for each position)
	//4. posObserve size is L (minor data_structure for pseudocount, record the expected observation)
	int Compute_Sequence_Weights(                                  //-> return success or not
		vector <string> &MSA_in,                               //-> input, canonical MSA
		int *aligned_blocks_left, int *aligned_blocks_right,   //-> input, aligned block
		double *pseudo_expno, int *column_residue_num,         //-> input, table of expectations
		int *posDistrib, int *posDistinct,int *posNum,         //-> input, temporary for pseudo_count
		double *match_weights, double *PosObserve,             //-> output for match_weight and pseudo_count
		double *column_weights);                               //-> output for column_weights

	//-> [3] Compute Freq Ratios
	//[note]: back_prob_single is the matrix dependent amino acid probability, size is 20
	//        std_prob is the background amino acid probability, size is 20
	//        back_freq_ratio is the matrix dependent pairwise amino acid frequency ratio, size is 20*20
	void Compute_Freq_Ratio(
		int query_len,                                         //-> input, query length
		double *std_prob, double *tmp_prob,                    //-> input, standard frequency
		double *back_prob_single, double *back_freq_ratio,     //-> input, matrix based frequency
		double *match_weights, double *PosObserve,             //-> input, match_weight and pseudo_count
		double *freq_ratios,double *freq_entro,                //-> output for frequency ratio, entro and pseudo
		double *pseudocounts);                                 //-> output for pseudocount

	//-> [4] Calculate Column_Weight_P
	void Calc_Column_Weight_P(
		int query_len, double *column_weights,                 //-> input, original column weight without pseudo
		int *aligned_blocks_left, int *aligned_blocks_right,   //-> input, aligned block
		int *posDistinct,double *pseudocounts,                 //-> input, pseudo count related
		double *column_weights_p);                             //-> output, column weight considering pseudo

	//-> [5] Compute PSSM Matrix (initial PSSM, before scaling)
	int Compute_PSSM_Matrix(                                       //-> return success or not
		string &query,                                         //-> input, query sequence
		double *std_prob, double *back_freq_ratio,             //-> input, matrix based frequency
		double *freq_ratios,double ideal_lambda,               //-> input, previously calculated frequency 
		int *PSSM_ini);                                        //-> output for initial PSSM (un-scaled)

	//-> [6] Scaling PSSM Matrix
	int Scaling_PSSM_Matrix(                                       //-> return success or not
		int query_len,                                         //-> input, query length
		double ideal_lambda,double gap_std_K, double ideal_K,  //-> input, initial ideal_lambda
		double *std_prob,                                      //-> input, standard background probability
		int *PSSM_ini,                                         //-> input, un-scaled initial PSSM
		int *PSSM_fin,                                         //-> output for final PSSM (scaled)
		double &lambda,double &H,double &K,double &gap_K);     //-> output for the current KARLIN parameter

	//---------- final process -----------//
	//[note]: the MSA input must be valid and upper case.
	int MSA_To_PSSM(
		vector <string> &MSA_in,                               //-> input, MSA must be valid and upper case, with first be query
		int block_len);                                        //-> intpu, if block length < block_len, then use blastpgp

};
//====class: BLAST_PSSM====//over
