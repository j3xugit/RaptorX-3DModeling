#pragma once
#include "blast_util.h"
#include "blast_mat.h"


//====class: BLAST_Stat ====//
//=> calculate Karlin-Altschul parameters
class BLAST_Stat
{
public:
	BLAST_Stat(void);
	~BLAST_Stat(void);

//----------- macros -----------------//
public:
	int BLAST_AASIZE;     // should be 20  //-> the most imporatn macros, for AA size
	int BLAST_SCORE_MIN;  // should be -32767
	int BLAST_SCORE_MAX;  // should be 32767
	/****************************************************************************
	For more accuracy in the calculation of K, set K_SUMLIMIT to 0.00001.
	For high speed in the calculation of K, use a K_SUMLIMIT of 0.001
	Note:  statistical significance is often not greatly affected by the value
	of K, so high accuracy is generally unwarranted.
	*****************************************************************************/
	double BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT; //   (1.e-5) 
		/**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd */
	int BLAST_KARLIN_LAMBDA_ITER_DEFAULT;     //   17      
		/**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)*/
	double BLAST_KARLIN_LAMBDA0_DEFAULT;   // 0.5             
		/**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */
	double BLAST_KARLIN_K_SUMLIMIT_DEFAULT; //0.0001          
		/**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */
	int BLAST_KARLIN_K_ITER_MAX; // 100                    
		/**< upper limit on iterations for BlastKarlinLHtoK */

//----------- variables --------------//
public:
	// statistics related parameters
	//-> ideal value
	double Ideal_Lambda;  // Ideal values (only dependent on standard background )    //-> [fixed]
	double Ideal_H;       // Ideal H                                                  //-> [fixed]
	double Ideal_K;       // Ideal K                                                  //-> [fixed]
	//-> standard valud
	double Std_Lambda;    // Standard values (only dependent on query sequence )      //-> [fixed]
	double Std_H;         // Standard H                                               //-> [fixed]
	double Std_K;         // Standard K                                               //-> [fixed]
	//-> gapped value
	double Gap_Lambda;    // Gap values (only dependent on gap_open, gap_extend )     //-> [fixed]
	double Gap_H;         // Gap H                                                    //-> [fixed]
	double Gap_K;         // Gap K                                                    //-> [fixed]
	//-> gap parameter
	int Gap_Open;         // gap open penalty                                         //-> default is 11
	int Gap_Extend;       // gap extend penalty                                       //-> default is 1
	// substitution matrix
	string Name_Matrix;   // name of the matrix, (default is BLOSUM62)                //-> [fixed]
	int Bit_Scale_Factor; // Used to obtain scores in bit units,                      //-> [fixed]
	int *Subs_Matrix;     // score of the matrix, (size is 20*20)                     //-> [fixed]
	double *Freq_Ratios;  // frequency ratios of the matrix, (size is 20*20) //-> used for pseudocount, [fixed]
	double *Freq_Matrix;  // frequency of the matrix, (size is 20*20)        //-> used for pseudocount, [fixed]
	double *Freq_Back;    // matrix based background frequency, (size is 20) //-> used for pseudocount, [fixed]
	double *Freq_Std;     // standard background frequency, (size is 20)     //-> for all other usage,  [fixed]

//----------- functions --------------//
public:
	//-- matrix related operations ---//
	int Init_BLAST_Stat(string &name,int gap_open,int gap_extend);           //-> calculate all fixed data once
	int Init_QUERY_Stat(string &query);                                      //-> calculate standard KARLIN parameter
	int Blast_Get_Matrix(string &name);                                      //-> get all matrix data
	int Calc_DataStruc_MinMax(int *input,int size,
		int &ret_min,int &ret_max);

	//-- score frequency operations --//
	//-> calculate Lambda_K for PSSM
	int Calc_PSSM_LambdaK(
		int *PSSM, int query_len, double *std_prob,                            //-> input PSSM, query and standard prob
		double gap_std_K, double ideal_K,                                      //-> input ideal and standard gap K
		double &lambda, double &H, double &K,double &gap_K);                   //-> output current psi and gap KARLIN parameter 
	//-> calculate score frequency from PSSM
	void Calc_ScoreFreq_With_PSSM(
		int *pssm, int length,                //-> input, PSSM matrix (size is L*20)
		int sfp_ori_min, int sfp_ori_max,     //-> input, PSSM related max/min score
		double* std_prob,                     //-> input, background frequency
		int &sfp_obs_min, int &sfp_obs_max,   //-> output, observed max/min score
		double &score_aver,                   //-> output, score average
		double *sprob);                       //-> output, score frequency
	//-> calculate score frequency from matrix
	void Calc_ScoreFreq_With_Matrix(
		int *matrix,                          //-> input, substitution matrix
		int sfp_ori_min, int sfp_ori_max,     //-> input, matrix related max/min score
		double* rfp1, double* rfp2,           //-> input, frequency from 1st/2nd sequence, with 2nd fixed to background
		int &sfp_obs_min, int &sfp_obs_max,   //-> output, observed max/min score
		double &score_aver,                   //-> output, score average
		double *sprob);                       //-> output, score frequency

	//-- parameter calculation functions --//
	int Blast_KarlinBlkParameterCalc(double *sfp_sprob, double sfp_score_avg,
		int sfp_obs_min, int sfp_obs_max,double &lambda, double &H, double &K);
	double BlastKarlinLHtoK(double *sfp_sprob, double sfp_score_avg, 
		int sfp_obs_min, int sfp_obs_max,double lambda, double H,
		double sumlimit_, int iterlimit_);
	double NlmKarlinLambdaNR(double* probs, int d, int low, int high, 
		double lambda0, double tolx, int itmax, int maxNewton, int * itn );
};
//====class: BLAST_Stat====//over
