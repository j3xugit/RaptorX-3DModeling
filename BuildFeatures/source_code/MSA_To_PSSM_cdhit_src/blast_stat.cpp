#include "blast_stat.h"



/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
BLAST_Stat::BLAST_Stat(void)
{
	// assign macros
	BLAST_AASIZE=20;
	BLAST_SCORE_MIN=-32767;
	BLAST_SCORE_MAX=32767;
	// macros used in calculation for Larlin parameters
	BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT = (1.e-5); /**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd */
	BLAST_KARLIN_LAMBDA_ITER_DEFAULT = 17;          /**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)*/
	BLAST_KARLIN_LAMBDA0_DEFAULT = 0.5;             /**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */
	BLAST_KARLIN_K_SUMLIMIT_DEFAULT = 0.0001;       /**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */
	BLAST_KARLIN_K_ITER_MAX = 100;                  /**< upper limit on iterations for BlastKarlinLHtoK */
	// KARLIN variables
	//-> ideal
	Ideal_Lambda=-1;
	Ideal_H=-1;
	Ideal_K=-1;
	//-> standard
	Std_Lambda=-1;
	Std_H=-1;
	Std_K=-1;
	//-> gap
	Gap_Lambda=-1;
	Gap_H=-1;
	Gap_K=-1;
	// gap parameter
	Gap_Open=11;
	Gap_Extend=1;
	// create data structure
	Subs_Matrix=new int[BLAST_AASIZE*BLAST_AASIZE];
	Freq_Matrix=new double[BLAST_AASIZE*BLAST_AASIZE];
	Freq_Ratios=new double[BLAST_AASIZE*BLAST_AASIZE];
	Freq_Back=new double[BLAST_AASIZE];
	Freq_Std=new double[BLAST_AASIZE];
}

/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
BLAST_Stat::~BLAST_Stat(void)
{
	// delete data structure
	delete [] Subs_Matrix;
	delete [] Freq_Matrix;
	delete [] Freq_Ratios;
	delete [] Freq_Back;
	delete [] Freq_Std;
}


//=================  matrix related operations ====================//

//---- BLAST_Stat init data ------//
int BLAST_Stat::Init_BLAST_Stat(
	string &name,                            //-> input, matrix name
	int gap_open,int gap_extend)             //-> input, gap parameter
{
	//========== calculate matrix parameter =======//
	int retv;
	//[1] get matrix related data
	retv=Blast_Get_Matrix(name);
	if(retv!=0)return -1;
	//[2] calculate original max and min
	int score_ori_min=0;
	int score_ori_max=0;
	retv=Calc_DataStruc_MinMax(Subs_Matrix,BLAST_AASIZE*BLAST_AASIZE,
		score_ori_min,score_ori_max);
	if(retv!=0)return -1;
	//[3] create score probability
	double *score_prob0=new double[score_ori_max-score_ori_min+1]; //this is only temporary structure !!
	double *score_prob=score_prob0-score_ori_min;

	//========= calculate KARLIN parameter ========//
	double lambda,H,K;
	//[4] calc ideal parameter
	int score_obs_min, score_obs_max;
	double score_aver;
	Calc_ScoreFreq_With_Matrix(Subs_Matrix,score_ori_min,score_ori_max,
		Freq_Std,Freq_Std,score_obs_min,score_obs_max,score_aver,score_prob);
	retv=Blast_KarlinBlkParameterCalc(score_prob,score_aver,
		score_obs_min,score_obs_max,lambda,H,K);
	delete [] score_prob0;
	if(retv!=0)return -1;
	Ideal_Lambda=lambda;
	Ideal_H=H;
	Ideal_K=K;
	//[5] calc gap parameter
	Gap_Open=gap_open;
	Gap_Extend=gap_extend;
	retv=Blast_GetMatrixRelatedGapValues(lambda,H,K,gap_open,gap_extend,name.c_str());
	if(retv!=0)return -1;
	Gap_Lambda=lambda;
	Gap_H=H;
	Gap_K=K;

	//final return
	return 0;
}

//---- BLAST_Stat query KARLIN parameter ------//
int BLAST_Stat::Init_QUERY_Stat(string &query)
{
	//[1] calc query alphabet probability
	double *temp_prob=new double[BLAST_AASIZE];
	for(int i=0;i<BLAST_AASIZE;i++)temp_prob[i]=0;
	int sum=0;
	for(int i=0;i<(int)query.length();i++)
	{
		int code=BLAST_AA_To_ARND_Mapping(query[i]);
		if(code<0)continue;
		temp_prob[code]++;
		sum++;
	}
	if(sum==0)
	{
		fprintf(stderr,"query %s error !! \n",query.c_str());
		delete [] temp_prob;
		return -1;
	}
	for(int i=0;i<BLAST_AASIZE;i++)temp_prob[i]/=(double)sum;

	//[2] calculate original max and min
	int score_ori_min=0;
	int score_ori_max=0;
	int retv=Calc_DataStruc_MinMax(Subs_Matrix,BLAST_AASIZE*BLAST_AASIZE,
		score_ori_min,score_ori_max);

	//[3] create score probability
	double *score_prob0=new double[score_ori_max-score_ori_min+1]; //this is only temporary structure !!
	double *score_prob=score_prob0-score_ori_min;

	//[4] calc standard parameter
	double lambda,H,K;
	int score_obs_min, score_obs_max;
	double score_aver;
	Calc_ScoreFreq_With_Matrix(Subs_Matrix,score_ori_min,score_ori_max,
		temp_prob,Freq_Std,score_obs_min,score_obs_max,score_aver,score_prob);
	retv=Blast_KarlinBlkParameterCalc(score_prob,score_aver,
		score_obs_min,score_obs_max,lambda,H,K);
	delete [] temp_prob;
	delete [] score_prob0;
	if(retv!=0)return -1;

	//[5] assign standard parameter
	Std_Lambda=lambda;
	Std_H=H;
	Std_K=K;

	//final return
	return 0;
}

//--- Get matrix with a given name. Default should be BLOSUM62 -----//
int BLAST_Stat::Blast_Get_Matrix(string &name)
{
	//get matrix related data
	if(Blast_FrequencyDataIsAvailable(name.c_str())!=0)
	{
		fprintf(stderr,"Unsupported matrix %s !!",name.c_str());
		return -1;
	}
	if(Blast_GetSubstituteForMatrix(Subs_Matrix,name.c_str())!=0)
	{
		fprintf(stderr,"subs_matrix load error in matrix %s !!",name.c_str());
		return -1;
	}
	if(Blast_GetJointProbsForMatrix(Freq_Matrix,name.c_str())!=0)
	{
		fprintf(stderr,"freq_matrix load error in matrix %s !!",name.c_str());
		return -1;
	}
	if(Blast_GetFreqRatioForMatrix(Freq_Ratios,name.c_str())!=0)
	{
		fprintf(stderr,"freq_matrix load error in matrix %s !!",name.c_str());
		return -1;
	}
	if(Blast_GetMatrixBackgroundFreq(Freq_Back,name.c_str())!=0)
	{
		fprintf(stderr,"freq_back load error in matrix %s !!",name.c_str());
		return -1;
	}
	//get standard background frequency
	Blast_GetStandardBackgroundFreq(Freq_Std);
	//get name
	Name_Matrix=name;
	if(name!="BLOSUM45")Bit_Scale_Factor=2;
	else Bit_Scale_Factor=3;
	//normalize
	BLAST_FreqNormalize(Freq_Std,BLAST_AASIZE,1.0);
	BLAST_FreqNormalize(Freq_Back,BLAST_AASIZE,1.0);
	BLAST_FreqNormalize(Freq_Matrix,BLAST_AASIZE*BLAST_AASIZE,1.0);
	//return 0
	return 0;
}

//---- misc function ---//
int BLAST_Stat::Calc_DataStruc_MinMax(int *input,int size,
	int &ret_min,int &ret_max)
{
	//init
	int min_score = BLAST_SCORE_MAX;    /* minimum score in pssm */
	int max_score = BLAST_SCORE_MIN;    /* maximum score in pssm */
	//calc
	for(int i=0;i<size;i++)
	{
		int kScore = input[i];
		if (kScore <= BLAST_SCORE_MIN || kScore >= BLAST_SCORE_MAX) continue;
		max_score = kScore>max_score?kScore:max_score;
		min_score = kScore<min_score?kScore:min_score;
	}
	//check
	if(min_score==BLAST_SCORE_MAX)
	{
		fprintf(stderr,"min_score bad in Calc_DataStruc_MinMax !!");
		return -1;
	}
	if(max_score==BLAST_SCORE_MIN)
	{
		fprintf(stderr,"max_score bad in Calc_DataStruc_MinMax !!");
		return -1;
	}
	//return
	ret_min=min_score;
	ret_max=max_score;
	return 0;
}

//=================  score frequency calculation ====================//

//---- calculate Lambda_K for PSSM -----//
int BLAST_Stat::Calc_PSSM_LambdaK(
	int *PSSM, int query_len, double *std_prob,               //-> input PSSM, query and standard prob
	double gap_std_K, double ideal_K,                         //-> input ideal and standard gap K
	double &lambda, double &H, double &K,double &gap_K)       //-> output current psi and gap KARLIN parameter 
{
	int retv;
	//[1] calculate original max and min
	int score_ori_min, score_ori_max;
	retv=Calc_DataStruc_MinMax(PSSM,query_len*BLAST_AASIZE,
		score_ori_min,score_ori_max);
	if(retv!=0)return -1;
	//[2] create score probability
	double *score_prob0=new double[score_ori_max-score_ori_min+1]; //this is only temporary structure !!
	double *score_prob=score_prob0-score_ori_min;
	//[3] calc ideal parameter
	int score_obs_min, score_obs_max;
	double score_aver;
	Calc_ScoreFreq_With_PSSM(PSSM, query_len, score_ori_min,score_ori_max,
		std_prob,score_obs_min,score_obs_max,score_aver,score_prob);
	retv=Blast_KarlinBlkParameterCalc(score_prob,score_aver,
		score_obs_min,score_obs_max,lambda,H,K);
	//[4] calculate gap parameter
	gap_K=K*gap_std_K/ideal_K;
	//[4] delete score probability
	delete [] score_prob0;
	return retv;
}

/** Compute the probabilities for each score in the PSSM. */
void BLAST_Stat::Calc_ScoreFreq_With_PSSM(
	int *pssm, int length,                //-> input, PSSM matrix (size is L*20)
	int sfp_ori_min, int sfp_ori_max,     //-> input, PSSM related max/min score
	double* std_prob,                     //-> input, background frequency
	int &sfp_obs_min, int &sfp_obs_max,   //-> output, observed max/min score
	double &score_aver,                   //-> output, score average
	double *sprob)                        //-> output, score frequency
	                                      //-> Note: since score could be negative, this pointer must be shifted !!
{
	int  score;
	double score_sum;
	int  index1;
	int  index2;
	
	//-> 1. initialize score probability. Note, since score could be negative, this pointer must be shifted !!
	for (score = sfp_ori_min; score <= sfp_ori_max; score++) sprob[score] = 0.0;
	//-> 2. calculate score probability.
	for (index1=0; index1<length; index1++)
	{
		for (index2=0; index2<BLAST_AASIZE; index2++)
		{
			score = pssm[index1*BLAST_AASIZE+index2];
			if (score >= sfp_ori_min) 
			{
				sprob[score] += std_prob[index2] / length;
			}
		}
	}
	//-> 3. calculate sum of the probability.
	score_sum = 0.;
	sfp_obs_min = BLAST_SCORE_MIN;
	sfp_obs_max = BLAST_SCORE_MIN;
	for (score = sfp_ori_min; score <= sfp_ori_max; score++)
	{
		if (sprob[score] > 0.) 
		{
			score_sum += sprob[score];
			sfp_obs_max = score;
			if (sfp_obs_min == BLAST_SCORE_MIN) sfp_obs_min = score;
		}
	}
	//-> 4. calculate score average	
	score_aver = 0.0;
	if (score_sum > 0.0001)
	{
		for (score = sfp_obs_min; score <= sfp_obs_max; score++) 
		{
			sprob[score] /= score_sum;
			score_aver += score * sprob[score];
		}
	}
}


/** Calculates the score frequencies with substitution matrix */
void BLAST_Stat::Calc_ScoreFreq_With_Matrix(
	int *matrix,                          //-> input, substitution matrix
	int sfp_ori_min, int sfp_ori_max,     //-> input, matrix related max/min score
	double* rfp1, double* rfp2,           //-> input, frequency from 1st/2nd sequence, with 2nd fixed to background
	int &sfp_obs_min, int &sfp_obs_max,   //-> output, observed max/min score
	double &score_aver,                   //-> output, score average
	double *sprob)                        //-> output, score frequency
	                                      //-> Note: since score could be negative, this pointer must be shifted !!
{
	int  score;
	double score_sum;
	int  index1, index2;
	
	//-> 1. initialize score probability. Note, since score could be negative, this pointer must be shifted !!
	for (score = sfp_ori_min; score <= sfp_ori_max; score++) sprob[score] = 0.0;
	//-> 2. calculate score probability.
	for (index1=0; index1<BLAST_AASIZE; index1++)
	{
		for (index2=0; index2<BLAST_AASIZE; index2++)
		{
			score = matrix[index1*BLAST_AASIZE+index2];
			if (score >= sfp_ori_min) 
			{
				sprob[score] += rfp1[index1] * rfp2[index2];
			}
		}
	}
	//-> 3. calculate sum of the probability.
	score_sum = 0.;
	sfp_obs_min = BLAST_SCORE_MIN;
	sfp_obs_max = BLAST_SCORE_MIN;
	for (score = sfp_ori_min; score <= sfp_ori_max; score++)
	{
		if (sprob[score] > 0.) 
		{
			score_sum += sprob[score];
			sfp_obs_max = score;
			if (sfp_obs_min == BLAST_SCORE_MIN) sfp_obs_min = score;
		}
	}
	//-> 4. calculate score average	
	score_aver = 0.0;
	if (score_sum > 0.0001)
	{
		for (score = sfp_obs_min; score <= sfp_obs_max; score++) 
		{
			sprob[score] /= score_sum;
			score_aver += score * sprob[score];
		}
	}
}


//===================== major calculation functions ==================//


/** Computes the parameters lambda, H K for use in calculating the
 * statistical significance of high-scoring segments or subalignments (see
 * comment on blast_stat.c for more details).
 * @param kbp object containing Lambda, H, and K as well as scoring information [in|out]
 * @param sfp array of probabilities for all scores [in]
 * @return zero on success, 1 on error.
 */
int BLAST_Stat::Blast_KarlinBlkParameterCalc(
	double *sfp_sprob, double sfp_score_avg,   //-> input, score probability
	int sfp_obs_min, int sfp_obs_max,          //-> input, score range
	double &lambda, double &H, double &K)      //-> output, Karlin parameters
{
	int i;
	int d;
	int low=sfp_obs_min;   /* Lowest score (must be negative)  */
	int high=sfp_obs_max;  /* Highest score (must be positive) */
	double score_avg=sfp_score_avg;  /* Expected score must be negative */
	double *sprob=sfp_sprob;
	//ini check
	if(low>=0)
	{
		fprintf(stderr,"Lowest score %d must be negative \n",low);
		return -1;
	}
	if(high<=0)
	{
		fprintf(stderr,"Highest score %d must be positive \n",high);
		return -1;
	}
	if(score_avg>=0.)
	{
		fprintf(stderr,"Expected score %lf must be positive \n",score_avg);
		return -1;
	}

	//----- 1. Calculate the parameter Lambda -----//
	//-> 1.1 Find greatest common divisor of all scores */
	for (i = 1, d = -low; i <= high-low && d > 1; ++i) 
	{
		if (sprob[i+low] != 0.0)
		{
			d = BLAST_Gcd(d, i);
		}
	}
	//-> 1.2 Check d is positive or not
	if(d<=0)
	{
		fprintf(stderr,"Greatest common divisor %d must be positive \n",d);
		return -1;
	}
	//-> 1.3 Calculate Lambda
	int itn;
	lambda = NlmKarlinLambdaNR( sprob, d, low, high,
	                       BLAST_KARLIN_LAMBDA0_DEFAULT,
	                       BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
	                 20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn );
	//-> 1.4 Check lambda is negative or not
	if(lambda<0.)
	{
		fprintf(stderr,"Lambda %lf should be positive \n",lambda);
		return -1;
	}

	//----- 2. Calculate the parameter H -----//
	int score;
	double etonlam, sum, scale;
	etonlam = exp( - lambda );
	//-> 2.1 Calculate sum
	sum = low * sprob[low];
	for( score = low + 1; score <= high; score++ ) {
		sum = score * sprob[score] + etonlam * sum;
	}
	//-> 2.2 Calculate scale
	scale = BLAST_Powi( etonlam, high );
	if( scale > 0.0 ) 
	{
		H = lambda * sum/scale;
	} 
	else   /* Underflow of exp( -lambda * high ) */
	{ 
		H = lambda * exp( lambda * high + log(sum) );
	}
	//-> 2.3 Check H is negative or not
	if(H<0.)
	{
		fprintf(stderr,"H %lf should be positive \n",H);
		return -1;
	}

	//----- 3. Calculate the parameter K -----//
	//-> 3.1 Calculate K
	K = BlastKarlinLHtoK(sprob, score_avg, low, high, lambda, H,
		BLAST_KARLIN_K_SUMLIMIT_DEFAULT, BLAST_KARLIN_K_ITER_MAX);
	//-> 3.2 Check K is negative or not
	if(K<0.)
	{
		fprintf(stderr,"K %lf should be positive \n",K);
		return -1;
	}

	//return
	return 0;
}


/**
 * Find positive solution to 
 *
 *     sum_{i=low}^{high} exp(i lambda) * probs[i] = 1.
 * 
 * Note that this solution does not exist unless the average score is
 * negative and the largest score that occurs with nonzero probability
 * is positive.
 * 
 * @param probs         probabilities of a score occurring 
 * @param d             the gcd of the possible scores. This equals 1 if
 *                      the scores are not a lattice
 * @param low           the lowest possible score that occurs with
 *                      nonzero probability
 * @param high          the highest possible score that occurs with
 *                      nonzero probability.
 * @param lambda0       an initial guess for lambda
 * @param tolx          the tolerance to which lambda must be computed
 * @param itmax         the maximum number of times the function may be
 *                      evaluated
 * @param maxNewton     the maximum permissible number of Newton
 *                      iterations; after that the computation will proceed
 *                      by bisection.
 * @param *itn          the number of iterations needed to compute Lambda,
 *                      or itmax if Lambda could not be computed.
 *
 * Let phi(lambda) =  sum_{i=low}^{high} exp(i lambda) - 1. Then
 * phi(lambda) may be written
 *
 *     phi(lamdba) = exp(u lambda) f( exp(-lambda) )
 *
 * where f(x) is a polynomial that has exactly two zeros, one at x = 1
 * and one at x = exp(-lamdba).  It is simpler to solve this problem
 * in x = exp(-lambda) than it is to solve it in lambda, because we
 * know that for x, a solution lies in [0,1], and because Newton's
 * method is generally more stable and efficient for polynomials than
 * it is for exponentials.
 * 
 * For the most part, this function is a standard safeguarded Newton
 * iteration: define an interval of uncertainty [a,b] with f(a) > 0
 * and f(b) < 0 (except for the initial value b = 1, where f(b) = 0);
 * evaluate the function and use the sign of that value to shrink the
 * interval of uncertainty; compute a Newton step; and if the Newton
 * step suggests a point outside the interval of uncertainty or fails
 * to decrease the function sufficiently, then bisect.  There are
 * three further details needed to understand the algorithm:
 *
 * 1)  If y the unique solution in [0,1], then f is positive to the left of
 *     y, and negative to the right.  Therefore, we may determine whether
 *     the Newton step -f(x)/f'(x) is moving toward, or away from, y by
 *     examining the sign of f'(x).  If f'(x) >= 0, we bisect instead
 *     of taking the Newton step.
 * 2)  There is a neighborhood around x = 1 for which f'(x) >= 0, so
 *     (1) prevents convergence to x = 1 (and for a similar reason
 *     prevents convergence to x = 0, if the function is incorrectly
 *     called with probs[high] == 0).
 * 3)  Conditions like  fabs(p) < lambda_tolerance * x * (1-x) are used in
 *     convergence criteria because these values translate to a bound
 *     on the relative error in lambda.  This is proved in the
 *     "Blast Scoring Parameters" document that accompanies the BLAST
 *     code.
 *
 * The iteration on f(x) is robust and doesn't overflow; defining a
 * robust safeguarded Newton iteration on phi(lambda) that cannot
 * converge to lambda = 0 and that is protected against overflow is
 * more difficult.  So (despite the length of this comment) the Newton
 * iteration on f(x) is the simpler solution.
 */
double BLAST_Stat::NlmKarlinLambdaNR(
	double* probs, int d, int low, int high, double lambda0,  //-> input, score probability and range
	double tolx, int itmax, int maxNewton, int * itn )        //-> input, parameters
{
	int k;
	double x0, x, a = 0, b = 1;
	double f = 4;  /* Larger than any possible value of the poly in [0,1] */
	int isNewton = 0; /* we haven't yet taken a Newton step. */
	
	x0 = exp( -lambda0 );
	x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;
	
	for( k = 0; k < itmax; k++ ) /* all iteration indices k */
	{
		int i;
		double g, fold = f;
		int wasNewton = isNewton; /* If true, then the previous step was a */
		/* Newton step */
		isNewton  = 0;            /* Assume that this step is not */
	
		/* Horner's rule for evaluating a polynomial and its derivative */
		g = 0;
		f = probs[low];
		for( i = low + d; i < 0; i += d ) 
		{
			g = x * g + f;
			f = f * x + probs[i];
		}
		g = x * g + f;
		f = f * x + probs[0] - 1;
		for( i = d; i <= high; i += d ) 
		{
			g = x * g + f;
			f = f * x + probs[i];
		}
		/* End Horner's rule */

		//--- check for f ----//
		if( f > 0 ) 
		{
			a = x; /* move the left endpoint */
		} 
		else if( f < 0 ) 
		{ 
			b = x; /* move the right endpoint */
		} 
		else  /* f == 0 */
		{ 
			break; /* x is an exact solution */
		}
		if( b - a < 2 * a * ( 1 - b ) * tolx ) 
		{
			/* The midpoint of the interval converged */
			x = (a + b) / 2; 
			break;
		}

		//--- check for Newton step ---//
		if( k >= maxNewton ||
			/* If convergence of Newton's method appears to be failing; or */
			( wasNewton && fabs( f ) > .9 * fabs(fold) ) ||  
			/* if the previous iteration was a Newton step but didn't decrease 
			* f sufficiently; or */
			g >= 0 
			/* if a Newton step will move us away from the desired solution */
		) /* then bisect*/
		{ 
			x = (a + b)/2;
		} 
		else /* try a Newton step */
		{
			double p = - f/g;
			double y = x + p;
			if( y <= a || y >= b ) /* The proposed iterate is not in (a,b) */
			{ 
				x = (a + b)/2;
			} 
			else  /* The proposed iterate is in (a,b). Accept it. */
			{ 
				isNewton = 1;
				x = y;
				if( fabs( p ) < tolx * x * (1-x) ) break; /* Converged */
			} /* else the proposed iterate is in (a,b) */
		} /* else try a Newton step. */ 
	} /* end for all iteration indices k */
	*itn = k; 
	return -log(x)/d;
}



/** The following procedure computes K. The input includes Lambda, H,
 *  and an array of probabilities for each score.
 *  There are distinct closed form for three cases:
 *  1. high score is 1 low score is -1
 *  2. high score is 1 low score is not -1
 *  3. low score is -1, high score is not 1
 *
 * Otherwise, in most cases the value is computed as:
 * -exp(-2.0*outerSum) / ((H/lambda)*(exp(-lambda) - 1)
 * The last term (exp(-lambda) - 1) can be computed in two different
 * ways depending on whether lambda is small or not.
 * outerSum is a sum of the terms
 * innerSum/j, where j is denoted by iterCounter in the code.
 * The sum is truncated when the new term innersum/j i sufficiently small.
 * innerSum is a weighted sum of the probabilities of
 * of achieving a total score i in a gapless alignment,
 * which we denote by P(i,j).
 * of exactly j characters. innerSum(j) has two parts
 * Sum over i < 0  P(i,j)exp(-i * lambda) +
 * Sum over i >=0  P(i,j)
 * The terms P(i,j) are computed by dynamic programming.
 * An earlier version was flawed in that ignored the special case 1
 * and tried to replace the tail of the computation of outerSum
 * by a geometric series, but the base of the geometric series
 * was not accurately estimated in some cases.
 *
 * @param sfp object holding scoring frequency information [in]
 * @param lambda a Karlin-Altschul parameter [in]
 * @param H a Karlin-Altschul parameter [in]
 * @return K, another Karlin-Altschul parameter
 */
double BLAST_Stat::BlastKarlinLHtoK(
	double *sfp_sprob, double sfp_score_avg,   //-> input, score probability
	int sfp_obs_min, int sfp_obs_max,          //-> input, score range
	double lambda, double H,                   //-> input, calculated Lambda and H
	double sumlimit_, int iterlimit_)          //-> input, parameters
{
	/*The next array stores the probabilities of getting each possible
	  score in an alignment of fixed length; the array is shifted
	  during part of the computation, so that
	  entry 0 is for score 0.  */
	double         *alignmentScoreProbabilities = NULL;
	int            low;    /* Lowest score (must be negative) */
	int            high;   /* Highest score (must be positive) */
	int            range;  /* range of scores, computed as high - low*/
	double          K;      /* local copy of K  to return*/
	int             i;   /*loop index*/
	int             iterCounter; /*counter on iterations*/
	int            divisor; /*candidate divisor of all scores with
	                           non-zero probabilities*/
	/*highest and lowest possible alignment scores for current length*/
	int            lowAlignmentScore, highAlignmentScore;
	int            first, last; /*loop indices for dynamic program*/
	double          innerSum;
	double          oldsum, oldsum2;  /* values of innerSum on previous
	                                     iterations*/
	double          outerSum;        /* holds sum over j of (innerSum
	                                    for iteration j/j)*/
	
	double          score_avg; /*average score*/
	/*first term to use in the closed form for the case where
	  high == 1 or low == -1, but not both*/
	double          firstTermClosedForm;  /*usually store H/lambda*/
	int             iterlimit; /*upper limit on iterations*/
	double          sumlimit; /*lower limit on contributions
	                            to sum over scores*/
	
	/*array of score probabilities reindexed so that low is at index 0*/
	double         *probArrayStartLow;
	
	/*pointers used in dynamic program*/
	double         *ptrP, *ptr1, *ptr2, *ptr1e;
	double          expMinusLambda; /*e^^(-Lambda) */
	
	if (lambda <= 0. || H <= 0.) 
	{
		/* Theory dictates that H and lambda must be positive, so
		 * return -1 to indicate an error */
		return -1.;
	}
	
	/*Karlin-Altschul theory works only if the expected score
	  is negative*/
	if (sfp_score_avg >= 0.0) return -1.;
	
	low   = sfp_obs_min;
	high  = sfp_obs_max;
	range = high - low;
	
	probArrayStartLow = &sfp_sprob[low];
	/* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
	   Karlin&Altschul (1990) */
	for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) 
	{
		if (probArrayStartLow[i] != 0.0) divisor = BLAST_Gcd(divisor, i);
	}
	
	high   /= divisor;
	low    /= divisor;
	lambda *= divisor;
	range = high - low;
	
	firstTermClosedForm = H/lambda;
	expMinusLambda      = exp((double) -lambda);

	//---- closed form return ----//
	if (low == -1 && high == 1) 
	{
		K = (sfp_sprob[low*divisor] - sfp_sprob[high*divisor]) *
			(sfp_sprob[low*divisor] - sfp_sprob[high*divisor]) / sfp_sprob[low*divisor];
		return(K);
	}
	if (low == -1 || high == 1) 
	{
		if (high != 1) 
		{
			score_avg = sfp_score_avg / divisor;
			firstTermClosedForm = (score_avg * score_avg) / firstTermClosedForm;
		}
		return firstTermClosedForm * (1.0 - expMinusLambda);
	}

	//---- iterative calculation ----//
//	sumlimit  = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
//	iterlimit = BLAST_KARLIN_K_ITER_MAX;
	sumlimit = sumlimit_;
	iterlimit = iterlimit_;

	alignmentScoreProbabilities = new double[iterlimit*range + 1];
	if (alignmentScoreProbabilities == NULL) return -1.; //new failed !!

	outerSum = 0.;
	lowAlignmentScore = highAlignmentScore = 0;
	alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

	for (
		iterCounter = 0;
		((iterCounter < iterlimit) && (innerSum > sumlimit));
		outerSum += innerSum /= ++iterCounter) 
	{
		first = last = range;
		lowAlignmentScore  += low;
		highAlignmentScore += high;
		/*dynamic program to compute P(i,j)*/
		for (
			ptrP = alignmentScoreProbabilities + (highAlignmentScore-lowAlignmentScore); 
			ptrP >= alignmentScoreProbabilities;
			*ptrP-- =innerSum) 
		{
			ptr1  = ptrP - first;
			ptr1e = ptrP - last;
			ptr2  = probArrayStartLow + first;
			for (innerSum = 0.; ptr1 >= ptr1e; ) 
			{
				innerSum += *ptr1  *  *ptr2;
				ptr1--;
				ptr2++;
			}
			if (first) --first;
			if (ptrP - alignmentScoreProbabilities <= range) --last;
		}
		
		/* Horner's rule */
		innerSum = *++ptrP;
		for( i = lowAlignmentScore + 1; i < 0; i++ ) innerSum = *++ptrP + innerSum * expMinusLambda;
		innerSum *= expMinusLambda;

		for (; i <= highAlignmentScore; ++i) innerSum += *++ptrP;
		oldsum2 = oldsum;
		oldsum  = innerSum;
	}
	
	K = -exp((double)-2.0*outerSum) /
		(firstTermClosedForm*BLAST_Expm1(-(double)lambda));
	
	// delete alignmentScoreProbabilities
	if (alignmentScoreProbabilities != NULL)delete [] alignmentScoreProbabilities;
	return K;
}
