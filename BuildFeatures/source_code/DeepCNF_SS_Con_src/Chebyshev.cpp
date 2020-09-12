#include "Chebyshev.h"

//========== constructor and destructor ============//
Chebyshev::Chebyshev(double beta)
{
	cheby_beta=beta;
}
Chebyshev::~Chebyshev(void)
{
}

//-> sigmoid function
double Chebyshev::sigmoid(double x)
{
	return (double)1 / ((double)1 + exp(-cheby_beta * x));
}

//-> Chebyshev expansion for step function
void Chebyshev::cheby_step_function(int maxdeg,double *cheby_coeff)
{
	int i;
	double _2_PI = 0.6366197723675813430755350534900574L;
	double PI_2 = 1.5707963267948966192313216916397514L;
	cheby_coeff[0] = (double)0.5;
	cheby_coeff[1] = _2_PI;
	for(i = 2; i <= maxdeg; i++)
	{
		double c;
		// i-th Chebyshev coefficient of the step function 
		if(i % 2 == 0) 
		{
			c = (double)0;
		}
		else
		{
			c = (double)1 / i;
			if(i % 4 == 3)c = -c;
			c /= PI_2;
		}
		cheby_coeff[i] = c;
	}
}

//-> Chebyshev expansion for any function [  double (*f)(double) ]
void Chebyshev::cheby_expand(int maxdeg,double *cheby_coeff)
{
	int i,j;
	double *ftab=new double[maxdeg + 1];
	double fac;
	double PI=3.1415926535897932384626433832795029L;
	//init
	for(i = 0; i <= maxdeg; i++)
	{
		cheby_coeff[i] = (double)0;
		ftab[i] = (double)0;
	}
	fac = (double)2 / ((double)maxdeg + (double)1);
	//calc
	for(i = 0; i <= maxdeg; i++)
	{
		double y;
		y = cos(PI * ((double)i + (double)0.5) / (double)(maxdeg + 1));
		ftab[i] = sigmoid(y);
	}
	for(i = 0; i <= maxdeg; i++)
	{
		double sum = (double)0;
		for(j = 0; j <= maxdeg; j++)
		{
			sum += ftab[j] * cos(PI * (double)i * ((double)j + (double)0.5) / (double)(maxdeg + 1));
		}
		cheby_coeff[i] = fac * sum;
	}
	//delete
	delete [] ftab;
	cheby_coeff[0] /= 2;
}

//-> convert Chebyshev expansion into a standard polynomial. 
int Chebyshev::cheby_convert(int maxdeg,double *cheby_exp_coeff,double *gamma_coeff)
{
	int i,j;
	double *prev_coeff=new double[maxdeg + 1];
	double *cheby_coeff1=new double[maxdeg + 1];
	double *cheby_coeff2=new double[maxdeg + 1];
	//init
	for(i = 0; i <= maxdeg; i++)
	{
		prev_coeff[i] = (double)0;
		cheby_coeff1[i] = (double)0;
		cheby_coeff2[i] = (double)0;
	}
	//calc
	cheby_coeff1[1] = (double)1;
	cheby_coeff2[0] = (double)1;
	prev_coeff[0] = cheby_exp_coeff[0];
	prev_coeff[1] = cheby_exp_coeff[1];
	for(i = 2; i <= maxdeg; i++)
	{
		double *swap_tmp;
		//-> compute coefficients of the next Chebyshev polynomial 
		cheby_coeff2[0] = -cheby_coeff2[0];
		for(j = 1; j <= i; j++)
		{
			cheby_coeff2[j] = 2 * cheby_coeff1[j-1] - cheby_coeff2[j];
		}
		swap_tmp = cheby_coeff2;
		cheby_coeff2 = cheby_coeff1;
		cheby_coeff1 = swap_tmp;
		//-> update the polynomial 
		for(j = 0; j <= i; j++)
		{
			prev_coeff[j] += cheby_coeff1[j] * cheby_exp_coeff[i];
		}
	}

	//------ precompute gamma coefficients ------//
	int idx=0;
	for(i = 0; i <= maxdeg; i++)
	{
		double ck = prev_coeff[i];
		long nck = (long)1;
		for(j = 0; j <= i; j++)
		{
			double alpha_kl = ck * (double)nck;
			if((i - j) % 2 != 0)alpha_kl = -alpha_kl;
			gamma_coeff[idx++] = alpha_kl;
			nck = nck * (i - j) / (j + 1);
		}
	}
	//delete
	delete [] prev_coeff;
	delete [] cheby_coeff2;
	delete [] cheby_coeff1;
	//return
	return idx;
}
