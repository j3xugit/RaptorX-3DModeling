#pragma once
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

//====== Chebyshev class ======//
class Chebyshev
{
public:
	Chebyshev(double beta=-1);
	~Chebyshev();

//---- variables ----//
public:
	double cheby_beta;            //-> parameter beta for sigmoid function

//---- functions ----//
public:
	double sigmoid(double x);
	void cheby_step_function(int maxdeg,double *cheby_coeff);
	void cheby_expand(int maxdeg,double *cheby_coeff);
	int cheby_convert(int maxdeg,double *cheby_exp_coeff,double *gamma_coeff);
};
