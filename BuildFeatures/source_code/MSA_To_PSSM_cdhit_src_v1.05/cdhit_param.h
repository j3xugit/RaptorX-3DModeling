#pragma once
#include "cdhit_util.h"
#include "cdhit_option.h"


//====== struct: WorkingParam ======//start
struct WorkingParam
{
	double aa1_cutoff;
	double aas_cutoff; /* or aa2 */
	double aan_cutoff;
	int    len_upper_bound;
	int    len_lower_bound;

	WorkingParam( double a1=0, double a2=0, double an=0 ){
		Set( a1, a2, an );
	}
	void Set( double a1=0, double a2=0, double an=0 ){
		aa1_cutoff = a1;
		aas_cutoff = a2;
		aan_cutoff = an;
		len_upper_bound = 0;
		len_lower_bound = 0;
	}

	int len_eff;
	int aln_cover_flag;
	int min_aln_lenS;
	int min_aln_lenL;
	int required_aa1;
	int required_aas; /* or aa2 */
	int required_aan;

	void ControlShortCoverage( int len, const Options & option );
	void ControlLongCoverage( int len, const Options & option );
	void ComputeRequiredBases( int NAA, int ss, const Options & option );
};
//====== struct: WorkingParam ======//end

