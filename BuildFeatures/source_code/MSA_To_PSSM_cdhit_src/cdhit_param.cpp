#include "cdhit_param.h"


//----------- struct WorkingParam---------//
void WorkingParam::ControlShortCoverage( int len, const Options & options )
{
	len_eff = len;
	aln_cover_flag = 0;
	if ((options.short_coverage > 0.0) || (options.min_control>0) ) { // has alignment coverage control
		aln_cover_flag = 1;
		min_aln_lenS = (int) (double(len) * options.short_coverage);
		if ( len-options.short_control > min_aln_lenS) min_aln_lenS = len-options.short_control;
		if ( options.min_control > min_aln_lenS) min_aln_lenS = options.min_control;
	}
	if (options.global_identity == 0) len_eff = min_aln_lenS; //global_identity==0
}
void WorkingParam::ControlLongCoverage( int len2, const Options & options )
{
	if (aln_cover_flag) {
		min_aln_lenL = (int) (double(len2) * options.long_coverage);
		if ( len2-options.long_control > min_aln_lenL) min_aln_lenL = len2-options.long_control;
		if ( options.min_control > min_aln_lenL) min_aln_lenL = options.min_control;
	}
}
void WorkingParam::ComputeRequiredBases( int NAA, int ss, const Options & options )
{
	// d: distance, fraction of errors;
	// e: number of errors;
	// g: length of the maximum gap;
	// m: word length;
	// n: sequence length;
	// alignment length = n - g + 1;
	// d = e / (n - g + 1);
	// e >= 1, so that, g <= n + 1 - 1/d
	// word count = (n - g - m + 1) - (e - 1)*m;
	//            = (n - g - m + 1) - (d*(n - g + 1) - 1)*m
	//            = (n - g + 1) - d*m*(n - g + 1)
	//            = (n - g + 1)*(1 - d*m)
	// minimum word count is reached when g == n + 1 - 1/d
	// so, minimum word count = 1/d - m.
	// if g == band_width: word count = (n - band + 1)*(1 - d*m);
	if( options.useDistance ){
		int invd = int( 1.0 / (options.distance_thd + 1E-9) );
		int ks = len_eff - ss + 1;
		int kn = len_eff - NAA + 1;
		int ks2 = invd - ss;
		int kn2= invd - NAA;
		//if( ks3 > ks2 ) ks2 = ks3;
		//if( kn3 > kn2 ) kn2 = kn3;
		required_aa1 = required_aas = (ks2 < ks ? ks2 : ks);
		required_aan = kn2 < kn ? kn2 : kn;
		if( required_aa1 <=0 ) required_aa1 = required_aas = 1;
		if( required_aan <=0 ) required_aan = 1;
		//required_aa1 = required_aas = required_aan = 0;
		return;
	}
#if 0
#endif
#if 1
	// (N-K)-K*(1-C)*N = C*K*N-(K-1)*N-K = (C*K-K+1)*N-K
	//required_aa1 = int((ss*aa1_cutoff-ss+1)*len_eff) - ss;
	required_aa1 = (len_eff - ss) - int(ss * ceil( (1.0 - aa1_cutoff) * len_eff ));
	if( required_aa1 < 0 ) required_aa1 = 0;
	required_aas = required_aa1;
	//required_aan = int((NAA*aa1_cutoff-NAA+1)*len_eff) - NAA;
	required_aan = (len_eff - NAA) - int(NAA * ceil( (1.0 - aa1_cutoff) * len_eff ));
	//printf( "%i %i\n", required_aa1, required_aan );
	if( required_aan < 0 ) required_aan = 0;
#if 0
	if( required_aan ){
		int part = NAA + required_aan - 1;
		int len_eff2 = len_eff - part;
		required_aa1 = (part - ss + 1) + (len_eff2 - ss) - int(ss * ceil( (1.0 - aa1_cutoff) * len_eff2 ));
		if( required_aa1 < 0 ) required_aa1 = 0;
		required_aas = required_aa1;
	}
#endif

	int aa1_old = int (aa1_cutoff* (double) len_eff) - ss + 1;
	int aas_old = int (aas_cutoff* (double) len_eff);
	int aan_old = int (aan_cutoff* (double) len_eff);

	double thd = options.cluster_thd;
	//double rest = (len_eff - ss) / double(len_eff * ss);
	double rest = (len_eff - NAA) / double(len_eff * NAA);
	double thd0 = 1.0 - rest;
	double fnew = 0;
	double fold = 1;
	if( thd > thd0 ){
		fnew = (thd - thd0) / rest;
		fold = 1.0 - fnew;
	}
	//printf( "%g %g %g\n", thd, thd0, fnew );

	required_aa1 = (int)(fnew*required_aa1 + fold*aa1_old);
	required_aas = (int)(fnew*required_aas + fold*aas_old);
	required_aan = (int)(fnew*required_aan + fold*aan_old);
#else
	required_aa1 = int (aa1_cutoff* (double) len_eff) - ss + 1;
	//printf( "%g %i\n", aa1_cutoff, len_eff );
	required_aas = (aa1_cutoff > 0.95) ?
		len_eff-ss  +1-(len_eff-required_aa1)*ss   :
		int (aas_cutoff* (double) len_eff);
	required_aan = (aa1_cutoff > 0.95) ?
		len_eff-NAA+1-(len_eff-required_aa1)*NAA :
		int (aan_cutoff* (double) len_eff);
#endif
}
