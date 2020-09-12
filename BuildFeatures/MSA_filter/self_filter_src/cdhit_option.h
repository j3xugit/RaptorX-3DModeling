#pragma once
#include "cdhit_util.h"
#include "cdhit_scomat.h"

//-> version define
#ifndef CDHIT_VERSION
#define CDHIT_VERSION  "4.5.4"
#endif

//-> OpenMP define
#ifndef NO_OPENMP
	#include <omp.h>
	#define WITH_OPENMP " (+OpenMP)"
#else
	#define WITH_OPENMP ""
	#define omp_set_num_threads(T) (T = T)
	#define omp_get_thread_num() 0
#endif


//====== struct: Options ======//start
struct Options
{
	int     NAA;
	int     NAAN;
	int     NAA_top_limit;

	size_t  max_memory; // -M: 400,000,000 in bytes
	int     min_length; // -l: 10 bases
	bool    cluster_best;  // -g: 0, the first; 1, the best
	bool    global_identity; // -G:
	bool    store_disk; // -B:
	int     band_width; // -b: 20
	double  cluster_thd; // -c
	double  distance_thd; // -D
	double  diff_cutoff; // -s: 0.0
	double  diff_cutoff2; // -s2: 1.0
	int     diff_cutoff_aa; // -S: 999999
	int     diff_cutoff_aa2; // -S2: 0
	int     tolerance; // -t: 2
	double  long_coverage; // -aL:
	int     long_control; // -AL:
	double  short_coverage; // -aS:
	int     short_control; // -AS:
	int     min_control; // -A:
	double  long_unmatch_per; // -uL
	double  short_unmatch_per; // -uS
	int     unmatch_len; // -U
	int     max_indel; // -D
	int     print;
	int     des_len;
	int     frag_size;
	int     option_r;
	int     threads;
	int     verbose;

	int     max_length;
	unsigned long long    max_entries;
	size_t  mem_limit;
	size_t  total_letters;
	size_t  total_desc;

	bool    has2D;
	bool    isEST;
	bool    is454;
	bool    useIdentity;
	bool    useDistance;

	string  input;
	string  input2;
	string  output;

	Options(){
		useIdentity = false;
		useDistance = false;
		has2D = false;
		isEST = false;
		is454 = false;
		NAA = 5;
		NAA_top_limit = 5;
		cluster_thd = 0.9;
		distance_thd = 0.0;
		max_memory = 800000000;
		min_length = 10;
		cluster_best = false;
		global_identity = true;
		store_disk = false;
		band_width = 20;
		diff_cutoff = 0.0;
		diff_cutoff2 = 1.0;
		diff_cutoff_aa = 999999;
		diff_cutoff_aa2 = 0;
		tolerance = 2;
		long_coverage = 0.0;
		long_control = 99999999;
		short_coverage = 0.0;
		short_control = 99999999;
		long_unmatch_per = 1.0;
		short_unmatch_per = 1.0;
		unmatch_len = 99999999;
		min_control = 0;
		max_indel = 1;
		print = 0;
		option_r  = 1;
		frag_size = 0;
		des_len = 20;
		threads = 1;
		verbose = 0;
		max_length = 0;
		max_entries = 0;
		mem_limit = 100000000;
		total_letters = 0;
		total_desc = 0;
	};

	//-> ScoreMatrix
	ScoreMatrix mat;

	//-> set option functions
	bool SetOptionCommon( const char *flag, const char *value );
	bool SetOption( const char *flag, const char *value );
	bool SetOption2D( const char *flag, const char *value );
	bool SetOptionEST( const char *flag, const char *value );
	bool SetOptions( int argc, char *argv[], bool twodata=false, bool est=false );

	//-> validate and print
	void Validate();
	void Print();
};
//====== struct: Options ======//end
