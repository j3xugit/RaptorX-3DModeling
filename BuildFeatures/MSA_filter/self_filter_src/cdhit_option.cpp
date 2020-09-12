#include "cdhit_option.h"


bool Options::SetOptionCommon( const char *flag, const char *value )
{
	int intval = atoi( value );
	if      (strcmp(flag, "-i" ) == 0) input = value;
	else if (strcmp(flag, "-o" ) == 0) output = value;
	else if (strcmp(flag, "-M" ) == 0) max_memory  = atol(value) * 1000000;
	else if (strcmp(flag, "-l" ) == 0) min_length  = intval;
	else if (strcmp(flag, "-c" ) == 0) cluster_thd  = atof(value), useIdentity = true;
	else if (strcmp(flag, "-D" ) == 0) distance_thd  = atof(value), useDistance = true;
	else if (strcmp(flag, "-b" ) == 0) band_width  = intval;
	else if (strcmp(flag, "-n" ) == 0) NAA       = intval;
	else if (strcmp(flag, "-d" ) == 0) des_len   = intval;
	else if (strcmp(flag, "-s" ) == 0) diff_cutoff  = atof(value);
	else if (strcmp(flag, "-S" ) == 0) diff_cutoff_aa  = intval;
	else if (strcmp(flag, "-B" ) == 0) store_disk  = intval;
	else if (strcmp(flag, "-p" ) == 0) print  = intval;
	else if (strcmp(flag, "-g" ) == 0) cluster_best  = intval;
	else if (strcmp(flag, "-G" ) == 0) global_identity  = intval;
	else if (strcmp(flag, "-aL") == 0) long_coverage = atof(value);
	else if (strcmp(flag, "-AL") == 0) long_control = intval;
	else if (strcmp(flag, "-aS") == 0) short_coverage = atof(value);
	else if (strcmp(flag, "-AS") == 0) short_control = intval;
	else if (strcmp(flag, "-A" ) == 0) min_control  = intval;
	else if (strcmp(flag, "-uL") == 0) long_unmatch_per = atof(value);
	else if (strcmp(flag, "-uS") == 0) short_unmatch_per = atof(value);
	else if (strcmp(flag, "-U") == 0) unmatch_len = intval;
	else if (strcmp(flag, "-tmp" ) == 0) CDHIT_temp_dir  = value;
	else if (strcmp(flag, "-T" ) == 0){
#ifndef NO_OPENMP
		int cpu = omp_get_num_procs();
		threads  = intval;
		if( threads > cpu ){
			threads = cpu;
			printf( "Warning: total number of CPUs in the system is %i\n", cpu );
		}
		if( threads == 0 ){
			threads = cpu;
			printf( "total number of CPUs in the system is %i\n", cpu );
		}
		if( threads != intval ) printf( "Actual number of CPUs to be used: %i\n\n", threads );
#else
		printf( "Option -T is ignored: multi-threading with OpenMP is NOT enabled!\n" );
#endif
	}else return false;
	return true;
}
bool Options::SetOption( const char *flag, const char *value )
{
	if( is454 ){
		if( strcmp(flag, "-s") == 0 ) return false;
		else if( strcmp(flag, "-S") == 0 ) return false;
		else if( strcmp(flag, "-G") == 0 ) return false;
		else if( strcmp(flag, "-A") == 0 ) return false;
		else if( strcmp(flag, "-r") == 0 ) return false;
		else if( strcmp(flag, "-D") == 0 ){ max_indel = atoi(value); return true; }
	}
	if( SetOptionCommon( flag, value ) ) return true;
	if (strcmp(flag, "-t" ) == 0) tolerance = atoi(value);
	else if (strcmp(flag, "-F" ) == 0) frag_size = atoi(value);
	else if (has2D && SetOption2D( flag, value ) ) return true;
	else if (isEST && SetOptionEST( flag, value ) ) return true;
	else return false;
	return true;
}
bool Options::SetOption2D( const char *flag, const char *value )
{
	if( SetOptionCommon( flag, value ) ) return true;
	if (strcmp(flag, "-i2" ) == 0) input2 = value;
	else if (strcmp(flag, "-s2") == 0) diff_cutoff2 = atof(value);
	else if (strcmp(flag, "-S2") == 0) diff_cutoff_aa2 = atoi(value);
	else return false;
	return true;
}
bool Options::SetOptionEST( const char *flag, const char *value )
{
	NAA_top_limit = 12;
	if( SetOptionCommon( flag, value ) ) return true;
	if (strcmp(flag, "-r" ) == 0) option_r  = atoi(value); 
	else if (strcmp(flag, "-gap") == 0) mat.gap = MAX_SEQ * atoi(value);
	else if (strcmp(flag, "-gap-ext") == 0) mat.ext_gap = MAX_SEQ * atoi(value);
	else if (strcmp(flag, "-match") == 0) mat.set_match( atoi(value) );
	else if (strcmp(flag, "-mismatch") == 0) mat.set_mismatch( atoi(value) );
	else return false;
	return true;
}
bool Options::SetOptions( int argc, char *argv[], bool twod, bool est )
{
	printf( "================================================================\n" );
	printf( "Program: CD-HIT, V" CDHIT_VERSION WITH_OPENMP "\n");
	printf( "Command:" );
	int i;
	int n = 9;
	for(i=0; i<argc; i++){
		n += strlen( argv[i] ) + 1;
		if( n >= 64 ){
			printf( "\n         %s", argv[i] );
			n = strlen( argv[i] ) + 9;
		}else{
			printf( " %s", argv[i] );
		}
	}
	printf( "\n\n" );
	printf( "Started: \n");
	printf( "================================================================\n" );
	printf( "                            Output                              \n" );
	printf( "----------------------------------------------------------------\n" );
	has2D = twod;
	isEST = est;
	for (i=1; i+1<argc; i+=2) if ( SetOption( argv[i], argv[i+1] ) == 0) return false;
	if( i < argc ) return false;

//	atexit( CleanUpTempFiles );  //-> atexit means the program will execute CleanUpTempFiles when exit.
	return true;
}
void Options::Validate()
{
	if( useIdentity && useDistance ) CDHIT_bomb_error( "can not use both identity cutoff and distance cutoff" );
	if( useDistance ){
		if ((distance_thd > 1.0) || (distance_thd < 0.0)) CDHIT_bomb_error("invalid distance threshold");
	}else if( isEST ){
		if ((cluster_thd > 1.0) || (cluster_thd < 0.8)) CDHIT_bomb_error("invalid clstr threshold, should >=0.8");
	}else{
		if ((cluster_thd > 1.0) || (cluster_thd < 0.4)) CDHIT_bomb_error("invalid clstr");
	}

	if (band_width < 1 ) CDHIT_bomb_error("invalid band width");
	if (NAA < 2 || NAA > NAA_top_limit) CDHIT_bomb_error("invalid word length");
	if (des_len < 0 ) CDHIT_bomb_error("too short description, not enough to identify sequences");
	if (! isEST && (tolerance < 0 || tolerance > 5) ) CDHIT_bomb_error("invalid tolerance");
	if ((diff_cutoff<0) || (diff_cutoff>1)) CDHIT_bomb_error("invalid value for -s");
	if (diff_cutoff_aa<0) CDHIT_bomb_error("invalid value for -S");
	if( has2D ){
		if ((diff_cutoff2<0) || (diff_cutoff2>1)) CDHIT_bomb_error("invalid value for -s2");
		if (diff_cutoff_aa2<0) CDHIT_bomb_error("invalid value for -S2");
	}
	if (global_identity == 0) print = 1;
	if (short_coverage < long_coverage) short_coverage = long_coverage;
	if (short_control > long_control) short_control = long_control;
	if ((global_identity == 0) && (short_coverage == 0.0) && (min_control == 0))
		CDHIT_bomb_error("You are using local identity, but no -aS -aL -A option");
	if (frag_size < 0) CDHIT_bomb_error("invalid fragment size");

#if 0
	if( useDistance ){
		/* when required_aan becomes zero */
		if( distance_thd * NAA >= 1 )
			CDHIT_bomb_warning( "word length is too long for the distance cutoff" );
	}else{
		/* when required_aan becomes zero */
		if( cluster_thd <= 1.0 - 1.0 / NAA )
			CDHIT_bomb_warning( "word length is too long for the identity cutoff" );
	}
#endif

	const char *message = "Your word length is %i, using %i may be faster!\n";
	if ( ! isEST && tolerance ) {
		int i, clstr_idx = (int) (cluster_thd * 100) - naa_stat_start_percent;
		int tcutoff = naa_stat[tolerance-1][clstr_idx][5-NAA];

		if (tcutoff < 5 )
			CDHIT_bomb_error("Too low cluster threshold for the word length.\n"
					"Increase the threshold or the tolerance, or decrease the word length.");
		for ( i=5; i>NAA; i--) {
			if ( naa_stat[tolerance-1][clstr_idx][5-i] > 10 ) {
				printf( message, NAA, i );
				break;
			}
		}
	} else if( isEST ) {
		if      ( cluster_thd > 0.9  && NAA < 8 ) printf( message, NAA, 8 );
		else if ( cluster_thd > 0.87 && NAA < 5 ) printf( message, NAA, 5 );
		else if ( cluster_thd > 0.80 && NAA < 4 ) printf( message, NAA, 4 );
		else if ( cluster_thd > 0.75 && NAA < 3 ) printf( message, NAA, 3 );
	} else {
		if      ( cluster_thd > 0.85 && NAA < 5 ) printf( message, NAA, 5 );
		else if ( cluster_thd > 0.80 && NAA < 4 ) printf( message, NAA, 4 );
		else if ( cluster_thd > 0.75 && NAA < 3 ) printf( message, NAA, 3 );
	}

	if ( (min_length + 1) < NAA ) CDHIT_bomb_error("Too short -l, redefine it");
}
void Options::Print()
{
	printf( "isEST = %i\n", isEST );
	printf( "has2D = %i\n", has2D );
	printf( "NAA = %i\n", NAA );
	printf( "NAA_top_limit = %i\n", NAA_top_limit );
	printf( "min_length = %i\n", min_length );
	printf( "cluster_best = %i\n", cluster_best );
	printf( "global_identity = %i\n", global_identity );
	printf( "cluster_thd = %g\n", cluster_thd );
	printf( "diff_cutoff = %g\n", diff_cutoff );
	printf( "diff_cutoff_aa = %i\n", diff_cutoff_aa );
	printf( "tolerance = %i\n", tolerance );
	printf( "long_coverage = %g\n", long_coverage );
	printf( "long_control = %i\n", long_control );
	printf( "short_coverage = %g\n", short_coverage );
	printf( "short_control = %i\n", short_control );
	printf( "frag_size = %i\n", frag_size );
	printf( "option_r = %i\n", option_r );
	printf( "print = %i\n", print );
}



