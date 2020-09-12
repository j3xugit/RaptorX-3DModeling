#pragma once
#include "cdhit_util.h"
#include "cdhit_option.h"
#include "cdhit_seq.h"


//====== struct: WorkingBuffer ======//start
struct WorkingBuffer
{
	Vector<int>  taap;
	Vector<int>  word_encodes;
	Vector<int>  word_encodes_backup;
	Vector<INTs> word_encodes_no;
	Vector<INTs> aap_list;
	Vector<INTs> aap_begin;
	//Vector<IndexCount>  indexCounts;
	NVector<IndexCount>  lookCounts;
	NVector<uint32_t>    indexMapping;
	MatrixInt64  score_mat;
	MatrixInt    back_mat;
	Vector<int>  diag_score;
	Vector<int>  diag_score2;
	Vector<int> aan_list_comp;
	char seqi_comp[MAX_SEQ];
	int total_bytes;

	WorkingBuffer( int frag=0, const Options & options=Options() ){
		Set( frag, options );
	}
	void Set( int frag, const Options & options ){
		bool est = options.isEST;
		int m = MAX_UAA*MAX_UAA;
		int max_len = options.max_length;
		int band = max_len*max_len;
		if( est ) m = m * m;
		if( band > options.band_width ) band = options.band_width;
		taap.resize( m );
		aap_list.resize( max_len );
		aap_begin.resize( m );
		//indexCounts.resize( max_len );
		word_encodes.resize( max_len );
		word_encodes_no.resize( max_len );
		word_encodes_backup.resize( max_len );
		/* each table can not contain more than MAX_TABLE_SEQ representatives or fragments! */
		if( frag > MAX_TABLE_SEQ ) frag = MAX_TABLE_SEQ;
		lookCounts.Resize( frag + 2 );
		indexMapping.Resize( frag + 2 );
		diag_score.resize( MAX_DIAG );
		diag_score2.resize( MAX_DIAG );
		aan_list_comp.resize( max_len );
		total_bytes = max_len;
		total_bytes += taap.size()*sizeof(int);
		total_bytes += word_encodes.size()*sizeof(int);
		total_bytes += word_encodes_backup.size()*sizeof(int);
		total_bytes += diag_score.size()*sizeof(int);
		total_bytes += diag_score2.size()*sizeof(int);
		total_bytes += aan_list_comp.size()*sizeof(int);
		total_bytes += word_encodes_no.size()*sizeof(INTs);
		total_bytes += aap_list.size()*sizeof(INTs);
		total_bytes += aap_begin.size()*sizeof(INTs);
		total_bytes += indexMapping.Size()*sizeof(uint32_t);
		//total_bytes += indexCounts.size()*sizeof(IndexCount);
		total_bytes += lookCounts.Size()*sizeof(IndexCount);
		total_bytes += max_len*(band*sizeof(int)+sizeof(VectorInt));
		total_bytes += max_len*(band*sizeof(int)+sizeof(VectorInt64));
	}

	int EncodeWords( Sequence *seq, int NA, bool est = false );
	void ComputeAAP( const char *seqi, int size );
	void ComputeAAP2( const char *seqi, int size );
};

