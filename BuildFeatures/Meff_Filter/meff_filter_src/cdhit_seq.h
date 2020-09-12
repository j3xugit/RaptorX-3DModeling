#pragma once
#include "cdhit_util.h"
#include "cdhit_option.h"


//====== struct: Sequence ======//start
struct Sequence
{
	// real sequence, if it is not stored swap file:
	char *data;
	// length of the sequence:
	int   size;
	int   bufsize;

	//uint32_t stats;

	// if swap != NULL, the sequence is stored in file.
	// swap is opened as temporary file, which will be deleted automatically
	// after the program is finished:
	FILE *swap;
	// stream offset of the sequence:
	int   offset;

	// stream offset of the description string in the database:
	size_t   des_begin;
	// length of the description:
	int   des_length;
	// length of the description in quality score part:
	int   des_length2;
	// length of data in fasta file, including line wrapping:
	int   dat_length;

	char *identifier;

	// index of the sequence in the original database:
	int   index;
	short state;
	int   cluster_id;
	float identity;
	float distance;
	int   coverage[4];

	Sequence();
	Sequence( const Sequence & other );
	~Sequence();

	void Clear();

	void operator=( const char *s );
	void operator+=( const char *s );

	void Resize( int n );
	void Reserve( int n );

	void Swap( Sequence & other );
	void Format();

	void ConvertBases();

	void SwapIn();
	void SwapOut();
	void PrintInfo( int id, FILE *fin, FILE *fout, const Options & options, char *buf );
};
//====== struct: Sequence ======//end

