#pragma once
#include "cdhit_util.h"
#include "cdhit_seq.h"


//====== class: WordTable ======//start
class WordTable
{
	private:
	public:
		Vector<NVector<IndexCount> > indexCounts; // hold index and word counts of seqs
		Vector<Sequence*>            sequences;
		int     NAA;                // length of word
		int     NAAN;               // rows of table
		char    is_aa;              // aa is for prot
		size_t  size;
		int     frag_count;

	public:
		WordTable( int naa=0, int naan=0 );
		void Init(int, int);
		void Clear();
		void SetDNA();
		int  AddWordCounts( NVector<IndexCount> & counts, Sequence *seq, bool skipN=false);
		int  AddWordCountsFrag( NVector<IndexCount> & counts, int frag, int frag_size, int repfrag );

		int  AddWordCounts(int aan_no, Vector<int> & word_encodes, 
				Vector<INTs> & word_encodes_no, int idx, bool skipN=false);
		int AddWordCountsFrag( int aan_no, Vector<int> & word_encodes, 
				Vector<INTs> & word_encodes_no, int frag, int frag_size );
		int CountWords(int aan_no, Vector<int> & aan_list, Vector<INTs> & aan_list_no, 
				NVector<IndexCount> & lookCounts, NVector<uint32_t> & indexMapping,
				bool est=false, int min=0);
		void PrintAll();
}; // END class INDEX_TBL
//====== class: WordTable ======//end

