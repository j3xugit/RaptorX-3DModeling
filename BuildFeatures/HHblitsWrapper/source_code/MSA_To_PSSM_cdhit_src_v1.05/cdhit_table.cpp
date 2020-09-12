#include "cdhit_table.h"


//----------- class Sequence---------//
WordTable::WordTable( int naa, int naan )
{
	NAA      = 0;
	NAAN     = 0;
	is_aa    = 1;
	size = 0;
	frag_count = 0;
	Init( naa, naan );
}

void WordTable::SetDNA()
{
	is_aa = 0;
}

void WordTable::Init(int naa, int naan)
{
	NAA  = naa;
	NAAN = naan;
	indexCounts.resize( NAAN );
}

void WordTable::Clear()
{
	int i;
	size = 0;
	frag_count = 0;
	sequences.clear();
	for (i=0; i<NAAN; i++) indexCounts[i].Clear();
}

int WordTable::AddWordCounts( NVector<IndexCount> & counts, Sequence *seq, bool skipN)
{
	int aan_no = counts.Size();
	int i, j, k, idx = sequences.size();
	for (i=0; i<aan_no; i++) {
		IndexCount ic = counts[i];
		if ( (k=ic.count) ) {
			j = ic.index;
			if ( skipN && j<0) continue; // for those has 'N'
			NVector<IndexCount> & row = indexCounts[j];
			int cap = row.capacity;
			ic.index = idx;
			row.Append( ic );
			size += row.capacity - cap;
		}
	}
	sequences.Append( seq );
	return OK_FUNC;
}
int WordTable::AddWordCountsFrag( NVector<IndexCount> & counts, int frag, int frag_size, int repfrag )
{
	return 0;
}
int WordTable::AddWordCounts(int aan_no, Vector<int> & word_encodes, Vector<INTs> & word_encodes_no, int idx, bool skipN)
{
	int i, j, k;
	//printf( "seq %6i: ", idx );
	for (i=0; i<aan_no; i++) {
		if ( (k=word_encodes_no[i]) ) {
			j = word_encodes[i];
			if( skipN && j<0) continue; // for those has 'N'
			NVector<IndexCount> & row = indexCounts[j];
			int cap = row.capacity;
			row.Append( IndexCount( idx, k ) );
			size += row.capacity - cap;
			//if( k >1 ) printf( " %3i", k );
		}
	}
	//printf( "\n" );
	return OK_FUNC;
}

int WordTable::AddWordCountsFrag( int aan_no, Vector<int> & word_encodes,
    Vector<INTs> & word_encodes_no, int frag, int frag_size )
{
	int i, j, k, i1, k1, fra;

	for (i=0; i<frag; i++) {
		k = (i+1)*frag_size < aan_no ? (i+1)*frag_size-1: aan_no-1;
		//quick_sort(&word_encodes[0], i*frag_size, k);
		std::sort( word_encodes.begin() + i*frag_size, word_encodes.begin() + k + 1 );
	}
	for(j=aan_no-1; j; j--) {
		if (word_encodes[j] == word_encodes[j-1]) {
			word_encodes_no[j-1] += word_encodes_no[j];
			word_encodes_no[j]=0;
		}
	}
	// END check_word_encodes

	for (i=0; i<aan_no; i+=frag_size) {
		k = frag_size < (aan_no-i) ? frag_size : (aan_no -i);
		fra = i / frag_size;
		//AddWordCounts(k, word_encodes+i, word_encodes_no+i, NR90f_no+fra);
		for (i1=i; i1<i+k; i1++) {
			if ( (k1=word_encodes_no[i1]) ) {
				j = word_encodes[i1];
				NVector<IndexCount> & row = indexCounts[j];
				int cap = row.capacity;
				row.Append( IndexCount( frag_count + fra, k1 ) );
				size += row.capacity - cap;
			}
		}
	}
	frag_count += frag;

	return 0;
}

int WordTable::CountWords(int aan_no, Vector<int> & word_encodes, Vector<INTs> & word_encodes_no,
    NVector<IndexCount> &lookCounts, NVector<uint32_t> & indexMapping, 
	bool est, int min)
{
	int  j, k, j0, j1, k1;
	IndexCount tmp;

	IndexCount *ic = lookCounts.items;
	for(j=0; j<lookCounts.size; j++, ic++) indexMapping[ ic->index ] = 0;
	lookCounts.size = 0;

	int *we = & word_encodes[0];
	j0 = 0;
	if( est ) while( *we <0 ) j0++, we++; // if met short word has 'N'
	INTs *wen = & word_encodes_no[j0];
	//printf( "\nquery : " );
	for (; j0<aan_no; j0++, we++, wen++) {
		j  = *we;
		j1 = *wen;
		//if( j1 >1 ) printf( " %3i", j1 );
		if( j1==0 ) continue;
		NVector<IndexCount> & one = indexCounts[j];
		k1 = one.Size();
		IndexCount *ic = one.items;

		int rest = aan_no - j0 + 1;
		for (k=0; k<k1; k++, ic++){
			int c = ic->count < j1 ? ic->count : j1;
			uint32_t *idm = indexMapping.items + ic->index;
			if( *idm ==0 ){
				if( rest < min ) continue;
				IndexCount *ic2 = lookCounts.items + lookCounts.size;
				lookCounts.size += 1;
				*idm = lookCounts.size;
				ic2->index = ic->index;
				ic2->count = c;
			}else{
				lookCounts[ *idm - 1 ].count += c;
			}
		}
	}
	//printf( "%6i %6i\n", S, lookCounts.size );
	lookCounts[ lookCounts.size ].count = 0;
	//printf( "\n\n" );
	return OK_FUNC;
}

void WordTable::PrintAll()
{
	int  i, j, k;
	int cols = 0;
	long long total_words = 0;
	k = 0;
	for (i=0; i<NAAN; i++) {
		int size = indexCounts[i].Size();
		if ( size == 0 ) continue;
		cols++;
		cout << k << "\t" << i << "\tsize:" << size << "\t";
		for (j=0; j<size; j++) {
			cout << indexCounts[i][j].index << "," << indexCounts[i][j].count << " ";
			total_words += indexCounts[i][j].count;
		}
		cout << endl;
		k++;
	}

	cout << "total cols: " << cols << " total words: " << total_words << endl;
}

