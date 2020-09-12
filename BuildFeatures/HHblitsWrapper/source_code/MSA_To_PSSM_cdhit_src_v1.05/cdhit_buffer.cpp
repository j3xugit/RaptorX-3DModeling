#include "cdhit_buffer.h"


//----------- struct WorkingBuffer---------//
int WorkingBuffer::EncodeWords( Sequence *seq, int NAA, bool est )
{
	char *seqi = seq->data;
	int len = seq->size;
	// check_word_encodes
	int aan_no = len - NAA + 1;
	int i, j, i0, i1;
	int skip = 0;
	unsigned char k, k1;
	for (j=0; j<aan_no; j++) {
		char *word = seqi + j;
		int encode = 0;
		for (k=0, k1=NAA-1; k<NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
		word_encodes[j] = word_encodes_backup[j] = encode;
	}

	if( est ){
		for (j=0; j<len; j++){
			if ( seqi[j] == 4 ) {                      // here N is 4
				i0 = (j-NAA+1 > 0) ? j-NAA+1 : 0;
				i1 = j < aan_no ? j : aan_no - 1;
				for (i=i0; i<=i1; i++) word_encodes[i]=-1;
			}
		}
		for (j=0; j<aan_no; j++) skip += (word_encodes[j] == -1);
	}

	std::sort( word_encodes.begin(), word_encodes.begin() + aan_no );
	for(j=0; j<aan_no; j++) word_encodes_no[j]=1;
	for(j=aan_no-1; j; j--) {
		if (word_encodes[j] == word_encodes[j-1]) {
			word_encodes_no[j-1] += word_encodes_no[j];
			word_encodes_no[j]=0;
		}
	}
	return skip;
	// END check_word_encodes
}
void WorkingBuffer::ComputeAAP( const char *seqi, int size )
{
	int len1 = size - 1;
	int sk, j1, mm, c22;
	for (sk=0; sk<NAA2; sk++) taap[sk] = 0;
	for (j1=0; j1<len1; j1++) {
		c22= seqi[j1]*NAA1 + seqi[j1+1];
		taap[c22]++;
	}
	for (sk=0,mm=0; sk<NAA2; sk++) {
		aap_begin[sk] = mm; mm+=taap[sk]; taap[sk] = 0;
	}
	for (j1=0; j1<len1; j1++) {
		c22= seqi[j1]*NAA1 + seqi[j1+1];
		aap_list[aap_begin[c22]+taap[c22]++] =j1;
	}
}
void WorkingBuffer::ComputeAAP2( const char *seqi, int size )
{
	int len1 = size - 3;
	int sk, j1, mm, c22;
	for (sk=0; sk<NAA4; sk++) taap[sk] = 0;
	for (j1=0; j1<len1; j1++) {
		if ((seqi[j1]==4) || (seqi[j1+1]==4) || (seqi[j1+2]==4) || (seqi[j1+3]==4)) continue; //skip N
		c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
		taap[c22]++;
	}
	for (sk=0,mm=0; sk<NAA4; sk++) {
		aap_begin[sk] = mm;  mm += taap[sk];  taap[sk] = 0;
	}
	for (j1=0; j1<len1; j1++) {
		if ((seqi[j1]==4) || (seqi[j1+1]==4) || (seqi[j1+2]==4) || (seqi[j1+3]==4)) continue; //skip N
		c22 = seqi[j1]*NAA3 + seqi[j1+1]*NAA2 + seqi[j1+2]*NAA1 + seqi[j1+3];
		aap_list[aap_begin[c22]+taap[c22]++] =j1;
	}
}
