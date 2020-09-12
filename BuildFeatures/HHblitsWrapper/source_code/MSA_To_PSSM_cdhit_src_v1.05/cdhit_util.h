#pragma once
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <string>
#include <valarray>
#include <vector>
#include <map>
using namespace std;


//--------- basic definition ---------//
#define MAX_AA 23
#define MAX_NA 6
#define MAX_UAA 21
#define MAX_SEQ 655360
#define MAX_DIAG (MAX_SEQ<<1)              // MAX_DIAG be twice of MAX_SEQ
#define MAX_GAP MAX_SEQ                    // MAX_GAP <= MAX_SEQ
#define MAX_DES 300000
#define MAX_LINE_SIZE 300000
#define MAX_FILE_NAME 1280
#define MAX_SEG 50
#define MAX_BIN_SWAP 2E9
#define MAX_TABLE_SIZE 200000000
#define CLOCK_TICKS 100
#define FAILED_FUNC 1
#define OK_FUNC 0

#define IS_REP 1
#define IS_REDUNDANT 2
#define IS_PROCESSED 16
#define IS_MINUS_STRAND 32

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

typedef unsigned int UINT4;
typedef unsigned short UINT2;


//if the longset sequence is longer than 65535, I use INT4
#define LONG_SEQ
#ifdef LONG_SEQ
#define INTs UINT4
#else
#define INTs UINT2
#endif

//define the maximal table sequence
#define MAX_TABLE_SEQ (1<<22)

//--- for dynamic programming traceback ----//
enum { DP_BACK_NONE=0, DP_BACK_LEFT_TOP=1, DP_BACK_LEFT=2, DP_BACK_TOP=3 };



//====== struct: IndexCount ======//start
struct IndexCount
{
	int index;
	int count;
	IndexCount( int i=0, int c=0 ){ index = i, count = c; }
};
//====== struct: IndexCount ======//end


// the parent containter must guarantee continuous memory allocation.
// std::valarray could be used instead of std::vector.
template<class TYPE>
class Vector : public vector<TYPE>
{
	public:
		Vector() : vector<TYPE>(){}
		Vector( size_t size ) : vector<TYPE>( size ){}
		Vector( size_t size, const TYPE & deft ) : vector<TYPE>( size, deft ){}

		void Append( const TYPE & item ){
			size_t n = this->size();
			if( n + 1 >= this->capacity() ) this->reserve( n + n/5 + 1 );
			this->push_back( item );
		}
};

// for primitive types only
template<class TYPE>
class NVector
{
	public:
		TYPE   *items;
		int     size;
		int     capacity;

		NVector(){ size = capacity = 0; items = NULL; }
		NVector( int n, const TYPE & v=TYPE() ){ 
			size = capacity = 0; items = NULL; 
			Resize( n, v );
		}
		NVector( const NVector & other ){
			size = capacity = 0; items = NULL; 
			if( other.items ){
				Resize( other.size );
				memcpy( items, other.items, other.size * sizeof(TYPE) );
			}
		}

		~NVector(){ if( items ) free( items ); }

		int  Size()const{ return size; }
		void Clear(){
			if( items ) free( items );
			size = capacity = 0; items = NULL; 
		}

		void Resize( int n, const TYPE & value=TYPE() ){
			if( n == size && capacity > 0 ) return;
			int i;
			// When resize() is called, probably this is the intended size,
			// and will not be changed frequently.
			if( n != capacity ){
				capacity = n;
				items = (TYPE*)realloc( items, capacity*sizeof(TYPE) );
			}
			for(i=size; i<n; i++ ) items[i] = value;
			size = n;
		}
		void Append( const TYPE & item ){
			if( size + 1 >= capacity ){
				capacity = size + size/5 + 1;
				items = (TYPE*)realloc( items, capacity*sizeof(TYPE) );
			}
			items[size] = item;
			size ++;
		}

		TYPE& operator[]( const int i ){
			//if( i <0 or i >= size ) printf( "out of range\n" );
			return items[i];
		}
		TYPE& operator[]( const int i )const{
			//if( i <0 or i >= size ) printf( "out of range\n" );
			return items[i];
		}
};

//----- typedef definition --------//
typedef NVector<int>      VectorInt;
typedef Vector<VectorInt> MatrixInt;
typedef NVector<int64_t>   VectorInt64;
typedef Vector<VectorInt64> MatrixInt64;
typedef NVector<INTs>      VectorIntX;
typedef Vector<VectorIntX> MatrixIntX;



//=========== temporary files related =========//
extern const char *CDHIT_temp_dir;
extern FILE* CDHIT_OpenTempFile( const char *dir = NULL );
extern void CDHIT_CleanUpTempFiles();


//=========== AA matrix related =============//
extern const char CDHIT_aa[];
extern int CDHIT_aa2idx[];
extern int CDHIT_BLOSUM62[];
extern int CDHIT_na2idx[];
extern int CDHIT_BLOSUM62_na[];
extern void CDHIT_setaa_to_na();


//=========== bomb error and warning ===========//
extern void CDHIT_bomb_error(const char *message);
extern void CDHIT_bomb_error(const char *message, const char *message2);
extern void CDHIT_bomb_warning(const char *message);
extern void CDHIT_bomb_warning(const char *message, const char *message2);


//------- format sequence ----//
extern void CDHIT_format_seq(char *seq);

//------- quick sort -----//
extern void PartialQuickSort( IndexCount *data, int first, int last, int partial );


//============ data structure =========//
//------- AA array data structure ---------//
extern void InitNAA( int max );
extern int NAA1 ;
extern int NAA2 ;
extern int NAA3 ;
extern int NAA4 ;
extern int NAA5 ;
extern int NAA6 ;
extern int NAA7 ;
extern int NAA8 ;
extern int NAA9 ;
extern int NAA10;
extern int NAA11;
extern int NAA12;
extern int NAAN_array[13];

//----- statistics data ---//
extern int naa_stat_start_percent;
extern int naa_stat[5][61][4];

