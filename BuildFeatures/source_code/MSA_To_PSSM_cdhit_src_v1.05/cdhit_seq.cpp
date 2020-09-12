#include "cdhit_seq.h"


//----------- struct Sequence---------//
Sequence::Sequence()
{
	memset( this, 0, sizeof( Sequence ) );
	distance = 2.0;
}
Sequence::Sequence( const Sequence & other )
{
	//printf( "new: %p  %p\n", this, & other );
	memcpy( this, & other, sizeof( Sequence ) );
	distance = 2.0;
	if( other.data ){
		size = bufsize = other.size;
		data = new char[size+1];
		//printf( "data: %p  %p\n", data, other.data );
		data[size] = 0;
		memcpy( data, other.data, size );
		//for (i=0; i<size; i++) data[i] = other.data[i];
	}
	if( other.identifier ){
		int len = strlen( other.identifier );
		identifier = new char[len+1];
		memcpy( identifier, other.identifier, len );
		identifier[len] = 0;
	}
}
Sequence::~Sequence()
{
	//printf( "delete: %p\n", this );
	if( data ) delete[] data;
	if( identifier ) delete[] identifier;
}

void Sequence::Clear()
{
	if( data ) delete[] data;
	/* do not set size to zero here, it is need for writing output */
	bufsize = 0;
	data = NULL;
}

void Sequence::operator=( const char *s )
{
	size = 0; // avoid copying;
	Resize( strlen( s ) );
	strcpy( data, s );
}
void Sequence::operator+=( const char *s )
{
	int m = size, n = strlen( s );
	Reserve( m + n );
	memcpy( data+m, s, n );
}
void Sequence::Resize( int n )
{
	int m = size < n ? size : n;
	size = n;
	if( size != bufsize ){
		char *old = data;
		bufsize = size;
		data = new char[ bufsize + 1 ];
		if ( data == NULL ) CDHIT_bomb_error( "Memory" );
		if ( old ){
			memcpy( data, old, m );
			delete []old;
		}
		if( size ) data[size] = 0;
	}
}
void Sequence::Reserve( int n )
{
	int m = size < n ? size : n;
	size = n;
	if( size > bufsize ){
		char *old = data;
		bufsize = size + size/5 + 1;
		data = new char[ bufsize + 1 ];
		if ( data == NULL ) CDHIT_bomb_error( "Memory" );
		if ( old ){
			memcpy( data, old, m );
			delete []old;
		}
		if( size ) data[size] = 0;
	}
}
void Sequence::ConvertBases()
{
	int i;
	for(i=0; i<size; i++) data[i] = CDHIT_aa2idx[data[i] - 'A'];
}

void Sequence::Swap( Sequence & other )
{
	Sequence tmp;
	memcpy( & tmp, this, sizeof( Sequence ) );
	memcpy( this, & other, sizeof( Sequence ) );
	memcpy( & other, & tmp, sizeof( Sequence ) );
	memset( & tmp, 0, sizeof( Sequence ) );
}
void Sequence::Format()
{
	int i, j=0;
	for (i=0; i<size; i++){
		char ch = data[i];
		if ( isalpha( ch ) ) data[j++] = toupper( ch );
	}
	data[j] = 0;
	size = j;
}

void Sequence::SwapIn()
{
	if( data ) return;
	if( swap == NULL ) CDHIT_bomb_error( "Can not swap in sequence" );
	Resize( size );
	fseek( swap, offset, SEEK_SET );
	if( fread( data, 1, size, swap ) ==0 ) CDHIT_bomb_error( "Can not swap in sequence" );
	data[size] = 0;
}
void Sequence::SwapOut()
{
	if( swap && data ){
		delete[] data;
		bufsize = 0;
		data = NULL;
	}
}
void Sequence::PrintInfo( int id, FILE *fin, FILE *fout, const Options & options, char *buf )
{
	const char *tag = options.isEST ? "nt" : "aa";
	bool print = options.print != 0;
	bool strand = options.isEST;
#if 0
	if( options.des_len ==0 ){
		len = 0;
		while( len < MAX_DES && ! isspace( buf[len] ) ) len ++;
	}
	buf[ len ] = 0;
#endif
	fprintf( fout, "%i\t%i%s, >%s...", id, size, tag, identifier+1 );
	if( identity ){
		int *c = coverage;
		fprintf( fout, " at " );
		if (print) fprintf( fout, "%i:%i:%i:%i/", c[0], c[1], c[2], c[3] );
		if (strand) fprintf( fout, "%c/", (state & IS_MINUS_STRAND) ? '-' : '+' );
		fprintf( fout, "%.2f%%", identity*100 );
		if( options.useDistance ) fprintf( fout, "/%.2f%%", distance*100 );
		fprintf( fout, "\n" );
	}else{
		fprintf( fout, " *\n" );
	}
}
