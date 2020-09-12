#include "cdhit_scomat.h"




//class function definition
const char aa[] = {"ARNDCQEGHILKMFPSTWYVBZX"};
//{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,20};
int aa2idx[] = {0, 2, 4, 3, 6, 13,7, 8, 9,20,11,10,12, 2,20,14,
                5, 1,15,16,20,19,17,20,18, 6};
    // idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
    //          Q  R  S  T  U  V  W  X  Y  Z
    // so  aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
    // and aa2idx['M'-'A'] => 12

int BLOSUM62[] = {
  4,                                                                  // A
 -1, 5,                                                               // R
 -2, 0, 6,                                                            // N
 -2,-2, 1, 6,                                                         // D
  0,-3,-3,-3, 9,                                                      // C
 -1, 1, 0, 0,-3, 5,                                                   // Q
 -1, 0, 0, 2,-4, 2, 5,                                                // E
  0,-2, 0,-1,-3,-2,-2, 6,                                             // G
 -2, 0, 1,-1,-3, 0, 0,-2, 8,                                          // H
 -1,-3,-3,-3,-1,-3,-3,-4,-3, 4,                                       // I
 -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,                                    // L
 -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,                                 // K
 -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5,                              // M
 -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,                           // F
 -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,                        // P
  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4,                     // S
  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,                  // T
 -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11,               // W
 -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,            // Y
  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,         // V
 -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,      // B
 -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,   // Z
  0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1 // X
//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
};


int na2idx[] = {0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 3, 3, 4, 4, 4, 4, 4};
    // idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P
    //          Q  R  S  T  U  V  W  X  Y  Z
    // so aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0,
    // and aa2idx['M'-'A'] => 4
int BLOSUM62_na[] = {
  2,               // A
 -2, 2,            // C
 -2,-2, 2,         // G
 -2,-2,-2, 2,      // T
 -2,-2,-2, 1, 2,   // U
 -2,-2,-2,-2,-2, 1 // N
//A  C  G  T  U  N
//0  1  2  3  3  4
};

void setaa_to_na()
{
	int i;
	for (i=0; i<26; i++) aa2idx[i]   = na2idx[i];
} // END void setaa_to_na


/////////////////
ScoreMatrix::ScoreMatrix()
{
	init();
}

void ScoreMatrix::init() 
{
	set_gap( -11, -1 );
	set_matrix( BLOSUM62 );
}

void ScoreMatrix::set_gap(int gap1, int ext_gap1)
{
	gap = MAX_SEQ * gap1;
	ext_gap = MAX_SEQ * ext_gap1;
}

void ScoreMatrix::set_matrix(int *mat1)
{
	int i, j, k;
	k = 0;
	for ( i=0; i<MAX_AA; i++)
		for ( j=0; j<=i; j++)
			matrix[j][i] = matrix[i][j] = MAX_SEQ * mat1[ k++ ];
}

void ScoreMatrix::set_to_na()
{
	set_gap( -6, -1 );
	set_matrix( BLOSUM62_na );
}
void ScoreMatrix::set_match( int score )
{
	int i;
	for ( i=0; i<MAX_AA; i++) matrix[i][i] = score;
	matrix[3][4] = matrix[4][3] = score;
}
void ScoreMatrix::set_mismatch( int score )
{
	int i, j;
	for ( i=0; i<MAX_AA; i++)
		for ( j=0; j<i; j++)
			matrix[j][i] = matrix[i][j] = MAX_SEQ * score;
	matrix[3][4] = matrix[4][3] = matrix[3][3];
}
