#pragma once
#include "cdhit_util.h"


////////// Class definition //////////
class ScoreMatrix { //Matrix
	private:

	public:
		int matrix[MAX_AA][MAX_AA];
		int gap, ext_gap;

		ScoreMatrix();
		void init();
		void set_gap(int gap1, int ext_gap1);
		void set_matrix(int *mat1);
		void set_to_na();
		void set_match( int score );
		void set_mismatch( int score );
}; // END class ScoreMatrixG

//-------- extern functions -------//
extern void setaa_to_na();
