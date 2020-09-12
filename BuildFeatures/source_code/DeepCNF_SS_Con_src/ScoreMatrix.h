#pragma once
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "Score.h"
using namespace std;


class ScoreMatrix
{

	int layers;
	int rows;
	int cols;
	Score *data;  

public:
	ScoreMatrix (){data=0; rows=cols=0;};
	ScoreMatrix (int rows, int cols){data=0;resize(rows, cols);};
	ScoreMatrix (const ScoreMatrix &m) : rows(m.rows), cols(m.cols), data(m.data){};
	~ScoreMatrix(){delete data;};

	void Fill (const Score &value){  for (int i = 0; i < rows * cols; i++) data[i] = value;};

	inline void resize(int rows, int cols)
	{
		if (0==rows || 0==cols)
		{
			cerr << " ERROR: " << rows << " rows; " << cols << " cols." << endl;
			exit(0);
		}
		this->rows = rows; this->cols=cols;
		if (data) delete data;
		data = new Score[rows*cols];
	}

	inline Score &operator() (int row, int col)
	{
		if (row>=rows || col>=cols)
		{
			cerr << "ERROR: ScoreMatrix " << row << " vs. " << rows << "; " << col << " vs. " << cols<< endl;
			exit(0);
		}
		return data[(row * cols + col)];
	}

	inline Score operator() (int row, int col) const 
	{
		if (row>=rows || col>=cols)
		{
			cerr << "ERROR: ScoreMatrix" << row << " vs. " << rows << "; " << col << " vs. " << cols<< endl;
			exit(0);
		}
		return data[(row * cols + col)];
	}
};
