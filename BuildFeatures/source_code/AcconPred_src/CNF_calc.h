#pragma once
#include "ScoreMatrix.h"
#include "Score.h"
//#include "LBFGS.h"

using namespace std;

#define DUMMY -1

class bCNF_Model;

class SEQUENCE_cnf
{
public:
	SEQUENCE_cnf(int len, bCNF_Model* pModel);
	~SEQUENCE_cnf();

	bCNF_Model* m_pModel;

	int length_seq;

	int* obs_label;
	Score *_features;
	Score **obs_feature;

	Score Partition;
	int *predicted_label;
	ScoreMatrix *forward;
	ScoreMatrix *backward;

	ScoreMatrix gates;
	ScoreMatrix arrVi;
	void ComputeGates();
	void ComputeVi();

	void ComputeViterbi();
	void ComputeForward();
	void ComputeBackward();
	void CalcPartition();

	void ComputePartition();
	Score ComputeScore(int leftState, int currState, int pos);

	void makeFeatures();
	Score* getFeatures(int pos);

	int GetObsState(int pos);
	void ComputeGradient(bool bCalculateGate);
	void MAP();	
	void MAP(vector <vector <double> > &output); //calculate probability of each states
	
	void ComputeTestAccuracy();

	Score Obj();
	Score Obj_p();
};

/*
class _LBFGS: public Minimizer
{
public:
	_LBFGS(bCNF_Model* pModel): Minimizer(false) {m_pModel = pModel;};

	bCNF_Model* m_pModel;

	void ComputeGradient(vector<double> &g, const vector<double> &x, bool bCalculateGate);
	double ComputeFunction(const vector<double> &x);
	void Report(const vector<double> &theta, int iteration, double objective, double step_length);
	void Report(const string &s);
};
*/

class bCNF_Model
{
public:
	bCNF_Model();
	~bCNF_Model();
public:
	int num_states;
	int num_data;
	int num_tst;
	int num_gates;
	int dim_one_pos;
	int dim_features;
	int window_size;
	int num_params;
	int totalPos;
	int totalCorrect;

	string model_file;

	double apw;
	int ct;

	vector<SEQUENCE_cnf*> trainData,testData;

	double* grad;
	double* weights;
	double* grad_sum;
	double* reg;

	// get the transition prob, num_states( for DUMMY ) + num_states*num_states
	inline Score GetWeightT(int leftState, int currState){return weights[num_states + leftState*num_states + currState];}

	// get the label-based weight, num_states*num_gates
	inline Score GetWeightL(int currState, int gate){ return weights[num_states + num_states*num_states + currState*num_gates + gate];}

	// get the Gate weight, num_gates*dim_features
	inline Score GetWeightG(int gate, int dim){ return weights[num_states*(num_states+num_gates + 1) + dim_features*gate + dim];}
	void SetSeed();
	void SetParameters(int w_size, int n_states, int n_gates, int n_local);
	void Initialize(string model_dir, int w_size, int n_states, int n_gates, int n_local, string output_file, string input_f);
	void LoadData(string input);
	Score Gate(Score sum){ return (Score)1.0/(1.0+exp(-(double)sum)); }
	Score GetGateOutput(int gate, Score* features){
		Score output = 0;
		//dim_features
		int weightGStart = num_states*(num_states+num_gates + 1) + dim_features*gate;
		for(int j=0;j<dim_features;j++){
			//output+=GetWeightG(gate,j)*features[j];
			output+=weights[weightGStart++]*features[j];
		}
		return Gate(output);
	}

	//just for predicting 
	void cnf_model_load(string &file,int states,int windows,int feat_num,int gates);  //model load
	void cnf_model_load(double *weight,int states,int windows,int feat_num,int gates);  //model load
	void LoadData_II(string &file);     //only for predicting //added by WS //__2012_06_30__//
	void LoadData_III(vector <string> &input); //__2012_06_35__//
	void LoadData_IV(vector <string> &input);  //__2012_08_20__//

};
