#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <limits.h>   //-> http://www.cplusplus.com/reference/climits/
#include "Score.h"
#include "ScoreMatrix.h"
#include "Chebyshev.h"
#include "DeepCNF_Misc.h"
using namespace std;

// ----- definitions ------//
#ifdef _FEAT_FLOAT
	typedef float FeatScore;
#else
	typedef double FeatScore;
#endif
// typedef
typedef long long int U_INT;   //-> for vector size


class DeepCNF_Model;

//====== DeepCNF_Seq class ======//
class DeepCNF_Seq
{
public:
	DeepCNF_Seq(int len, DeepCNF_Model* pModel);
	DeepCNF_Seq();
	~DeepCNF_Seq();
	int Gate_Function;            //-> 0: for sigmoid, 1: for tanh (by default)
	int sequence_length;          //-> input sequence length
	DeepCNF_Model* m_pModel;   //-> records the model parameters

//---- variables ----//
public:
	//-> model structure
	int dim_one_pos; //-> initial dimension of features (without considering bias)
	int dim_states;  //-> final dimension of states
	int tot_layer;   //-> total number of layers ( = wind_layer+2)

	//-> general data structure
	int *observed_label;
	int *predicted_label;
	ScoreMatrix *forward;
	ScoreMatrix *backward;
	ScoreMatrix *forward2;
	ScoreMatrix *backward2;

	//-> current sequence related
	U_INT *feat_dim;          //-> layer-wise feature dimension
	FeatScore **feat_orig;    //-> original feature
	FeatScore **feat_rec;     //-> expanded feature 


//----- vice functions ----//
public:
	//-> init data structure
	void Feature_Structre_Create(
		int layer_number,int *layer_windows,int *layer_nodes,int sequence_length);   //-> input for layer and window

	//-> gate function
	double AffineOut(FeatScore* features,double *weights,U_INT dim_features);
	double GateOut(double sum);
	double GateDiff(double sum);

	//-> feed-forward and back-propagation (the core function of DeepCNF)
	void FeedForward(	
		U_INT * feat_dim, FeatScore ** feat_orig, FeatScore ** feat_rec,     //-> input/output for feature dimension
		double ** layer_weight, int sequence_length,                       //-> input for layer dimension and sequence_length
		int tot_layer, int * nodes_dim, int * windows_dim);                //-> input for window dimension and total layer
	void BackPropagate(	
		U_INT * feat_dim, FeatScore ** feat_orig, FeatScore ** feat_rec,     //-> input/output for feature dimension
		double ** layer_weight, double ** layer_grad,int sequence_length,  //-> output for gradient
		int tot_layer, int * nodes_dim, int * windows_dim);                //-> input for layer and window dimension

	//-> calculate score function
	double ComputeScore(int state_num,int leftState, int currState, int pos, 
		double *state_weight, FeatScore *header_node);
	double CalcPartition(int state_num, int sequence_length, ScoreMatrix* forward );
	int GetState(int pos, int sequence_length, int *state_label);

	//-> MAP and Viterbi function
	//--> Viterbi assignment
	void ComputeViterbi(int state_num, int sequence_length, 
		double *state_weight, FeatScore *header_node, int* predicted_label);     //-> compute Viterbi assignment
	//--> forward
	void ComputeForward(int state_num,int sequence_length, 
		double *state_weight, FeatScore *header_node, ScoreMatrix* forward);     //-> compute forward probability
	//--> backward
	void ComputeBackward(int state_num,int sequence_length, 
		double *state_weight, FeatScore *header_node, ScoreMatrix* backward);    //-> compute backward probability
	//--> MAP assignment
	void ComputeMAP(int state_num, int sequence_length,
		ScoreMatrix* forward, ScoreMatrix* backward, int* predicted_label);   //-> compute MAP assignment
	//--> MAP probability
	void ComputeMAP(int state_num, int sequence_length,
		ScoreMatrix* forward, ScoreMatrix* backward, vector <vector <double> > &output);   //-> compute MAP assignment (output probability)

	//-> objective and gradient function
	double ComputeObj(int state_num,int sequence_length, 
		ScoreMatrix* forward, double *state_weight, FeatScore *header_node, int *state_label);
	void ComputeGrad(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, double *state_weight, FeatScore *header_node, int *state_label,
		double *state_grad, FeatScore *state_error);

	//-> compute test accuracy
	void ComputeTestAccuracy(int sequence_length,int *obs_label,
		int *predicted_label,double &total_num,double &correct_num);
	void ComputeTestAccuracy(double &total_num,double &correct_num);
	void ComputeTestAccuracy_Weight(int sequence_length,int *obs_label, 
		int *predicted_label,double &total_num,double &correct_num);
	void ComputeTestAccuracy_Weight(double &total_num,double &correct_num);

	//------ for maximal accuracy ------//
	void ComputeForward2(int state_num,int sequence_length, double *state_weight, FeatScore *header_node, 
		ScoreMatrix* forward, ScoreMatrix* forward2, vector <double> &q_function, vector <int> &chosen_label);
	void ComputeBackward2(int state_num,int sequence_length, double *state_weight, FeatScore *header_node, 
		ScoreMatrix* backward, ScoreMatrix* backward2, vector <double> &q_function, vector <int> &chosen_label);
	double ComputeObj2(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, double *state_weight, FeatScore *header_node, int *state_label);
	void ComputeGrad2(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
		double *state_weight, FeatScore *header_node, int *state_label, double *state_grad, FeatScore *state_error);

	//------ for maximize AUC -------//
	//-> objective function related
	double Compute_q_function(vector <vector <double> > &post_prob, 
		vector <double> &q_function,double &tot_num,
		int power_value, int *state_label,int chosen_label,
		int CHOOSEorNOT,int DERIVEorNOT);
	void ComputeObj3(vector <vector <double> > &post_prob,
		int *state_label,int chosen_label,int power_value_l,int power_value_u,
		vector <double> &q_function_s, vector <double> &q_function_v,
		double &s_value, double &v_value, double &s_value_n,double &v_value_n, int GRADorNOT);
	void ComputeObj3_Total(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, double *state_weight, FeatScore *header_node, int *state_label,
		int max_degree,double *s_value_out,double *v_value_out,double *s_value_out_n,double *v_value_out_n);
	//-> gradient function related
	void ComputeGrad3_Part(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
		double *state_weight, FeatScore *header_node, int chosen_label,vector <double> &state_grad, vector <double> &state_error,
		vector <double> &q_function, double func, double Partition);
	void ComputeGrad3(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
		double *state_weight, FeatScore *header_node, int *state_label, double *state_grad, FeatScore *state_error,
		int chosen_label, int power_value_l,int power_value_u,double gamma_value,
		double s_value_sum, double v_value_sum, double s_value_sum_n, double v_value_sum_n );
	void ComputeGrad3_Total(int state_num,int sequence_length, 
		ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
		double *state_weight, FeatScore *header_node, int *state_label, double *state_grad, FeatScore *state_error,
		int max_degree, double *gamma_coeff,double *s_value_sum, double *v_value_sum,double *s_value_sum_n, double *v_value_sum_n);


//----- main functions ----//
public:
	//-> read in features
	void Read_Init_Features(vector <string> &input);
	void Read_Init_Features(vector <vector <FeatScore> > &input);
	void Calc_Forward_Backward(void);
	void Calc_Forward_Backward2(void);  //-> for maximal accuracy

	//-> for prediction
	//--> gate output
	void Gate_Output(vector <vector <FeatScore> > &input,vector <vector <FeatScore> > &output);
	void Gate_Output(vector <string> &input,vector <vector <FeatScore> > &output);
	void Gate_Output(vector <vector <FeatScore> > &output);
	//--> MAP probability
	void MAP_Probability(vector <vector <FeatScore> > &input,vector <vector <double> > &output);
	void MAP_Probability(vector <string> &input,vector <vector <double> > &output);
	void MAP_Probability(vector <vector <double> > &output);
	//--> MAP assign
	void MAP_Assign(vector <vector <FeatScore> > &input,vector <int> &output_label);
	void MAP_Assign(vector <string> &input,vector <int> &output_label);
	void MAP_Assign(vector <int> &output_lab);
	void MAP_Assign(int *output_lab);
	void MAP_Assign(void);
	//--> Viterbi assign
	void Viterbi_Assign(vector <vector <FeatScore> > &input,vector <int> &output_label);
	void Viterbi_Assign(vector <string> &input,vector <int> &output_label);
	void Viterbi_Assign(vector <int> &output_lab);
	void Viterbi_Assign(int *output_lab);
	void Viterbi_Assign(void);

	//-> for training
	double Obj_Calc(void);
	void Grad_Calc(void);
	double Obj_Calc2(void);
	void Grad_Calc2(void);
	void Obj_Calc3(void);
	void Grad_Calc3(void);
};



//====== DeepCNF_Model class ======//
class DeepCNF_Model
{
public:
	DeepCNF_Model(int layer_number, string &window_str, string &nodes_str, int n_states, int n_local,
		int Train_Or_Not_=0,int MAX_AUC_DEGREE=0,double MAX_AUC_BETA=-1,int Gate_Function_=1);
	DeepCNF_Model(int Train_Or_Not_=0,int MAX_AUC_DEGREE=0,double MAX_AUC_BETA=-1, int Gate_Function_=1);
	~DeepCNF_Model();
	int Train_Or_Not;  //-> training swith, default is OFF [0]

//---- variables ----//
public:
	//-> gate function
	int Gate_Function;        //-> default 1 (tahn)
	//-> for actual computation
	//--> basic dimension
	int dim_one_pos;          //-> initial dimension of features (without considering bias)
	int dim_states;           //-> final dimension of states
	U_INT total_param_count;  //-> total dimension of parameters
	//--> dimension related
	int wind_layer;           //-> window layer
	int tot_layer;            //-> total number of layers ( = wind_layer+2)
	int * windows_dim;        //-> record the windows for each layer ( = wind_layer)
	int * nodes_dim;          //-> record the nodes for each layer ( = wind_layer+2)
	//--> weight and gradient
	double ** layer_weight;   //-> record layer-wise weight
	double ** layer_grad;     //-> record layer-wise gradient
	U_INT * layer_count;      //-> this data_structure is critical for Multi<->One dimension transform

	//-> dummy data structure for optimization
	double* weights;
	double* grad;             //-> for training
	double* grad_sum;         //-> for training
	//-> label weight //__2015_04_15__//
	double *label_weight;     //-> for re-weighted label optimization
	//-> maximize AUC
	int max_degree;           //-> determine the gamma_coeff (default = 0)
	double beta_sigmoid;      //-> soft AUC with parameter beta for sigmoid function
	double *gamma_coeff;      //-> for maximize AUC optimization
	double *s_value_out;      //-> for outer s_value
	double *v_value_out;      //-> for outer v_value
	double *s_value_sum;      //-> for summation s_value
	double *v_value_sum;      //-> for summation v_value
	double *s_value_out_n;    //-> number of outer s_value
	double *v_value_out_n;    //-> number of outer v_value
	double *s_value_sum_n;    //-> number of summation s_value
	double *v_value_sum_n;    //-> number of summation v_value
	double *auc_for_label;    //-> auc value for each label
	


//----- functions ----//
public:
	//-> collect gradsum
	void Collect_GradSum(double *grad_sum);
	void Collect_GradSum(void);

	//-> init data structure
	void Weight_Structre_Create(
		int layer_number,string &windows_str,string &nodes_str,int n_states, int n_local);  //-> input for layer and window

	//-> for dimension transformation
	void MultiDim_To_OneDim(int tot_layer,U_INT * layer_count,double ** layer_weight,double *out_weight);
	void OneDim_To_MultiDim(int tot_layer,U_INT * layer_count,double ** layer_weight,double *in_weight);

	//-> for loading model
	void cnf_model_load(string &file);                 //model load from file
	void cnf_model_load(double *input,U_INT totnum);   //model load from data-structure

	//-> for predicting
	void Pred_Single_Chain(vector <string> &input,vector <vector <double> > &output);   //-> for predicting purpose

	//-> max AUC related
	int init_gamma_coeff(double *gamma_coeff,int maxdeg,double beta);
	void Initialize_SVvalue(void);
	void SumOver_SVvalue(void);
	double Calculate_AUC_Value(void);

};
