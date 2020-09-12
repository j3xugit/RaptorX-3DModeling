#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "CNF_calc.h"

/*
//---- constact_param ----//
int num_procs=1; 
int proc_id=0;
*/

//----- start -----//
SEQUENCE_cnf::SEQUENCE_cnf(int length_seq, bCNF_Model* pModel)
{
	m_pModel = pModel;
	this->length_seq = length_seq;
	forward = new ScoreMatrix(m_pModel->num_states,length_seq);
	backward = new ScoreMatrix(m_pModel->num_states,length_seq);
	obs_label = new int[length_seq];
	int df = m_pModel->dim_features*length_seq;
	
	_features = new Score[df];

	predicted_label = new int[length_seq];	
	obs_feature = new Score*[length_seq];
	for(int i=0;i<length_seq;i++)
		obs_feature[i] = new Score[m_pModel->dim_one_pos];
}


SEQUENCE_cnf::~SEQUENCE_cnf(){
	delete forward;
	delete backward;
	delete obs_label;
	delete _features;

	for(int i=0;i<length_seq;i++)
		delete obs_feature[i];
	delete obs_feature;
}

Score SEQUENCE_cnf::Obj_p(){ // predicted label : log-likelihood
	ComputeForward();
	ComputeBackward();
	CalcPartition();
	Score obj = -Partition;
	
	for(int t=0;t<length_seq;t++){
		int leftState = DUMMY;
		if(t>0) leftState = predicted_label[t-1];
		int currState = predicted_label[t];
		obj+=ComputeScore(leftState,currState,t);
	}
	cout << " predicted obj = " << obj << endl;
	return obj;
}

Score SEQUENCE_cnf::Obj(){ // observed label : log-likelihood`
	ComputeGates();
	ComputeVi();
	ComputeForward();
	//ComputeBackward();
	CalcPartition();
	Score obj = -Partition;

	for(int t=0;t<length_seq;t++){
		int leftState = GetObsState(t-1);
		int currState = GetObsState(t);
		obj+=ComputeScore(leftState,currState,t);
	}
	
	return obj;
}

void SEQUENCE_cnf::ComputeViterbi()
{
	// Viterbi Matrix
	ScoreMatrix best(m_pModel->num_states,length_seq);
	best.Fill((Score)LogScore_ZERO);
	// TraceBack Matrix
	ScoreMatrix traceback(m_pModel->num_states,length_seq);
	traceback.Fill(DUMMY);
	// compute the scores for the first position
	for(int i=0;i<m_pModel->num_states;i++){
		best(i,0)=ComputeScore(DUMMY,i,0);
	}

	// Compute the Viterbi Matrix
	for(int t=1;t<length_seq;t++){
		for(int currState=0;currState<m_pModel->num_states;currState++){
			for(int leftState=0;leftState<m_pModel->num_states;leftState++){
				Score new_score = ComputeScore(leftState,currState,t) + best(leftState,t-1);
				if(new_score > best(currState,t)){
					best(currState,t) = new_score;
					traceback(currState,t) = leftState;
				}
			}
		}
	}
	
	Score max_s = LogScore_ZERO;
	int last_state = 0;
	for(int i=0;i<m_pModel->num_states;i++)
		if(best(i,length_seq-1)>max_s)
			max_s = best(i,length_seq-1), last_state = i;
	//TraceBack
	for(int t=length_seq-1; t>=0;t--){
		predicted_label[t]=last_state;
		last_state=(int)traceback(last_state,t);
	}
}

void SEQUENCE_cnf::MAP(){ // Posterior Decoding (Marginal Probability Decoder)
	ComputeForward();
	ComputeBackward();
	for(int t=0;t<length_seq;t++){
		int idx = 0;
		Score maxS = LogScore_ZERO;
		for(int i=0;i<m_pModel->num_states;i++){
			Score s = 0;
			//if(t==0) 
			//	s = (*backward)(i,0); 
			//else 
				s = (*backward)(i,t) + (*forward)(i,t);
			if(s > maxS) maxS = s,idx = i;
		}
		predicted_label[t]=idx;
	}
}

void SEQUENCE_cnf::MAP(vector <vector <double> > &output){ // Posterior Decoding (Marginal Probability Decoder)
	ComputeForward();
	ComputeBackward();
	CalcPartition();
	for(int t=0;t<length_seq;t++){
		int idx = 0;
		Score maxS = LogScore_ZERO;
		double wstot=0;
		for(int i=0;i<m_pModel->num_states;i++){
			Score s = 0;
			//if(t==0) 
			//	s = (*backward)(i,0); 
			//else 
//				s = (*backward)(i,t) + (*forward)(i,t);
//			if(s > maxS) maxS = s,idx = i;
          double post=
                      exp(
                        (
                         (*forward)(i,t) +
                         (*backward)(i,t)
                          - Partition
                        ));
		  //assign
		  if(post>1)post=1;
		  if(post<0)post=0;
		  output[t][i]=post;
		  wstot+=post;
		}
//		predicted_label[t]=idx;
		//normalize
		if(wstot==0)for(int i=0;i<m_pModel->num_states;i++)output[t][i]=1.0/m_pModel->num_states;
		else for(int i=0;i<m_pModel->num_states;i++)output[t][i]=1.0*output[t][i]/wstot;
	}
}

void SEQUENCE_cnf::ComputeForward()
{
	forward->Fill(LogScore_ZERO);
	for(int i=0;i<m_pModel->num_states;i++){
		(*forward)(i,0)=ComputeScore(DUMMY,i,0);
	}

	for(int t=1;t<length_seq;t++){
		for(int currState=0;currState<m_pModel->num_states;currState++){
			for(int leftState=0;leftState<m_pModel->num_states;leftState++) {
				Score new_score = ComputeScore(leftState,currState,t);
				LogScore_PLUS_EQUALS((*forward)(currState,t),new_score + (*forward)(leftState,t-1));
			}
		}
	}
}

void SEQUENCE_cnf::ComputeBackward()
{
	backward->Fill(LogScore_ZERO);
	for(int i=0;i<m_pModel->num_states;i++){
		(*backward)(i,length_seq-1)=0;
	}

	for(int t=length_seq-2;t>=0;t--){
		for(int currState=0;currState<m_pModel->num_states;currState++){
			for(int rightState=0;rightState<m_pModel->num_states;rightState++){
				Score new_score = ComputeScore(currState,rightState,t+1);
				LogScore_PLUS_EQUALS((*backward)(currState,t),new_score + (*backward)(rightState,t+1));
			}
		}
	}
}

void SEQUENCE_cnf::ComputeGates()
{
	gates.resize(m_pModel->num_gates, length_seq);
	for (int pos=0;pos<length_seq; pos++) {
		Score* features=getFeatures(pos);
		for (int k=0;k<m_pModel->num_gates;k++) {
			gates(k,pos) = m_pModel->GetGateOutput(k,features);
		}
	}
}

void SEQUENCE_cnf::ComputeVi()
{
	arrVi.resize(m_pModel->num_states, length_seq);
	int num_states = m_pModel->num_states;
	int num_gates = m_pModel->num_gates;
	for (int pos=0;pos<length_seq; pos++)
		for (int currState=0;currState<m_pModel->num_states;currState++) {
			int weightLStart = num_states + num_states*num_states + currState*num_gates;
			Score score = 0;
			for(int k=0;k<m_pModel->num_gates;k++) {
				Score output = gates(k,pos);
				score += m_pModel->weights[weightLStart++]*output;
			}
			arrVi(currState, pos) = score;
		}
}

Score SEQUENCE_cnf::ComputeScore(int leftState, int currState, int pos){
	//compute num_gates neurons' output
	Score score = m_pModel->GetWeightT(leftState,currState)+arrVi(currState, pos);
	return score;
}

void SEQUENCE_cnf::makeFeatures(){
	//build features for local windows	
	for(int t=0;t<length_seq;t++){
		int pivot = t*m_pModel->dim_features;
		for(int i=0; i< m_pModel->window_size; i++){
			int offset = t+i-m_pModel->window_size/2;
			if(offset <0 || offset >=length_seq){
				for(int j=0;j<m_pModel->dim_one_pos;j++)
					_features[pivot] = 0, pivot++;
			} else {
				for(int j=0;j<m_pModel->dim_one_pos;j++)
					_features[pivot] = obs_feature[offset][j], pivot++;
			}
		}
		_features[pivot] = 1; // last one for the bias term
	}
}

void SEQUENCE_cnf::CalcPartition(){
	Partition = (Score)LogScore_ZERO;
	Score FP = (Score)LogScore_ZERO;
	//Score BP = (Score)LogScore_ZERO;
	for(int k=0;k<m_pModel->num_states;k++)
		LogScore_PLUS_EQUALS(Partition, (*forward)(k,length_seq-1));
	//for(int k=0;k<m_pModel->num_states;k++)
	//	LogScore_PLUS_EQUALS(BP, (*backward)(k,0) + (*forward)(k,0));
}


Score* SEQUENCE_cnf::getFeatures(int pos){
	int offset;
	offset = pos*m_pModel->dim_features;
	return _features+offset;
}

int SEQUENCE_cnf::GetObsState(int pos){
	if(pos<0 || pos>=length_seq) return DUMMY;
	return obs_label[pos];
}

int bypass;
double regularizer;
int num_test_set;


bCNF_Model::bCNF_Model(void)
{
	grad=NULL;
	weights=NULL;
	grad_sum=NULL;
	reg=NULL;
}
bCNF_Model::~bCNF_Model(void)
{
	if(grad!=NULL)delete [] grad;
	if(weights!=NULL)delete [] weights;
	if(grad_sum!=NULL)delete [] grad_sum;
	if(reg!=NULL)delete [] reg;
}


void bCNF_Model::SetParameters(int w_size, int n_states, int n_gates, int n_local)
{
	//setSeed();
	window_size = w_size; //best 13

	num_states = n_states;
	num_gates = n_gates; //best 20
	dim_one_pos = n_local;  // dim of local feature + 1 (for bias)

	dim_features = (window_size*dim_one_pos)+1;
	
	num_params = num_states*(1 + num_states + num_gates) + dim_features*num_gates;

	weights = new double[num_params];
	grad = new double[num_params];	
	grad_sum = new double[num_params];
	reg = new double[num_params];	
	double r=regularizer; // regularization coefficients
	int i;
	for(i=0;i<num_states*(1+num_states);i++) 
		reg[i] = r;
	for(;i<num_states*(1+num_states+num_gates);i++) 
		reg[i] = r;
	for(;i<num_params;i++)
		reg[i] = r;
}

//------- data_load -------//
void bCNF_Model::LoadData_II(string &file)
{
	ifstream fin;
	string buf,temp;
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",file.c_str());
		exit(-1);
	}
	vector<SEQUENCE_cnf*> DATA;
	DATA.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream trn_in(buf);
		int length_s=1; //only one node
		SEQUENCE_cnf *seq = new SEQUENCE_cnf(length_s, this);
		for(int j=0;j<length_s;j++){
			for(int k=0;k<dim_one_pos-1;k++){
				trn_in >> seq->obs_feature[j][k];
			}
			seq->obs_feature[j][dim_one_pos-1] = 1;
		}
		DATA.push_back(seq);
	}
	//push back
	testData.clear();
	for(int i=0;i<(int)DATA.size();i++)testData.push_back(DATA[i]);
}
void bCNF_Model::LoadData_III(vector <string> &input)
{
	vector<SEQUENCE_cnf*> DATA;
	DATA.clear();
	for(int i=0;i<(int)input.size();i++)
	{
		istringstream trn_in(input[i]);
		int length_s=1; //only one node
		SEQUENCE_cnf *seq = new SEQUENCE_cnf(length_s, this);
		for(int j=0;j<length_s;j++){
			for(int k=0;k<dim_one_pos-1;k++){
				trn_in >> seq->obs_feature[j][k];
			}
			seq->obs_feature[j][dim_one_pos-1] = 1;
		}
		DATA.push_back(seq);
	}
	//push back
	testData.clear();
	for(int i=0;i<(int)DATA.size();i++)testData.push_back(DATA[i]);
}

//----------- for normal chain CNF ------------//
void bCNF_Model::LoadData_IV(vector <string> &input)  //-> for single chain
{
	testData.clear();
	//load in
	int length_s=(int)input.size(); //single chain
	SEQUENCE_cnf *seq = new SEQUENCE_cnf(length_s, this);
	for(int j=0;j<length_s;j++)
	{
		istringstream trn_in(input[j]);
		for(int k=0;k<dim_one_pos-1;k++)
		{
			trn_in >> seq->obs_feature[j][k];
		}
		seq->obs_feature[j][dim_one_pos-1] = 1;
	}
	testData.push_back(seq);
}

//------- model_load --------//
void bCNF_Model::cnf_model_load(double *weight,int states,int windows,int feat_num,int gates)
{
	//assign parameter
	window_size = windows; //for single node, this is always one
	num_states = states;   //now is 20
	num_gates = gates;     //now is 20
	dim_one_pos = feat_num;  // dim of local feature + 1 (for bias)
	dim_features = (window_size*dim_one_pos)+1;
	num_params = num_states*(1 + num_states + num_gates) + dim_features*num_gates;
	//assign weight
	weights = new double[num_params];
	for(int i=0;i<num_params;i++)weights[i]=weight[i];
}

//the requirement is: states, windows, feat_num, gates
void bCNF_Model::cnf_model_load(string &file,int states,int windows,int feat_num,int gates)
{
	ifstream fin;
	string buf,temp;
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",file.c_str());
		exit(-1);
	}
	//skip header
/*
	//[1] num parameter
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s file format bad!\n",file.c_str());
		exit(-1);
	}
	{
		istringstream www(buf);
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		num_params=atof(temp.c_str());
	}
	//[2] num_gates
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s file format bad!\n",file.c_str());
		exit(-1);
	}
	{
		istringstream www(buf);
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		num_gates=atof(temp.c_str());
	}
	//[3] window_size
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s file format bad!\n",file.c_str());
		exit(-1);
	}
	{
		istringstream www(buf);
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		window_size=atof(temp.c_str());
	}
	//[4] dim_features
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s file format bad!\n",file.c_str());
		exit(-1);
	}
	{
		istringstream www(buf);
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		if(!(www>>temp))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
		dim_features=atof(temp.c_str());
	}
*/
	//assign parameter
	window_size = windows; //for single node, this is always one
	num_states = states;   //now is 20
	num_gates = gates;     //now is 20
	dim_one_pos = feat_num;  // dim of local feature + 1 (for bias)
	dim_features = (window_size*dim_one_pos)+1;
	num_params = num_states*(1 + num_states + num_gates) + dim_features*num_gates;

	//skip
	int i;
	for(i=0;i<6;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"%s file format bad!\n",file.c_str());
			exit(-1);
		}
	}
	//get feature
	vector<double> w_init;
	{
		istringstream www(buf);
		for(;;)
		{
			if(!(www>>temp))break;
			double value=atof(temp.c_str());
			w_init.push_back(value);
		}
	}
	//final assign
	if(w_init.size()!=num_params)
	{
		fprintf(stderr,"%s number of params not equal !! [%d!=%d]\n",file.c_str(),
			w_init.size(),num_params);
		exit(-1);
	}
	weights = new double[num_params];
	for(i=0;i<(int)w_init.size();i++)weights[i]=w_init[i];
}


/*
//------------ main -----------//
int main(int argc, char **argv)
{
	//---------------- CNF_calc ------------//
	{
		if(argc<4)
		{
			printf("CNF_calc <feat_file> <model_file> <out_file> \n");
			printf("   [states=20] [feat_num=131] [gates=20] \n");
			printf("[note]: feat_file must in NN format \n");
			exit(-1);
		}
		string feat_file=argv[1];
		string model_file=argv[2];
		string out_file=argv[3];
		int states=20;    //20 amino acid
		int feat_num=131; //real_feat + 1
		int gates=20;
		if(argc>4)
		{
			if(argv[4]!=">" && argv[4]!=">>")
			{
				states=atoi(argv[4]);
				if(argc>5)
				{
					if(argv[5]!=">" && argv[5]!=">>")
					{
						feat_num=atoi(argv[5]);
						if(argc>6)
						{
							if(argv[6]!=">" && argv[6]!=">>")
							{
								gates=atoi(argv[6]);
							}
						}
					}
				}
			}
		}

		//--- load model ---//
		bCNF_Model cnfModel;
		int windows=1;    //for single node, always 1
		cnfModel.cnf_model_load(model_file,states,windows,feat_num,gates);
		cnfModel.LoadData_II(feat_file);
		//-- calc --//
		int i,j;
		int size=(int)cnfModel.testData.size();
		vector <vector <double> > output_full (size, vector <double> (states));
		for(i=0;i<size;i++)
		{
			vector <vector <double> > output (1, vector <double> (states));
			cnfModel.testData[i]->makeFeatures();
			cnfModel.testData[i]->ComputeGates();
			cnfModel.testData[i]->ComputeVi();
			cnfModel.testData[i]->MAP(output);
			for(j=0;j<20;j++)output_full[i][j]=output[0][j];
		}
		//output
		FILE *fp=fopen(out_file.c_str(),"wb");
		for(i=0;i<size;i++)
		{
			for(j=0;j<states;j++)fprintf(fp,"%5.2f ",log(1.0*output_full[i][j]/_aa_freq_17W[j]));
			fprintf(fp,"\n");
		}
		fclose(fp);
		//exit
		exit(0);
	}
}
*/
