#include "DeepCNF.h"


//---- Deep Convolution CNF architecture ----//-> take layer=2 as example
/*
            iiiiiiiiii    -> input feature layer (input layer) [window 1] -> the bottom layer
            \/\/\/\/\/    -> Weight 1 (input_to_1)
             xxxxxxxx     -> neuron feature layer (layer 1) [window 2]
             \/\/\/\/     -> Weight 2 (1_to_2)
              xxxxxx      -> neuron feature layer (layer 2) [no window]   -> the header layer
              ||||||      -> Weight 2_to_label
              LLLLLL      -> label layer (output layer)
              ------      -> Weight label_to_label
*/

//========== constructor and destructor ============//
DeepCNF_Seq::DeepCNF_Seq(int len, DeepCNF_Model* pModel)
{
	//-> gate function
	Gate_Function = pModel->Gate_Function;   //-> default is tanh 
	
	//-> basic assign
	m_pModel = pModel;
	sequence_length = len;

	//-> model assign
	dim_one_pos = pModel->dim_one_pos;
	dim_states = pModel->dim_states;
	tot_layer = pModel->tot_layer;

	//-> general data structure
	observed_label = new int[len];
	predicted_label = new int[len];
	forward = new ScoreMatrix(pModel->dim_states,len);
	backward = new ScoreMatrix(pModel->dim_states,len);
	forward2 = new ScoreMatrix(pModel->dim_states,len);
	backward2 = new ScoreMatrix(pModel->dim_states,len);

	//-> current sequence related
	Feature_Structre_Create(pModel->tot_layer,pModel->windows_dim,pModel->nodes_dim,len);
}
DeepCNF_Seq::DeepCNF_Seq(void)
{
	//-> gate function
	Gate_Function = 1;   //-> default is tanh 

	//-> general data structure
	observed_label=NULL;
	predicted_label=NULL;
	forward=NULL;
	backward=NULL;
	forward2=NULL;
	backward2=NULL;

	//-> current sequence related
	feat_dim=NULL;
	feat_orig=NULL;
	feat_rec=NULL;
}
DeepCNF_Seq::~DeepCNF_Seq(void)
{
	//-> general data structure
	if(observed_label!=NULL)delete [] observed_label;
	if(predicted_label!=NULL)delete [] predicted_label;
	if(forward!=NULL)delete forward;
	if(backward!=NULL)delete backward;
	if(forward2!=NULL)delete forward2;
	if(backward2!=NULL)delete backward2;

	//-> current sequence related
	if(feat_dim!=NULL)delete [] feat_dim;
	if(feat_orig!=NULL)
	{
		for(int i=0;i<tot_layer;i++)delete [] feat_orig[i];
		delete [] feat_orig;
	}
	if(feat_rec!=NULL)
	{
		for(int i=0;i<tot_layer;i++)delete [] feat_rec[i];
		delete [] feat_rec;
	}
}


//=================== feature structure create for each input sequence ==================//
void DeepCNF_Seq::Feature_Structre_Create(
	int layer_number,int *layer_windows,int *layer_nodes,int sequence_length)   //-> input for layer and window
{
	// ------- feature related ----------//
	// feature dimension (we should include the 'bias' term !!! )
	U_INT *feature_dimension=new U_INT[layer_number];
	for(int i=0;i<layer_number-2;i++)feature_dimension[i]=(layer_nodes[i]+1)*(2*layer_windows[i]+1)+1; // (we include bias term here)
	feature_dimension[layer_number-2]=layer_nodes[layer_number-2];           //-> N for header layer (we exclude bias term here)
	feature_dimension[layer_number-1]=layer_nodes[layer_number-1];           //-> N+1 for state layer

	// record original feature
	FeatScore **feature_orig=new FeatScore*[layer_number];
	for(int i=0;i<layer_number;i++)feature_orig[i]=new FeatScore[layer_nodes[i]*sequence_length];
	// record expanded feature 
	FeatScore **feature_record=new FeatScore*[layer_number];
	for(int i=0;i<layer_number;i++)feature_record[i]=new FeatScore[feature_dimension[i]*sequence_length];

	//---- final assign ----//
	//-> feat out
	feat_dim=feature_dimension;       //-> layer_number+2   //-> this datastructure is very important
	feat_rec=feature_record;          //-> layer_number+2
	feat_orig=feature_orig;           //-> layer_number+2
}


//=========================== affine and gate function ======================//
//-> affine function
double DeepCNF_Seq::AffineOut(FeatScore* features,double *weights,U_INT dim_features)
{
	double output = 0.0;
	for(U_INT j=0;j<dim_features;j++)output+=weights[j]*features[j];
	return output;
}
//-> gate function
double DeepCNF_Seq::GateOut(double sum)
{
	if(Gate_Function==1)       //-> tanh 
	{
		double value=exp(-2.0*sum);
		return (double)( (1.0-value)/(1.0+value) );
	}
	else if(Gate_Function==2)  //-> sigmoid
	{
		double value=exp(-1.0*sum);
		return (double)( 1.0/(1.0+value) ); 
	}
	else                       //-> relu
	{
		if(sum>0)return sum;
		else return 0;
	}
}
//-> gate derivative
double DeepCNF_Seq::GateDiff(double sum)
{
	if(Gate_Function==1)      //-> tanh 
	{
		return 1.0-sum*sum;
	}
	else if(Gate_Function==2) //-> sigmoid
	{
		return (1.0-sum)*sum;
	}
	else                      //-> relu
	{
		if(sum>0)return 1;
		else return 0;
	}
}


//========================== feed-forward and back-propagation function ======================//

//------------ feed forward function -------------//
void DeepCNF_Seq::FeedForward(	
	U_INT * feat_dim, FeatScore **feat_orig, FeatScore **feat_rec,    //-> input/output for feature dimension
	double ** layer_weight, int sequence_length,                    //-> input for layer dimension and sequence_length
	int tot_layer, int * nodes_dim, int * windows_dim )             //-> input for window dimension and total layer
{
	int wind_layer=tot_layer-2;
	//--------- bottom feature layer -----------//
	{
		int l=0;
		for(int t=0;t<sequence_length;t++)
		{
			U_INT pivot=feat_dim[l]*t;
			for(int i=0; i< 2*windows_dim[l]+1; i++)
			{
				int offset = t+i-windows_dim[l];
				if(offset <0 || offset >=sequence_length)
				{
					for(int j=0;j<nodes_dim[l];j++)feat_rec[l][pivot] = 0, pivot++;
					feat_rec[l][pivot] = 0, pivot++; //-> bias term
				} 
				else 
				{
					for(int j=0;j<nodes_dim[l];j++)feat_rec[l][pivot] = feat_orig[l][offset*nodes_dim[l]+j], pivot++;
					feat_rec[l][pivot] = 1, pivot++; //-> bias term
				}
			}
			feat_rec[l][pivot] = 1; // last one for the bias term
		}
	}

	//--------- feedforward from 1 to N layers ------------//
	for(int l=1;l<=wind_layer;l++)
	{
		//----- calculate the original feature of the current layer from the previous layer of the expanded features ----//
		for(int t=0;t<sequence_length;t++)
		{
			U_INT pivot=nodes_dim[l]*t;
			FeatScore *features=feat_rec[l-1]+t*feat_dim[l-1];
			for (int j=0;j<nodes_dim[l];j++)
			{
				double *weights=layer_weight[l-1]+j*feat_dim[l-1];
				double sumover=AffineOut(features,weights,feat_dim[l-1]);
				feat_orig[l][pivot+j]=(FeatScore)GateOut(sumover);
			}
		}

		//----- consider windows to calculate the expanded features at the current layer -----//
		if(l<wind_layer) //--> for 1 to N-1 layers, considering windows
		{
			for(int t=0;t<sequence_length;t++)
			{
				U_INT pivot=feat_dim[l]*t;
				for(int i=0; i< 2*windows_dim[l]+1; i++)
				{
					int offset = t+i-windows_dim[l];
					if(offset <0 || offset >=sequence_length)
					{
						for(int j=0;j<nodes_dim[l];j++)feat_rec[l][pivot] = 0, pivot++;
						feat_rec[l][pivot] = 0, pivot++; //-> bias term
					} 
					else 
					{
						for(int j=0;j<nodes_dim[l];j++)feat_rec[l][pivot] = feat_orig[l][offset*nodes_dim[l]+j], pivot++;
						feat_rec[l][pivot] = 1, pivot++; //-> bias term
					}
				}
				feat_rec[l][pivot] = 1; // last one for the bias term
			}
		}
		else    //--> for the header layer (say N-th layer), without considering windows
		{
			for(int t=0;t<sequence_length;t++)
			{
				U_INT pivot=feat_dim[l]*t;
				for(int j=0;j<nodes_dim[l];j++)feat_rec[l][pivot] = feat_orig[l][t*nodes_dim[l]+j], pivot++;
			}
		}
	}

	//---- feedforward the state layer ----//
	{
		int l=wind_layer+1;
		for(int t=0;t<sequence_length;t++)
		{
			U_INT pivot=nodes_dim[l]*t;
			FeatScore *features=feat_rec[l-1]+t*feat_dim[l-1];
			for (int j=0;j<nodes_dim[l];j++)
			{
				double *weights=layer_weight[l-1]+j*feat_dim[l-1];
				double sumover=AffineOut(features,weights,feat_dim[l-1]);
				feat_orig[l][pivot+j]=(FeatScore)sumover;
				feat_rec[l][pivot+j]=(FeatScore)sumover;
			}
		}
	}
}

//------------ back propagation function -------------//
void DeepCNF_Seq::BackPropagate(	
	U_INT * feat_dim, FeatScore ** feat_orig, FeatScore ** feat_rec,           //-> input/output for feature dimension
	double ** layer_weight, double ** layer_grad,                            //-> output for gradient
	int sequence_length, int tot_layer, int * nodes_dim, int * windows_dim)  //-> input for layer and window dimension
{
	int wind_layer=tot_layer-2;
	//------- propagate state layer error -------//
	{
		int l=wind_layer+1;
		//-> init
		for(U_INT k=0;k<nodes_dim[l]*feat_dim[l-1];k++)layer_grad[l-1][k]=0;
		//-> calculate
		for(int t=0;t<sequence_length;t++)
		{
			U_INT pivot=nodes_dim[l]*t;
			FeatScore *features=feat_rec[l-1]+t*feat_dim[l-1];
			//--> calculate gradient
			for (int j=0;j<nodes_dim[l];j++)
			{
				FeatScore cur_error=feat_orig[l][pivot+j];        //-> final error should be put here !!
				double *grad=layer_grad[l-1]+j*feat_dim[l-1];
				for(U_INT k=0;k<feat_dim[l-1];k++)grad[k]+=cur_error*features[k];
			}
			//--> calculate expanded error
			vector <FeatScore> tmp_rec(feat_dim[l-1],0.0);
			for (int j=0;j<nodes_dim[l];j++)
			{
				FeatScore cur_error=feat_orig[l][pivot+j];
				double *weight=layer_weight[l-1]+j*feat_dim[l-1];
				for(U_INT k=0;k<feat_dim[l-1];k++)tmp_rec[k]+=(FeatScore)cur_error*weight[k]*GateDiff(features[k]);
			}
			for(U_INT k=0;k<feat_dim[l-1];k++)features[k]=tmp_rec[k];
		}
	}

	//--------- propagate from N to 1 layers ------------//
	for(int l=wind_layer;l>=1;l--)
	{
		// ------------- merge feat_rec to feat_orig to combine all the error ------------//
		//-> init
		for(U_INT k=0;k<nodes_dim[l]*sequence_length;k++)feat_orig[l][k]=0;
		for(U_INT k=0;k<nodes_dim[l]*feat_dim[l-1];k++)layer_grad[l-1][k]=0;
		//-> consider windows to calculate the expanded features at the current layer
		if(l<wind_layer) //--> for 1 to N-1 layers, considering windows
		{
			for(int t=0;t<sequence_length;t++)
			{
				U_INT pivot=feat_dim[l]*t;
				for(int i=0; i< 2*windows_dim[l]+1; i++)
				{
					int offset = t+i-windows_dim[l];
					if(offset <0 || offset >=sequence_length)
					{
						for(int j=0;j<nodes_dim[l];j++)pivot++;
						pivot++; //-> bias term
					} 
					else 
					{
						for(int j=0;j<nodes_dim[l];j++) feat_orig[l][offset*nodes_dim[l]+j] += feat_rec[l][pivot], pivot++;
						pivot++; //-> bias term
					}
				}
			}
		}
		else    //--> for the header layer (say N-th layer), without considering windows
		{
			for(int t=0;t<sequence_length;t++)
			{
				U_INT pivot=feat_dim[l]*t;
				for(int j=0;j<nodes_dim[l];j++) feat_orig[l][t*nodes_dim[l]+j] += feat_rec[l][pivot], pivot++;
			}
		}

		// ------------- given the merged error, propagate to the next layer ------------//
		for(int t=0;t<sequence_length;t++)
		{
			U_INT pivot=nodes_dim[l]*t;
			FeatScore *features=feat_rec[l-1]+t*feat_dim[l-1];
			//--> calculate gradient
			for (int j=0;j<nodes_dim[l];j++)
			{
				FeatScore cur_error=feat_orig[l][pivot+j];
				double *grad=layer_grad[l-1]+j*feat_dim[l-1];
				for(U_INT k=0;k<feat_dim[l-1];k++)grad[k]+=cur_error*features[k];
			}
			//--> calculate expanded error
			if(l>1)
			{
				vector <FeatScore> tmp_rec(feat_dim[l-1],0.0);
				for (int j=0;j<nodes_dim[l];j++)
				{
					FeatScore cur_error=feat_orig[l][pivot+j];
					double *weight=layer_weight[l-1]+j*feat_dim[l-1];
					for(U_INT k=0;k<feat_dim[l-1];k++)tmp_rec[k]+=(FeatScore)cur_error*weight[k]*GateDiff(features[k]);
				}
				for(U_INT k=0;k<feat_dim[l-1];k++)features[k]=tmp_rec[k];
			}
		}
	}
}



//=============== calculate Function and Error ==========//

//-------- compute score related ----------//
//-> compute score: the main function for Forward and Backward
double DeepCNF_Seq::ComputeScore(
	int state_num,int leftState, int currState, int pos, //-> input data
	double *state_weight, FeatScore *header_node)           //-> input parameter
{
	double state_transit_score=state_weight[state_num+leftState*state_num+currState];  //-> if leftState is -1, it will output dummy weight for state
	double header_node_score=header_node[pos*state_num+currState];
	double score=state_transit_score+header_node_score;
	return m_pModel->label_weight[currState]*score;  //-> consider the label weight, by default every label weight is 1.0 //__2015_04_15__//
}

//-> calculate partition function
double DeepCNF_Seq::CalcPartition(int state_num, int sequence_length, ScoreMatrix* forward )
{
	double Partition = (double)LogScore_ZERO;
	for(int k=0;k<state_num;k++)
	{
		LogScore_PLUS_EQUALS(Partition, (*forward)(k,sequence_length-1));
	}
	return Partition;
}

//------- get state -----//
int DeepCNF_Seq::GetState(int pos, int sequence_length, int *state_label)
{
	if(pos<0 || pos>=sequence_length) return -1;  //-> for dummy state
	return state_label[pos];
}



//========================== MAP and Viterbi function =========================//


//--------- compute viterbi -----------//
void DeepCNF_Seq::ComputeViterbi(int state_num, int sequence_length, 
		double *state_weight, FeatScore *header_node, int* predicted_label)
{
	//--> Viterbi Matrix
	ScoreMatrix best(state_num,sequence_length);
	best.Fill((double)LogScore_ZERO);
	//--> TraceBack Matrix
	ScoreMatrix traceback(state_num,sequence_length);
	traceback.Fill(-1);  //-> -1 for dummy
	//--> compute the scores for the first position
	for(int i=0;i<state_num;i++)
	{
		best(i,0)=ComputeScore(state_num,-1,i,0,state_weight,header_node);
	}
	//--> Compute the Viterbi Matrix
	for(int t=1;t<sequence_length;t++)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			for(int leftState=0;leftState<state_num;leftState++)
			{
				double new_score = ComputeScore(state_num,leftState,currState,t,state_weight,header_node) + best(leftState,t-1);
				if(new_score > best(currState,t))
				{
					best(currState,t) = new_score;
					traceback(currState,t) = leftState;
				}
			}
		}
	}
	//--> get maximal
	double max_s = LogScore_ZERO;
	int last_state = 0;
	for(int i=0;i<state_num;i++)
	{
		if(best(i,sequence_length-1)>max_s)max_s = best(i,sequence_length-1), last_state = i;
	}
	//--> TraceBack
	for(int t=sequence_length-1; t>=0;t--)
	{
		predicted_label[t]=last_state;
		last_state=(int)traceback(last_state,t);
	}
}

//-------- MAP related -------//
//-> compute forward: the main function for probability 
void DeepCNF_Seq::ComputeForward(int state_num,int sequence_length, 
	double *state_weight, FeatScore *header_node, ScoreMatrix* forward)
{
	//--> init
	forward->Fill(LogScore_ZERO);
	for(int i=0;i<state_num;i++)
	{
		(*forward)(i,0)=ComputeScore(state_num,-1,i,0,state_weight,header_node);
	}
	//--> fill-up
	for(int t=1;t<sequence_length;t++)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			for(int leftState=0;leftState<state_num;leftState++) 
			{
				double new_score = ComputeScore(state_num,leftState,currState,t,state_weight,header_node);
				LogScore_PLUS_EQUALS((*forward)(currState,t),new_score + (*forward)(leftState,t-1));
			}
		}
	}
}

//-> compute backward: the main function for probability 
void DeepCNF_Seq::ComputeBackward(int state_num,int sequence_length, 
	double *state_weight, FeatScore *header_node, ScoreMatrix* backward)
{
	//--> init
	backward->Fill(LogScore_ZERO);
	for(int i=0;i<state_num;i++)
	{
		(*backward)(i,sequence_length-1)=0;
	}
	//--> fill-up
	for(int t=sequence_length-2;t>=0;t--)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			for(int rightState=0;rightState<state_num;rightState++)
			{
				double new_score = ComputeScore(state_num,currState,rightState,t+1,state_weight,header_node);
				LogScore_PLUS_EQUALS((*backward)(currState,t),new_score + (*backward)(rightState,t+1));
			}
		}
	}
}

//-> compute MAP assignment
void DeepCNF_Seq::ComputeMAP(int state_num, int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, int* predicted_label)
{
	//--> TraceBack
	for(int t=0;t<sequence_length;t++)
	{
		int idx = 0;
		double maxS = LogScore_ZERO;
		for(int i=0;i<state_num;i++)
		{
			double s = 0;
			s = (*backward)(i,t) + (*forward)(i,t);
			if(s > maxS) maxS = s,idx = i;
		}
		predicted_label[t]=idx;
	}
}

//-> compute MAP probability
void DeepCNF_Seq::ComputeMAP(int state_num, int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, vector <vector <double> > &output)
{
	//--> calculate partition
	double Partition = CalcPartition(state_num,sequence_length,forward);
	//--> TraceBack
	output.clear();
	for(int t=0;t<sequence_length;t++)
	{
		//init
		vector <double> output_tmp(state_num,0);
		//calculate
		double tot_post=0;
		for(int i=0;i<state_num;i++)
		{
			double post=exp( (*forward)(i,t) +	(*backward)(i,t) - Partition );
			//assign
			if(post>1)post=1;
			if(post<0)post=0;
			output_tmp[i]=post;
			tot_post+=post;
		}
		//normalize
		if(tot_post==0)for(int i=0;i<state_num;i++)output_tmp[i]=1.0/state_num;
		else for(int i=0;i<state_num;i++)output_tmp[i]=1.0*output_tmp[i]/tot_post;
		//push_back
		output.push_back(output_tmp);
	}
}

// ====================== compute_obj and compute_gradient ======================//

//-> calculate objective function ( the only requirement is 'ScoreMatrix* forward', which should be calculated before. )
double DeepCNF_Seq::ComputeObj(int state_num,int sequence_length, 
	ScoreMatrix* forward, double *state_weight, FeatScore *header_node, int *state_label)
{
	//--> calculate partition
	double obj = -1.0*CalcPartition(state_num,sequence_length,forward);
	//--> calculate objective function
	for(int t=0;t<sequence_length;t++)
	{
		int leftState = GetState(t-1,sequence_length,state_label);
		int currState = GetState(t,sequence_length,state_label);
		obj+=ComputeScore(state_num,leftState,currState,t,state_weight,header_node);
	}
	//--> return objective function
	return obj;
}

//-> calculate gradient function ( the only requirement is 'ScoreMatrix* forward, backward', and the results are 'double *state_error, state_grad )
void DeepCNF_Seq::ComputeGrad(int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, double *state_weight, FeatScore *header_node, int *state_label,
	double *state_grad, FeatScore *state_error)
{
	//--> init gradient
	for(int i=0;i<(state_num+1)*state_num;i++)state_grad[i]=0;
	for(int i=0;i<sequence_length*state_num;i++)state_error[i]=0;
	//--> calculate gradient function
	double Partition = CalcPartition(state_num,sequence_length,forward);
	for(int t=0;t<sequence_length;t++)
	{
		//---> calculate observed error
		int leftState_ = GetState(t-1,sequence_length,state_label);
		int currState_ = GetState(t,sequence_length,state_label);
		state_grad[state_num + leftState_*state_num + currState_] += 1;  //--> for state_grad
		state_error[t*state_num + currState_] += 1;                      //--> for state_error
		//---> calculate expected error
		for(int currState=0;currState<state_num;currState++)
		{
			if(t==0)  // t==0
			{
				int leftState = -1;
				//calculate prob
				double prob = (*forward)(currState,t) + (*backward)(currState,t) - Partition;
				prob = exp(prob);
				//-> update state_grad
				state_grad[state_num + leftState*state_num + currState] -= prob;
				//-> update state_error
				state_error[t*state_num + currState] -= prob;
			}
			else      // t!=0
			{
				//--> for state_grad
				for(int leftState=0;leftState<state_num;leftState++)
				{
					//current score
					double curr_sco = ComputeScore(state_num,leftState,currState,t,state_weight,header_node);
					//calculate prob
					double prob = (*forward)(leftState,t-1) + curr_sco + (*backward)(currState,t) - Partition;
					prob = exp(prob);
					//-> update state_grad
					state_grad[state_num + leftState*state_num + currState] -= prob;
				}
				//--> for state_error
				{
					//calculate prob
					double prob = (*forward)(currState,t) + (*backward)(currState,t) - Partition;
					prob = exp(prob);
					//-> update state_error
					state_error[t*state_num + currState] -= prob;
				}
			}
		}
	}
}


//====================== maximize posterior probability ======================//
/*
cite the following paper:

"Training Conditional Random Fields for Maximum Labelwise Accuracy" -> 2007
by Samuel S. Gross, Olga Russakovsky, Chuong B. Do, and Serafim Batzoglou

*/

//-> ComputeForward2
void DeepCNF_Seq::ComputeForward2(int state_num,int sequence_length, double *state_weight, FeatScore *header_node, 
	ScoreMatrix* forward, ScoreMatrix* forward2, vector <double> &q_function, vector <int> &chosen_label)
{
	//--> init
	forward2->Fill(LogScore_ZERO);
	for(int i=0;i<state_num;i++)
	{
		if(i==chosen_label[0])
		{
			(*forward2)(i,0)=q_function[0]+(*forward)(i,0);
		}
	}
	//--> fill-up
	for(int t=1;t<sequence_length;t++)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			for(int leftState=0;leftState<state_num;leftState++) 
			{
				double new_score = ComputeScore(state_num,leftState,currState,t,state_weight,header_node);
				LogScore_PLUS_EQUALS((*forward2)(currState,t),new_score + (*forward2)(leftState,t-1));
				//-> for each observed label
				if(currState==chosen_label[t])
				{
					LogScore_PLUS_EQUALS((*forward2)(currState,t),new_score + q_function[t]+(*forward)(leftState,t-1));
				}
			}
		}
	}
}
//-> ComputeBackward2
void DeepCNF_Seq::ComputeBackward2(int state_num,int sequence_length, double *state_weight, FeatScore *header_node, 
	ScoreMatrix* backward, ScoreMatrix* backward2, vector <double> &q_function, vector <int> &chosen_label)
{
	//--> init
	backward2->Fill(LogScore_ZERO);
	//--> fill-up
	for(int t=sequence_length-2;t>=0;t--)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			for(int rightState=0;rightState<state_num;rightState++)
			{
				double new_score = ComputeScore(state_num,currState,rightState,t+1,state_weight,header_node);
				LogScore_PLUS_EQUALS((*backward2)(currState,t),new_score + (*backward2)(rightState,t+1));
				//-> for each observed label
				if(rightState==chosen_label[t+1])
				{
					LogScore_PLUS_EQUALS((*backward2)(currState,t),new_score + q_function[t+1]+(*backward)(rightState,t+1));
				}
			}
		}
	}
}

//-------------- for compute_function and compute_gradient of maximal accuracy -----------//
//-> ComputeFunction2
double DeepCNF_Seq::ComputeObj2(int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, 
	double *state_weight, FeatScore *header_node, int *state_label)
{
	//--> calculate partition
	double Partition = CalcPartition(state_num,sequence_length,forward);
	//--> calculate objective function
	double func = 0;
	for(int t=0;t<sequence_length;t++)
	{
		for(int i=0;i<state_num;i++)
		{
			if(i==state_label[t])
			{
				double prob = exp((*forward)(i,t)+(*backward)(i,t)-Partition);
				func+=prob;
			}
		}
	}
	return func;
}
//-> ComputeGradient2
void DeepCNF_Seq::ComputeGrad2(int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
	double *state_weight, FeatScore *header_node, int *state_label, double *state_grad, FeatScore *state_error)
{
	//--> init gradient
	for(int i=0;i<(state_num+1)*state_num;i++)state_grad[i]=0;
	for(int i=0;i<sequence_length*state_num;i++)state_error[i]=0;
	//--> compute forward2 and backward2
	//init q_function and chosen_label
	vector <double> q_function(sequence_length,0); //-> 0 means log(1.0)
	vector <int> chosen_label(sequence_length,0);
	for(int i=0;i<sequence_length;i++)chosen_label[i]=state_label[i];
	//calculate
	ComputeForward2(dim_states,sequence_length,
		state_weight,header_node,forward,forward2,q_function,chosen_label);
	ComputeBackward2(dim_states,sequence_length,
		state_weight,header_node,backward,backward2,q_function,chosen_label);
	//--> calculate gradient function
	double Partition = CalcPartition(state_num,sequence_length,forward);
	double func = ComputeObj2(state_num,sequence_length,forward,backward,state_weight,header_node,state_label);
	//--> calculate gradient for maximal accuracy
	for(int t=0;t<sequence_length;t++)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			if(t==0)  // t==0
			{
				int leftState = -1;
				//-> prepare probability
				{
					//prob0
					double prob0 = (*forward)(currState, t) + (*backward)(currState, t) - Partition;
					prob0 = exp(prob0);
					//prob1
					double prob1 = (*forward2)(currState, t) + (*backward)(currState, t) - Partition;
					prob1 = exp(prob1);
					//prob2
					double prob2 = (*forward)(currState, t) + (*backward2)(currState, t) - Partition;
					prob2 = exp(prob2);
					//prob
					double prob = prob0*func - prob1 - prob2;
					//-> update state_grad
					state_grad[state_num + leftState*state_num + currState] -= prob;
					//-> update state_error
					state_error[t*state_num + currState] -= prob;
				}
			}
			else      // t!=0
			{
				//--> for state_grad
				for(int leftState=0;leftState<state_num;leftState++)
				{
					//-> prepare probability
					//current score
					double curr_sco = ComputeScore(state_num,leftState,currState,t,state_weight,header_node);
					//prob0
					double prob0 = (*forward)(leftState, t-1) + curr_sco + (*backward)(currState, t) - Partition;
					prob0 = exp(prob0);
					//prob1
					double prob1 = (*forward2)(leftState, t-1) + curr_sco + (*backward)(currState, t) - Partition;
					prob1 = exp(prob1);
					//prob2
					double prob2 = (*forward)(leftState, t-1) + curr_sco + (*backward2)(currState, t) - Partition;
					prob2 = exp(prob2);
					//prob
					double prob = prob0*func - prob1 - prob2;
					//-> for each observed label
					if(currState==state_label[t])prob -= prob0;
					//-> update state_grad
					state_grad[state_num + leftState*state_num + currState] -= prob;
				}
				//--> for state_error
				{
					//-> prepare probability
					//prob0
					double prob0 = (*forward)(currState, t) + (*backward)(currState, t) - Partition;
					prob0 = exp(prob0);
					//prob1
					double prob1 = (*forward2)(currState, t) + (*backward)(currState, t) - Partition;
					prob1 = exp(prob1);
					//prob2
					double prob2 = (*forward)(currState, t) + (*backward2)(currState, t) - Partition;
					prob2 = exp(prob2);
					//prob
					double prob = prob0*func - prob1 - prob2;
					//-> update state_grad
					state_error[t*state_num + currState] -= prob;
				}
			}
		}
	}
}

//================== Max AUC ===============//__150517__//
//-> basic ideas
/*

1) for a given label t, optimize the following approximated AUC

sigma_miu, sigma_l, gamma_[miu,l] * S( P^l, D_t) * V( P^(miu-l), D_!t)

where S function is,

S( P, D_t ) = sigma_i, delta(i,t) P(y_i,t)

where P(y_i,t) is the posterior probability at position i with label t,
delta(i,t) is the indicator function that if the real label of the i-th position equals to ¦Ó, then the value is 1, and otherwise 0.

*/

//-> Compute q_function
double DeepCNF_Seq::Compute_q_function(
	vector <vector <double> > &post_prob, 
	vector <double> &q_function,double &tot_num,
	int power_value, int *state_label,int chosen_label,
	int CHOOSEorNOT,int DERIVEorNOT)
{
	double total_func=0;
	tot_num=0;
	//for each position
	q_function.resize(post_prob.size(),LogScore_ZERO);
	for(int t=0; t<(int)post_prob.size();t++)
	{
		double value=0;
		//calculate the value
		if(DERIVEorNOT==1) // derivative the function
		{
			if(power_value==0)value=0;
			else value=power_value*pow(post_prob[t][chosen_label],1.0*(power_value-1));
		}
		else               // original value
		{
			if(power_value==0)value=1;
			else value=pow(post_prob[t][chosen_label],1.0*power_value);
		}
		//choose the label or not
		if(CHOOSEorNOT==1) // choose the label, for s_value
		{
			if(state_label[t]==chosen_label)
			{
				if(value<exp(LogScore_ZERO))q_function[t]=LogScore_ZERO;
				else q_function[t]=log(value);
				if(DERIVEorNOT==1) total_func+=value*post_prob[t][chosen_label];
				else total_func+=value;
				tot_num++;
			}
			else q_function[t]=LogScore_ZERO;
		}
		else               // not choose the label,for v_value
		{
			if(state_label[t]==chosen_label)q_function[t]=LogScore_ZERO;
			else
			{
				if(value<exp(LogScore_ZERO))q_function[t]=LogScore_ZERO;
				else q_function[t]=log(value);
				if(DERIVEorNOT==1) total_func+=value*post_prob[t][chosen_label];
				else total_func+=value;
				tot_num++;
			}
		}
	}
	//return
	return total_func;
}

//-------------- for compute_function and compute_gradient of maximal AUC -----------//
//------------- Compute Function ------------//
//-> ComputeFunction3
//[note]: if GRADorNOT==1, then perform obj_for_grad; otherwisse, perform original obj
void DeepCNF_Seq::ComputeObj3( 
	vector <vector <double> > &post_prob,
	int *state_label,int chosen_label,
	int power_value_l,int power_value_u,
	vector <double> &q_function_s, vector <double> &q_function_v,
	double &s_value, double &v_value, 
	double &s_value_n,double &v_value_n,
	int GRADorNOT)
{
	//--> calculate s_value
	s_value=Compute_q_function(post_prob,q_function_s,s_value_n,
		power_value_l,state_label,chosen_label,1,GRADorNOT);
	//--> calculate v_value
	v_value=Compute_q_function(post_prob,q_function_v,v_value_n,
		power_value_u-power_value_l,state_label,chosen_label,0,GRADorNOT);
}
//-> compute objective total
void DeepCNF_Seq::ComputeObj3_Total(
	int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, 
	double *state_weight, FeatScore *header_node, int *state_label,
	int max_degree,
	double *s_value_out,double *v_value_out,
	double *s_value_out_n,double *v_value_out_n)
{
	//-> init value
	int total_num=(max_degree + 2)*(max_degree + 1)/2;
	int totnum=total_num*dim_states;
	for(int i=0;i<totnum;i++)s_value_out[i]=0;
	for(int i=0;i<totnum;i++)v_value_out[i]=0;
	for(int i=0;i<totnum;i++)s_value_out_n[i]=0;
	for(int i=0;i<totnum;i++)v_value_out_n[i]=0;
	//-> calculate posterior probability
	vector <vector <double> > post_prob;
	ComputeMAP(state_num, sequence_length, forward, backward, post_prob);
	//-> for each label
	for(int k=0;k<state_num;k++)
	{
		//--> for each degree
		int idx=0;
		for(int i = 0; i <= max_degree; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				//calculate object function given label and degree
				vector <double> q_function_s;
				vector <double> q_function_v;
				double s_value;
				double v_value;
				double s_value_n;
				double v_value_n;
				ComputeObj3(post_prob,state_label,k,j,i,q_function_s,q_function_v,
					s_value,v_value,s_value_n,v_value_n,0);
				//idx++
				s_value_out[idx*state_num+k]+=s_value;
				v_value_out[idx*state_num+k]+=v_value;
				s_value_out_n[idx*state_num+k]+=s_value_n;
				v_value_out_n[idx*state_num+k]+=v_value_n;
				idx++;
			}
		}
	}
}


//------------- Compute Gradient ------------//
//-> ComputeGradient3 part
void DeepCNF_Seq::ComputeGrad3_Part(int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
	double *state_weight, FeatScore *header_node, int chosen_label,
	vector <double> &state_grad, vector <double> &state_error,
	vector <double> &q_function, double func, double Partition)
{
	for(int t=0;t<sequence_length;t++)
	{
		for(int currState=0;currState<state_num;currState++)
		{
			if(t==0)  // t==0
			{
				int leftState = -1;
				//-> prepare probability
				{
					//prob0
					double prob0 = (*forward)(currState, t) + (*backward)(currState, t) - Partition;
					prob0 = exp(prob0);
					//prob1
					double prob1 = (*forward2)(currState, t) + (*backward)(currState, t) - Partition;
					prob1 = exp(prob1);
					//prob2
					double prob2 = (*forward)(currState, t) + (*backward2)(currState, t) - Partition;
					prob2 = exp(prob2);
					//prob
					double prob = prob0*func - prob1 - prob2;
					//-> update state_grad
					state_grad[state_num + leftState*state_num + currState] -= prob;
					//-> update state_error
					state_error[t*state_num + currState] -= prob;
				}
			}
			else      // t!=0
			{
				//--> for state_grad
				for(int leftState=0;leftState<state_num;leftState++)
				{
					//-> prepare probability
					//current score
					double curr_sco = ComputeScore(state_num,leftState,currState,t,state_weight,header_node);
					//prob0
					double prob0 = (*forward)(leftState, t-1) + curr_sco + (*backward)(currState, t) - Partition;
					prob0 = exp(prob0);
					//prob1
					double prob1 = (*forward2)(leftState, t-1) + curr_sco + (*backward)(currState, t) - Partition;
					prob1 = exp(prob1);
					//prob2
					double prob2 = (*forward)(leftState, t-1) + curr_sco + (*backward2)(currState, t) - Partition;
					prob2 = exp(prob2);
					//prob
					double prob = prob0*func - prob1 - prob2;
					//-> for each observed label
					if(currState==chosen_label)prob -= exp(q_function[t])*prob0;
					//-> update state_grad
					state_grad[state_num + leftState*state_num + currState] -= prob;
				}
				//--> for state_error
				{
					//-> prepare probability
					//prob0
					double prob0 = (*forward)(currState, t) + (*backward)(currState, t) - Partition;
					prob0 = exp(prob0);
					//prob1
					double prob1 = (*forward2)(currState, t) + (*backward)(currState, t) - Partition;
					prob1 = exp(prob1);
					//prob2
					double prob2 = (*forward)(currState, t) + (*backward2)(currState, t) - Partition;
					prob2 = exp(prob2);
					//prob
					double prob = prob0*func - prob1 - prob2;
					//-> update state_grad
					state_error[t*state_num + currState] -= prob;
				}
			}
		}
	}
}
//-> ComputeGradient3
void DeepCNF_Seq::ComputeGrad3(int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
	double *state_weight, FeatScore *header_node, int *state_label, double *state_grad, FeatScore *state_error,
	int chosen_label,int power_value_l,int power_value_u,double gamma_value,
	double s_value_sum, double v_value_sum,double s_value_sum_n, double v_value_sum_n )
{
	//--> init gradient
	vector <double> state_grad_s ((state_num+1)*state_num,0);
	vector <double> state_error_s (sequence_length*state_num,0);
	vector <double> state_grad_v ((state_num+1)*state_num,0);
	vector <double> state_error_v (sequence_length*state_num,0);

	//--> calculate q_function and parition_function
	//-> 1. calculate posterior probability
	vector <vector <double> > post_prob;
	ComputeMAP(state_num, sequence_length, forward, backward, post_prob);
	//-> 2. calculate q_function and func_value
	vector <double> q_function_s;
	vector <double> q_function_v;
	double s_value,v_value;
	double s_value_n,v_value_n;
	ComputeObj3(post_prob,state_label,chosen_label,
		power_value_l,power_value_u,q_function_s,q_function_v,s_value,v_value,s_value_n,v_value_n,1);
	//-> 3. calcualte partition funcrtion
	double Partition = CalcPartition(state_num,sequence_length,forward);

	//--> calculate gradient for maximal AUC
	//part0. chosen_label
	vector <int> chosen_label_vec(sequence_length,chosen_label);
	//part1. s_value
	//-> 1.1 calculate forward2 and backward2
	ComputeForward2(dim_states,sequence_length,
		state_weight,header_node,forward,forward2,q_function_s,chosen_label_vec);
	ComputeBackward2(dim_states,sequence_length,
		state_weight,header_node,backward,backward2,q_function_s,chosen_label_vec);
	//-> 1.2 calculate gradient
	ComputeGrad3_Part(state_num,sequence_length,forward,backward,forward2,backward2,
		state_weight,header_node,chosen_label,state_grad_s,state_error_s,
		q_function_s,s_value,Partition);
	//part2. v_value
	//-> 2.1 calculate forward2 and backward2
	ComputeForward2(dim_states,sequence_length,
		state_weight,header_node,forward,forward2,q_function_v,chosen_label_vec);
	ComputeBackward2(dim_states,sequence_length,
		state_weight,header_node,backward,backward2,q_function_v,chosen_label_vec);
	//-> 2.2 calculate gradient
	ComputeGrad3_Part(state_num,sequence_length,forward,backward,forward2,backward2,
		state_weight,header_node,chosen_label,state_grad_v,state_error_v,
		q_function_v,v_value,Partition);

	//---> calculate the final state_grad and state_error
	for(int i=0;i<(state_num+1)*state_num;i++)
	{
		state_grad[i]+=gamma_value*( state_grad_s[i]*v_value_sum + state_grad_v[i]*s_value_sum )/s_value_sum_n/v_value_sum_n;
	}
	for(int i=0;i<sequence_length*state_num;i++)
	{
		state_error[i]+=gamma_value*( state_error_s[i]*v_value_sum + state_error_v[i]*s_value_sum )/s_value_sum_n/v_value_sum_n;
	}
}
//-> compute gradient total
void DeepCNF_Seq::ComputeGrad3_Total(int state_num,int sequence_length, 
	ScoreMatrix* forward, ScoreMatrix* backward, ScoreMatrix* forward2, ScoreMatrix* backward2, 
	double *state_weight, FeatScore *header_node, int *state_label,
	double *state_grad, FeatScore *state_error,
	int max_degree, double *gamma_coeff,
	double *s_value_sum, double *v_value_sum,
	double *s_value_sum_n, double *v_value_sum_n )
{
	//-> init gradient
	for(int i=0;i<(state_num+1)*state_num;i++)state_grad[i]=0;
	for(int i=0;i<sequence_length*state_num;i++)state_error[i]=0;
	//-> for each label
	for(int k=0;k<state_num;k++)
	{
		//--> for each gamma_coefficient
		int idx=0;
		for(int i = 0; i <= max_degree; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				//get gamma_coefficient
				double gamma_value=gamma_coeff[idx];
				//calculate gradient and error_function
				ComputeGrad3(state_num,sequence_length,forward,backward,forward2,backward2,
					state_weight,header_node,state_label,state_grad,state_error,k,j,i,gamma_value,
					s_value_sum[idx*state_num+k],v_value_sum[idx*state_num+k],
					s_value_sum_n[idx*state_num+k],v_value_sum_n[idx*state_num+k]);
				//idx++
				idx++;
			}
		}
	}
}


//------ calculate test accuracy --------//
void DeepCNF_Seq::ComputeTestAccuracy(
	int sequence_length,int *observed_label, int *predicted_label,
	double &total_num,double &correct_num)
{
	total_num += sequence_length;
	for(int t=0; t<sequence_length;t++)
		if(observed_label[t]==predicted_label[t])correct_num++;
}
void DeepCNF_Seq::ComputeTestAccuracy(
	double &total_num,double &correct_num)
{
	ComputeTestAccuracy(sequence_length,observed_label,predicted_label,total_num,correct_num);
}
void DeepCNF_Seq::ComputeTestAccuracy_Weight(
	int sequence_length,int *observed_label, int *predicted_label,
	double &total_num,double &correct_num)
{
	for(int t=0; t<sequence_length;t++)
	{
		if(observed_label[t]==predicted_label[t])
		{
			correct_num+=m_pModel->label_weight[observed_label[t]];
		}
		total_num+=m_pModel->label_weight[observed_label[t]];
	}
}
void DeepCNF_Seq::ComputeTestAccuracy_Weight(
	double &total_num,double &correct_num)
{
	ComputeTestAccuracy_Weight(sequence_length,observed_label,predicted_label,total_num,correct_num);
}


//=============== main functions for prediction and training ==============//

//---- read in features --//
void DeepCNF_Seq::Read_Init_Features(vector <string> &input)
{
	for(int j=0;j<sequence_length;j++)
	{
		istringstream trn_in(input[j]);
		for(int k=0;k<dim_one_pos;k++)
		{
			FeatScore value;
			trn_in >> value;
			feat_orig[0][j*dim_one_pos+k]=value;
		}
	}
}
void DeepCNF_Seq::Read_Init_Features(vector <vector <FeatScore> > &input)
{
	for(int j=0;j<sequence_length;j++)
	{
		for(int k=0;k<dim_one_pos;k++)
		{
			double value=input[j][k];
			feat_orig[0][j*dim_one_pos+k]=value;
		}
	}
}

//---- calculate forward/backward ---//
void DeepCNF_Seq::Calc_Forward_Backward(void)
{
	//--> calculate
	FeedForward(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
	ComputeForward(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],forward);
	ComputeBackward(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],backward);
}

//---- for prediction ----//
//-> gate output
void DeepCNF_Seq::Gate_Output(vector <vector <FeatScore> > &input,vector <vector <FeatScore> > &output)
{
	//--> read init features
	Read_Init_Features(input);
	//--> feed-forward
	FeedForward(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
	//--> output lables
	output.clear();
	int cur_layer=m_pModel->tot_layer-2;
	for(int t=0;t<sequence_length;t++)
	{
		vector <FeatScore> tmp_rec;
		U_INT pivot=m_pModel->nodes_dim[cur_layer]*t;
		for(int k=0;k<m_pModel->nodes_dim[cur_layer];k++)
		{
			FeatScore value=feat_orig[cur_layer][pivot+k];
			tmp_rec.push_back(value);
		}
		output.push_back(tmp_rec);
	}
}
void DeepCNF_Seq::Gate_Output(vector <string> &input,vector <vector <FeatScore> > &output)
{
	//--> read init features
	Read_Init_Features(input);
	//--> feed-forward
	FeedForward(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
	//--> output lables
	output.clear();
	int cur_layer=m_pModel->tot_layer-2;
	for(int t=0;t<sequence_length;t++)
	{
		vector <FeatScore> tmp_rec;
		U_INT pivot=m_pModel->nodes_dim[cur_layer]*t;
		for(int k=0;k<m_pModel->nodes_dim[cur_layer];k++)
		{
			FeatScore value=feat_orig[cur_layer][pivot+k];
			tmp_rec.push_back(value);
		}
		output.push_back(tmp_rec);
	}
}
void DeepCNF_Seq::Gate_Output(vector <vector <FeatScore> > &output)
{
	//--> feed-forward
	FeedForward(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
	//--> output lables
	output.clear();
	int cur_layer=m_pModel->tot_layer-2;
	for(int t=0;t<sequence_length;t++)
	{
		vector <FeatScore> tmp_rec;
		U_INT pivot=m_pModel->nodes_dim[cur_layer]*t;
		for(int k=0;k<m_pModel->nodes_dim[cur_layer];k++)
		{
			FeatScore value=feat_orig[cur_layer][pivot+k];
			tmp_rec.push_back(value);
		}
		output.push_back(tmp_rec);
	}
}


//-> MAP probability
void DeepCNF_Seq::MAP_Probability(vector <vector <FeatScore> > &input,vector <vector <double> > &output)
{
	//--> read init features
	Read_Init_Features(input);
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,output);
}
void DeepCNF_Seq::MAP_Probability(vector <string> &input,vector <vector <double> > &output)
{
	//--> read init features
	Read_Init_Features(input);
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,output);
}
void DeepCNF_Seq::MAP_Probability(vector <vector <double> > &output)
{
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,output);
}

//-> MAP assignment
void DeepCNF_Seq::MAP_Assign(vector <vector <FeatScore> > &input,vector <int> &output_lab)
{
	//--> read init features
	Read_Init_Features(input);
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,predicted_label);
	output_lab.clear();
	for(int i=0;i<sequence_length;i++)output_lab.push_back(predicted_label[i]);
}
void DeepCNF_Seq::MAP_Assign(vector <string> &input,vector <int> &output_lab)
{
	//--> read init features
	Read_Init_Features(input);
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,predicted_label);
	output_lab.clear();
	for(int i=0;i<sequence_length;i++)output_lab.push_back(predicted_label[i]);
}
void DeepCNF_Seq::MAP_Assign(vector <int> &output_lab)
{
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,predicted_label);
	output_lab.clear();
	for(int i=0;i<sequence_length;i++)output_lab.push_back(predicted_label[i]);
}
void DeepCNF_Seq::MAP_Assign(int *output_lab)
{
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,predicted_label);
	for(int i=0;i<sequence_length;i++)output_lab[i]=predicted_label[i];
}
void DeepCNF_Seq::MAP_Assign(void)
{
	//--> calculate forward/backward
	Calc_Forward_Backward();
	//--> output lables
	ComputeMAP(dim_states,sequence_length,forward,backward,predicted_label);
}

//-> Viterbi assignment
void DeepCNF_Seq::Viterbi_Assign(vector <vector <FeatScore> > &input,vector <int> &output_lab)
{
	//--> read init features
	Read_Init_Features(input);
	//--> output lables
	ComputeViterbi(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],predicted_label);
	output_lab.clear();
	for(int i=0;i<sequence_length;i++)output_lab.push_back(predicted_label[i]);
}
void DeepCNF_Seq::Viterbi_Assign(vector <string> &input,vector <int> &output_lab)
{
	//--> read init features
	Read_Init_Features(input);
	//--> output lables
	ComputeViterbi(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],predicted_label);
	output_lab.clear();
	for(int i=0;i<sequence_length;i++)output_lab.push_back(predicted_label[i]);
}
void DeepCNF_Seq::Viterbi_Assign(vector <int> &output_lab)
{
	//--> output lables
	ComputeViterbi(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],predicted_label);
	output_lab.clear();
	for(int i=0;i<sequence_length;i++)output_lab.push_back(predicted_label[i]);
}
void DeepCNF_Seq::Viterbi_Assign(int *output_lab)
{
	//--> output lables
	ComputeViterbi(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],predicted_label);
	for(int i=0;i<sequence_length;i++)output_lab[i]=predicted_label[i];
}
void DeepCNF_Seq::Viterbi_Assign(void)
{
	//--> output lables
	ComputeViterbi(dim_states,sequence_length,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],predicted_label);
}

//=============== for training purpose =================//
//[note]: before calculate Obj_Calc and Grad_Calc, forward and/or backward should be prepared
//----- maximize log_probability ------//
double DeepCNF_Seq::Obj_Calc(void)
{
	return ComputeObj(dim_states,sequence_length,forward,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],observed_label);
}
void DeepCNF_Seq::Grad_Calc(void)
{
	//-> calculate overall gradient
	ComputeGrad(dim_states,sequence_length,forward,backward,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],observed_label,
		m_pModel->layer_grad[tot_layer-1],feat_orig[tot_layer-1]);
	//-> propagate back
	BackPropagate(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,m_pModel->layer_grad,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
}

//----- maximal accuracy -----//
double DeepCNF_Seq::Obj_Calc2(void)
{
	return ComputeObj2(dim_states,sequence_length,forward,backward,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],observed_label);
}
void DeepCNF_Seq::Grad_Calc2(void)
{
	//-> calculate overall gradient
	ComputeGrad2(dim_states,sequence_length,forward,backward,forward2,backward2,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],observed_label,
		m_pModel->layer_grad[tot_layer-1],feat_orig[tot_layer-1]);
	//-> propagate back
	BackPropagate(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,m_pModel->layer_grad,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
}

//----- maximize AUC ------//
void DeepCNF_Seq::Obj_Calc3(void)
{
	ComputeObj3_Total(dim_states,sequence_length,forward,backward,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],observed_label,
		m_pModel->max_degree,
		m_pModel->s_value_out,m_pModel->v_value_out,
		m_pModel->s_value_out_n,m_pModel->v_value_out_n);
}
//-> note that before computing gradient of maximize AUC, 
//        s/v_value_out and s/v_value_out_n must be pre-calculated 
//        from all samples for each label.
void DeepCNF_Seq::Grad_Calc3(void)
{
	//-> calculate overall gradient
	ComputeGrad3_Total(dim_states,sequence_length,forward,backward,forward2,backward2,
		m_pModel->layer_weight[tot_layer-1],feat_rec[tot_layer-1],observed_label,
		m_pModel->layer_grad[tot_layer-1],feat_orig[tot_layer-1],
		m_pModel->max_degree, m_pModel->gamma_coeff,
		m_pModel->s_value_out, m_pModel->v_value_out,
		m_pModel->s_value_out_n, m_pModel->v_value_out_n);
	//-> propagate back
	BackPropagate(feat_dim,feat_orig,feat_rec,
		m_pModel->layer_weight,m_pModel->layer_grad,sequence_length,
		m_pModel->tot_layer,m_pModel->nodes_dim,m_pModel->windows_dim);
}




//========================== Deep Convolution CNF model =========================//

//---- Deep Convolution CNF architecture ----//-> take layer=2 as example
/*
            iiiiiiiiii    -> input feature layer (input layer) [window 1] -> the bottom layer
            \/\/\/\/\/    -> Weight 1 (input_to_1)
             xxxxxxxx     -> neuron feature layer (layer 1) [window 2]
             \/\/\/\/     -> Weight 2 (1_to_2)
              xxxxxx      -> neuron feature layer (layer 2) [no window]   -> the header layer
              ||||||      -> Weight 2_to_label
              LLLLLL      -> label layer (output layer)
              ------      -> Weight label_to_label
*/

//========== constructor and destructor ============//
DeepCNF_Model::DeepCNF_Model(
	int layer_number,string &window_str, string &nodes_str, int n_states, int n_local,
	int Train_Or_Not_, int MAX_AUC_DEGREE, double MAX_AUC_BETA, int Gate_Function_)
{
	//-> gate function
	Gate_Function=Gate_Function_;

	//-> training switch
	Train_Or_Not=Train_Or_Not_;   //-> default is OFF [0]
	max_degree=MAX_AUC_DEGREE;    //-> default is OFF [0]
	beta_sigmoid=MAX_AUC_BETA;    //-> default is OFF [-1]

	//-> create weight parameter
	Weight_Structre_Create(layer_number,window_str,nodes_str,n_states,n_local);

	//-> training or not
	weights=new double[total_param_count];
	if(Train_Or_Not==1)
	{
		grad=new double[total_param_count];
		grad_sum=new double[total_param_count];
	}
	else
	{
		grad=NULL;
		grad_sum=NULL;
	}

	//-> label weight
	label_weight=new double[n_states];
	for(int i=0;i<n_states;i++)label_weight[i]=1;

	//-> max AUC or not
	if(MAX_AUC_DEGREE>0)
	{
		int total_num=(max_degree + 2)*(max_degree + 1)/2;
		int totnum=total_num*n_states;
		gamma_coeff=new double[total_num];
		s_value_out=new double[totnum];
		v_value_out=new double[totnum];
		s_value_sum=new double[totnum];
		v_value_sum=new double[totnum];
		s_value_out_n=new double[totnum];
		v_value_out_n=new double[totnum];
		s_value_sum_n=new double[totnum];
		v_value_sum_n=new double[totnum];
		auc_for_label=new double[n_states];
		//init gamma coefficient
		init_gamma_coeff(gamma_coeff,max_degree,beta_sigmoid);
	}
	else
	{
		gamma_coeff=NULL;
		s_value_out=NULL;
		v_value_out=NULL;
		s_value_sum=NULL;
		v_value_sum=NULL;
		s_value_out_n=NULL;
		v_value_out_n=NULL;
		s_value_sum_n=NULL;
		v_value_sum_n=NULL;
		auc_for_label=NULL;
	}
}
DeepCNF_Model::DeepCNF_Model(int Train_Or_Not_,
	int MAX_AUC_DEGREE,double MAX_AUC_BETA, int Gate_Function_)
{
	//-> gate function
	Gate_Function=Gate_Function_;

	//-> training switch
	Train_Or_Not=Train_Or_Not_;   //-> default is OFF [0]
	max_degree=MAX_AUC_DEGREE;    //-> default is OFF [0]
	beta_sigmoid=MAX_AUC_BETA;    //-> default is OFF [-1]
	
	//-> for actual computation
	layer_weight=NULL;
	layer_grad=NULL;
	layer_count=NULL;
	windows_dim=NULL;
	nodes_dim=NULL;

	//-> for optimization
	weights=NULL;
	grad=NULL;
	grad_sum=NULL;

	//-> for label weight
	label_weight=NULL;

	//-> for max AUC
	gamma_coeff=NULL;
	s_value_out=NULL;
	v_value_out=NULL;
	s_value_sum=NULL;
	v_value_sum=NULL;
	s_value_out_n=NULL;
	v_value_out_n=NULL;
	s_value_sum_n=NULL;
	v_value_sum_n=NULL;
	auc_for_label=NULL;
}
DeepCNF_Model::~DeepCNF_Model(void)
{
	//-> for actual computation
	if(layer_weight!=NULL)
	{
		for(int i=0;i<tot_layer;i++)delete [] layer_weight[i];
		delete [] layer_weight;
	}
	if(layer_grad!=NULL)
	{
		for(int i=0;i<tot_layer;i++)delete [] layer_grad[i];
		delete [] layer_grad;
	}
	if(layer_count!=NULL)delete [] layer_count;
	if(windows_dim!=NULL)delete [] windows_dim;
	if(nodes_dim!=NULL)delete [] nodes_dim;

	//-> for optimization
	if(weights!=NULL)delete [] weights;
	if(grad!=NULL)delete [] grad;
	if(grad_sum!=NULL)delete [] grad_sum;

	//-> for label weight
	if(label_weight!=NULL)delete [] label_weight;

	//-> for max AUC
	if(gamma_coeff!=NULL)delete [] gamma_coeff;
	if(s_value_out!=NULL)delete [] s_value_out;
	if(v_value_out!=NULL)delete [] v_value_out;
	if(s_value_sum!=NULL)delete [] s_value_sum;
	if(v_value_sum!=NULL)delete [] v_value_sum;
	if(s_value_out_n!=NULL)delete [] s_value_out_n;
	if(v_value_out_n!=NULL)delete [] v_value_out_n;
	if(s_value_sum_n!=NULL)delete [] s_value_sum_n;
	if(v_value_sum_n!=NULL)delete [] v_value_sum_n;
	if(auc_for_label!=NULL)delete [] auc_for_label;

}

//----- collect gradsum -----//
void DeepCNF_Model::Collect_GradSum(double *grad_sum)
{
	//MultiDim_To_OneDim
	MultiDim_To_OneDim(tot_layer,layer_count,layer_grad,grad);
	for(U_INT i=0;i<total_param_count;i++)grad_sum[i]+=grad[i];
}
void DeepCNF_Model::Collect_GradSum(void)
{
	//MultiDim_To_OneDim
	MultiDim_To_OneDim(tot_layer,layer_count,layer_grad,grad);
	for(U_INT i=0;i<total_param_count;i++)grad_sum[i]+=grad[i];
}

//================ weight structure create for the whole model ===============//
void DeepCNF_Model::Weight_Structre_Create(
	int layer_number,string &windows_str,string &nodes_str,int n_states, int n_local)  //-> input for layer and window
{
	// ----- string parser ---//
	vector <int> layer_windows;
	int wind_ret=Parse_Str(windows_str,layer_windows,',');
	vector <int> layer_nodes;
	int node_ret=Parse_Str(nodes_str,layer_nodes,',');
	// ----- init check ------//
	if(layer_number!=wind_ret || layer_number!=node_ret)
	{
		fprintf(stderr,"dimension not valid !! layer_number %d != wind_ret %d != node_ret %d \n",
			layer_number,wind_ret,node_ret);
		exit(-1);
	}

	// ------- dimension related ----------//
	//-> layber number
	int window_layer=layer_number;   //-> here is the layer having windows
	int total_layer=layer_number+2;  //-> we consider header layer and label layer
	//-> window dimension
	int *windows_dimension=new int[layer_number]; //-> only bottom and the N-1 layer have windows
	for(int i=0;i<layer_number;i++)windows_dimension[i]=layer_windows[i];
	//-> nodes dimension (consider the bottom layer, so the size is layer_number+1)
	int *nodes_dimension=new int[layer_number+2];
	nodes_dimension[0]=n_local;                                             //-> no bias here !!
	for(int i=1;i<=layer_number;i++)nodes_dimension[i]=layer_nodes[i-1];  //-> no bias here !!
	nodes_dimension[layer_number+1]=n_states;

	// -------- weight related ----------//
	// weight total count
	U_INT *weight_total_count=new U_INT[layer_number+2];
	//-> bottom weights (from 0 -> N-1)
	for(int i=0;i<layer_number;i++)weight_total_count[i]=( (nodes_dimension[i]+1)*(2*layer_windows[i]+1)+1 )* nodes_dimension[i+1];
	//-> gate_to_label weight
	weight_total_count[layer_number]=( nodes_dimension[layer_number] )*n_states;
	//-> label_to_label weight
	weight_total_count[layer_number+1]=n_states*n_states+n_states;  //-> the last 'num_states' is for dummy state
	// weight dimension
	//-> bottom weights
	double **layer_wise_weight=new double*[layer_number+2];
	for(int i=0;i<layer_number+2;i++)layer_wise_weight[i]=new double[weight_total_count[i]];
	// grad dimension
	//-> bottom weights
	double **layer_wise_grad=NULL;
	if(Train_Or_Not==1)
	{
		layer_wise_grad=new double*[layer_number+2];
		for(int i=0;i<layer_number+2;i++)layer_wise_grad[i]=new double[weight_total_count[i]];
	}

	//---- final assign ----//
	//-> basic dimension
	dim_one_pos=n_local;
	dim_states=n_states;
	total_param_count=0;
	for(int i=0;i<layer_number+2;i++)total_param_count+=weight_total_count[i];
	//-> dimension out
	wind_layer=window_layer;          //-> layer_number
	windows_dim=windows_dimension;    //-> layer_number
	tot_layer=total_layer;            //-> layer_number+2
	nodes_dim=nodes_dimension;        //-> layer_number+2
	//-> weight out
	layer_count=weight_total_count;   //-> layer_number+2   -> this data_structure is critical for Multi<->One dimension transform
	layer_weight=layer_wise_weight;   //-> layer_number+2
	layer_grad=layer_wise_grad;       //-> layer_number+2
}

//=============== parameter re-organization =============//
void DeepCNF_Model::MultiDim_To_OneDim(int tot_layer,U_INT * layer_count,double ** layer_weight,double *out_weight)
{
	U_INT pivot=0;
	for(int i=0;i<tot_layer;i++)
	{
		for(U_INT j=0;j<layer_count[i];j++)
		{
			out_weight[pivot]=layer_weight[i][j];
			pivot++;
		}
	}
}
void DeepCNF_Model::OneDim_To_MultiDim(int tot_layer,U_INT * layer_count,double ** layer_weight,double *in_weight)
{
	U_INT pivot=0;
	for(int i=0;i<tot_layer;i++)
	{
		for(U_INT j=0;j<layer_count[i];j++)
		{
			layer_weight[i][j]=in_weight[pivot];
			pivot++;
		}
	}
}

//==================== model_load ==================//
//-> from file
void DeepCNF_Model::cnf_model_load(string &file)
{
	ifstream fin;
	string buf,temp;
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",file.c_str());
		exit(-1);
	}
	//get 1st line
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"%s file format bad!\n",file.c_str());
		exit(-1);
	}
	//get feature
	U_INT totnum=0;
	vector <double> w_init;
	{
		istringstream www(buf);
		for(;;)
		{
			if(!(www>>temp))break;
			double value=atof(temp.c_str());
			w_init.push_back(value);
			totnum++;
		}
	}
	//check
	if(totnum!=total_param_count)
	{
		fprintf(stderr,"input_param %lld not equal to total_param_count %lld \n",
			totnum,total_param_count);
		exit(-1);
	}
	//assign
	for(U_INT i=0;i<totnum;i++)weights[i]=w_init[i];
	//OneDim_To_MultiDim
	OneDim_To_MultiDim(tot_layer,layer_count,layer_weight,weights);
}
//-> from other data-structure
void DeepCNF_Model::cnf_model_load(double *input,U_INT totnum)
{
	//check
	if(totnum!=total_param_count)
	{
		fprintf(stderr,"input_param %lld not equal to total_param_count %lld \n",
			totnum,total_param_count);
		exit(-1);
	}
	//assign
	for(U_INT i=0;i<totnum;i++)weights[i]=input[i];
	//OneDim_To_MultiDim
	OneDim_To_MultiDim(tot_layer,layer_count,layer_weight,weights);
}


//================ for prediction ==============//
void DeepCNF_Model::Pred_Single_Chain(
	vector <string> &input,vector <vector <double> > &output)
{
	int length_s=(int)input.size();
	DeepCNF_Seq *seq = new DeepCNF_Seq(length_s, this);
	seq->MAP_Probability(input,output);
	delete seq;
}


//================ for maximize AUC =============//__20150517__//
//-> init Chebyshev coefficients
int DeepCNF_Model::init_gamma_coeff(double *gamma_coeff,int maxdeg,double beta)
{
	//init judge beta
	double *cheby_coeff=new double[maxdeg+1];
	int idx;
	if(beta<=0) //-> use step function
	{
		Chebyshev chebyshev(beta);
		chebyshev.cheby_step_function(maxdeg,cheby_coeff);
		idx=chebyshev.cheby_convert(maxdeg,cheby_coeff,gamma_coeff);
	}
	else        //-> use sigmoid function
	{
		Chebyshev chebyshev(beta);
		chebyshev.cheby_expand(maxdeg,cheby_coeff);
		idx=chebyshev.cheby_convert(maxdeg,cheby_coeff,gamma_coeff);
	}
	//delete
	delete [] cheby_coeff;
	//return
	return idx;
}

//------ collect s_value out/sum -----//
//-> init s_value_out and v_value_out
void DeepCNF_Model::Initialize_SVvalue(void)
{
	int total_num=(max_degree + 2)*(max_degree + 1)/2;
	int totnum=total_num*dim_states;
	for(int i=0;i<totnum;i++)s_value_sum[i]=0;
	for(int i=0;i<totnum;i++)v_value_sum[i]=0;
	for(int i=0;i<totnum;i++)s_value_sum_n[i]=0;
	for(int i=0;i<totnum;i++)v_value_sum_n[i]=0;
}
//-> sum over s_value_out and v_value_out
void DeepCNF_Model::SumOver_SVvalue(void)
{
	int total_num=(max_degree + 2)*(max_degree + 1)/2;
	int totnum=total_num*dim_states;
	for(int i=0;i<totnum;i++)s_value_sum[i]+=s_value_out[i];
	for(int i=0;i<totnum;i++)v_value_sum[i]+=v_value_out[i];
	for(int i=0;i<totnum;i++)s_value_sum_n[i]+=s_value_out_n[i];
	for(int i=0;i<totnum;i++)v_value_sum_n[i]+=v_value_out_n[i];
}

//------ calculate AUC value ------//
//-> calculate AUC value based on s_value_sum and v_value_sum
double DeepCNF_Model::Calculate_AUC_Value(void)
{
	double ret_val=0;
	for(int k=0;k<dim_states;k++)
	{
		int idx=0;
		double total_value=0;
		for(int i = 0; i <= max_degree; i++)
		{
			for(int j = 0; j <= i; j++)
			{
				double gamma_value=gamma_coeff[idx];
				double upper_value=s_value_out[idx*dim_states+k]*v_value_out[idx*dim_states+k];
				double lower_value=s_value_out_n[idx*dim_states+k]*v_value_out_n[idx*dim_states+k];
				total_value+=gamma_value*upper_value/lower_value;
				idx++;
			}
		}
		auc_for_label[k]=total_value;
		ret_val+=total_value;
	}
	return ret_val;
}

