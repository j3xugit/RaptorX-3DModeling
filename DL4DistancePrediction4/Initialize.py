import numpy as np

def InitializeModelSpecs():
	modelSpecs = dict()
	modelSpecs['Version'] = 3.0

	## if set to True, then train the model by minimizing the loss of real data minus the loss of reference data (which can be intepreted as energy function)
	## otherwise, train the model by minimizing the loss of real data
	modelSpecs['TrainByRefLoss'] = False

	## here 0.95 and 0.05 are the ratio of the overall data used for training and validation
        modelSpecs['trainFile'] = '0.95'
        modelSpecs['validFile'] = '0.05'
        modelSpecs['predFile'] = None
        modelSpecs['checkpointFile'] = None

	modelSpecs['ModelID'] = ''
        modelSpecs['network'] = 'DilatedResNet2D'

	## w4responses stores the weight factors for responses to be predicted
        modelSpecs['responseStr'] = 'CbCb:25C'
	modelSpecs['responses'] = ['CbCb_Discrete25C']
	modelSpecs['w4responses'] = [ 1. ]

	## the ratio of top predicted contacts to be evaluated
	modelSpecs['topRatios'] = [ 0.5 ]

	## algorithms, number of epochs and learning rates
	modelSpecs['algStr'] = 'AdamW:20+0.0002:1+0.00004'
        modelSpecs['algorithm'] = 'AdamW'
        modelSpecs['numEpochs'] = [ 20, 1 ]
        modelSpecs['lrs'] = [np.float32(0.0002), np.float32(0.0002)/5 ]

	## algorithms used to train a model for prediction of variance
	modelSpecs['algorithm4var'] = modelSpecs['algorithm']
        modelSpecs['numEpochs4var'] = modelSpecs['numEpochs']
	modelSpecs['lrs4var'] = modelSpecs['lrs']

	## how often to validate the trained model
        modelSpecs['validation_frequency'] = 100

	## when to stop training when validation loss does not have improved in the past patience epochs
        modelSpecs['patience'] = 5


        ##default number of hidden units at 1d convolutional layer
        modelSpecs['conv1d_hiddens'] = [30, 35, 40, 45]
        modelSpecs['conv1d_repeats'] = [ 0,  0,  0,  0]
	modelSpecs['conv1d_hwsz'] = 7

        ## the number of hidden units at 2d convolutional layer
        modelSpecs['conv2d_hiddens'] = [50, 55, 60, 65, 70, 75]
        modelSpecs['conv2d_repeats'] = [4,  4,  4,  4,  4,  4 ]
	modelSpecs['conv2d_hwszs'] =   [1,  1,  1,  1,  1,  1 ]
	modelSpecs['conv2d_dilations'] = [1, 1, 2,  4,  2,  1 ]

        ## for the logistic regression at the final stage
        modelSpecs['logreg_hiddens'] = [ 80 ]

	## half window size for 1d and 2d convolution kernel
        modelSpecs['halfWinSize_seq'] = 7
        modelSpecs['halfWinSize_matrix'] = 2
        modelSpecs['activation'] = 'RELU'

	## for conversion of sequential features to pairwise features
        modelSpecs['seq2matrixMode'] = {}
        modelSpecs['seq2matrixMode']['SeqOnly' ] = [ 4, 6, 12 ]
        modelSpecs['seq2matrixMode']['OuterCat' ] = [ 70, 35 ]

	## for the compression of the original matrix input, which has too many channels when full precision matrix or covariance matrix is used as input
	modelSpecs['hiddens4MatrixCompress'] = [100, 50]

	## L2 regularization factor
        modelSpecs['L2reg'] = 0.0001

	## the lower and upper bounds of a minibatch in terms of the number of residue pairs
        modelSpecs['minibatchSize'] = 60000
        modelSpecs['maxbatchSize'] = 350*350
	## the maximum offset from the diagonal line in sampling a submatrix
	modelSpecs['boundingBoxOffset'] = None

	## input features
	modelSpecs['UseSequentialFeatures'] = True
	modelSpecs['UseOneHotEncoding'] = True

	modelSpecs['UseSS'] = True
	modelSpecs['UseACC'] = True
	modelSpecs['UsePSSM'] = True
	modelSpecs['UseDisorder'] = False
	modelSpecs['UseCCMZ'] = False
	modelSpecs['UseRawCCM'] = True
	modelSpecs['ECInfo'] = 0

	##OtherPairs include mutual information and contact potential
	modelSpecs['UseOtherPairs'] = True
	modelSpecs['UseMI'] = True
	modelSpecs['UseContactPotential'] = True

        modelSpecs['UsePriorDistancePotential'] = False
        modelSpecs['UsePSICOV'] = False

	## bias added for smaller distance
        modelSpecs['LRbias'] = 'mid'

        ## by All, we consider all-range residue pairs including those pairs (i, j) with abs(i-j)<6
        modelSpecs['rangeMode'] = 'All'

        modelSpecs['batchNorm'] = True
        modelSpecs['SeparateTrainByRange'] = False

	## UseSampleWeight is usually set to True so that weight zero can be set for some special residue pair, e.g. a pair of two sequentially-adjacent residues or residues without valid 3D coordinates
        modelSpecs['UseSampleWeight'] = True
	## use the same weight for all labels
	modelSpecs['NoWeight4Label'] = True
	## use the same weight for all ranges
	modelSpecs['NoWeight4Range'] = True

	## do not use sequence separation of two residues/atoms as features
	modelSpecs['NoOldPosFeatures'] = True
	modelSpecs['UseNewPosFeatures'] = False
	modelSpecs['UseRgAsSeqFeature'] = False

	## distance threshold for orientation, i.e., when the distance is larger than this threshold, do not use the relative orientation of these two residues
	modelSpecs['distThreshold4Orientation'] = np.float32(20)

	modelSpecs['UseTemplate'] = False
	modelSpecs['Attention'] = None

	modelSpecs['NumDataLoaders'] = 2
	modelSpecs['QSize'] = 400

	return modelSpecs

