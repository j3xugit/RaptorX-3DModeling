import numpy as np
#import theano.tensor as T

## ResNet1DV21 is almost same as ResNet1D, but the former has removed unused BatchNormLayers and unused model parameters
## ResNet1DV22 has two BatchNormLayers in a ResBlock while ResNet1D and ResNet1DV21 has only one
## Experiments indicate that ResNet1D and ResNet1DV21 seem to work better, but it may depend on training algorithms and learning rate
allNetworks = ['ResNet1D', 'ResNet1DV21', 'ResNet1DV22']

allLabelTypes = ['vonMise2d', 'vonMise2d2', 'vonMise2d4', 'Gauss1d', 'Gauss2d', 'Gauss2d2', 'Gauss2d4', 'Discrete2C', 'Discrete3C', 'Discrete8C', 'Discrete18C' ]
allLabelNames = ['PhiPsi', 'SS3', 'SS8', 'ACC', 'ThetaTau', 'DISO', 'CLE']

#the true repsonse is the combination of one element in allLabelNames and one element in allLabelTypes
def Response2LabelType(response):
	return response.split('_')[1]

def Response2LabelName(response):
	return response.split('_')[0]

## the number of dimensions for a predicted value
responseValueDims = dict()
responseValueDims['Gauss1d'] = 1
responseValueDims['Gauss2d'] = 2
responseValueDims['Gauss2d2'] = 2
responseValueDims['Gauss2d4'] = 2

responseValueDims['vonMise2d'] = 2
responseValueDims['vonMise2d2'] = 2
responseValueDims['vonMise2d4'] = 2

responseValueDims['Discrete2C'] = 1
responseValueDims['Discrete3C'] = 1
responseValueDims['Discrete8C'] = 1
responseValueDims['Discrete18C'] = 1

## the number of paramters for the probability distribution function of a response
responseProbDims = dict()
responseProbDims['Gauss1d']=2
responseProbDims['Gauss2d']=5
responseProbDims['Gauss2d2']=2
responseProbDims['Gauss2d4']=4

responseProbDims['vonMise2d']=5
responseProbDims['vonMise2d2']=2
responseProbDims['vonMise2d4']=4

responseProbDims['Discrete2C']=2
responseProbDims['Discrete3C']=3
responseProbDims['Discrete8C']=8
responseProbDims['Discrete18C']=18

allAlgorithms = ['SGDM', 'SGDM2', 'SGNA', 'Adam', 'AdamW', 'AMSGrad', 'AdamWAMS']
allEmbeddingModes = ['SeqOnly', 'Seq+SS' ]
allKeywords = ['UseOneHotEncoding', 'UsePSFM', 'UseTemplate', 'UseSequenceEmbedding', 'activation', 'batchNorm', 'w4disorder' ]


## encoding 20 amino acids, only used to represent a primary sequence as a L*20 matrix
AAOrders = { 'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E' : 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19 }
AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
AAVectors = np.zeros((26,20)).astype(np.int32)

for aa in AAs:
        index = ord(aa) - ord('A')
        AAVectors[index][ AAOrders[aa] ] = 1

## conduct one-hot encoding of a protein sequence, which consists of a bunch of amino acids. Each amino acid is an element in AAs
def SeqOneHotEncoding(sequence):
        seq2int = (np.array(map(ord, sequence)) - ord('A') ).astype(np.int32)
        return AAVectors[seq2int]

"""
def InitializeModelSpecs():
	modelSpecs = dict()
        modelSpecs['trainFile'] = None
        modelSpecs['validFile'] = None
        modelSpecs['predFile'] = None
        modelSpecs['checkpointFile'] = None

        modelSpecs['network'] = 'ResNet1D'
        modelSpecs['responseStr'] = 'PhiPsi_vonMise'
        modelSpecs['responses'] = ['PhiPsi_vonMise']
	modelSpecs['w4responses'] = [1.]

        modelSpecs['algorithm'] = 'Adam'
        modelSpecs['numEpochs'] = [19, 2]
        modelSpecs['lrs'] = [np.float32(0.0002),  np.float32(0.0002)/10 ]

	modelSpecs['algorithm4var'] = 'Adam'
        modelSpecs['numEpochs4var'] = modelSpecs['numEpochs']
        modelSpecs['lrs4var'] = modelSpecs['lrs']

	modelSpecs['algStr'] = 'Adam:21:0.00022'

        modelSpecs['minibatchSize'] = 200
        modelSpecs['maxbatchSize'] = 1500

	modelSpecs['validation_frequency'] = 100
	modelSpecs['patience'] = 5

        ##default number of hidden units at 1d convolutional layer
        modelSpecs['conv1d_hiddens'] = [80, 100, 120, 140]
        modelSpecs['conv1d_repeats'] = [ 0,  0,  0,  0]

        ## for the logistic regression at the final stage
        modelSpecs['logreg_hiddens'] = [ ]

        modelSpecs['halfWinSize_seq'] = 5
        modelSpecs['activation'] = T.nnet.relu

        modelSpecs['L2reg'] = np.float32(0.0001)
	modelSpecs['w4diso'] = 3.

        modelSpecs['batchNorm'] = True
        modelSpecs['UseSampleWeight'] = True
	modelSpecs['UsePSFM'] = False
	modelSpecs['UseTemplate'] = False
	modelSpecs['UseSequenceEmbedding'] = False
	modelSpecs['UseOneHotEncoding'] = False

	return modelSpecs
"""
