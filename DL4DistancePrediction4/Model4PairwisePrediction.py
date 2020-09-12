import os
import sys

import theano
import theano.tensor as T
from theano import printing
from theano.ifelse import ifelse

import numpy as np
from numpy import random as rng

import config
from Conv1d import Conv1DLayer
from NN4LogReg import NN4LogReg
from NN4Normal import NN4Normal
from ResNet4Distance import ResNet
from DilatedResNet4Distance import DilatedResNet
from utils import OuterConcatenate, MidpointFeature
import DistanceUtils

from config import Response2LabelName, Response2LabelType, ParseResponse, GetResponseProbDims, GetResponseValueDims

##this class calculates the convolution of 1D sequence and then convert it to a 2D matrix using OuterConcatenate
class Conv1D2Matrix:

    def __init__(self, rng, input, n_in, n_hiddens=[], halfWinSize=0, mask=None):
        ## input has shape (batchSize, seqLen, n_in) where n_in is the number of features at each position
        ## all the sequences aligned at the right hand side, mask is used for the left hand side positions
        ## mask has shape (batchSize, #positions_to_be_masked ) 

        ## n_hiddens specify the number of hidden neurons at all the hidden layers
	## halfWinSize specify the half window size of the convolution operation

        self.input = input
        self.mask = mask
	self.n_in = n_in
	self.halfWinSize = halfWinSize

	self.layers = []
        self.params = []
        self.paramL1 = 0
        self.paramL2 = 0

        n_out_at_last_layer = n_in
        output_at_last_layer = input

        for i in xrange(len(n_hiddens)):

            layer = Conv1DLayer(rng, input=output_at_last_layer, numOfInFeatures=n_out_at_last_layer, numOfOutFeatures=n_hiddens[i], halfWinSize=halfWinSize, mask=mask)
            self.layers.append(layer)
            n_out_at_last_layer = n_hiddens[i]
            output_at_last_layer = layer.output

            self.params += layer.params
            self.paramL1 += layer.paramL1
            self.paramL2 += layer.paramL2


	## conv_out has shape (batchSize, seqLen, n_out_at_last_layer ), it is the final convolution result
        conv_out = output_at_last_layer

        ## now we take the OuterConcatenate operation, the result has shape (batchSize, seqLen, seqLen, 2 * n_out_at_last_layer)
        self.output = OuterConcatenate (input = conv_out)
        self.n_out = 2 * n_out_at_last_layer

	##add the below two lines for convenience
	self.conv_out = conv_out
	self.n_conv_out = n_out_at_last_layer


def TestConv1D2Matrix():
    x = T.tensor3('x')
    conv1dlayer = Conv1D2Matrix( rng=np.random.RandomState(), input=x, n_in=2, n_hiddens=[3], halfWinSize=1)
    f = theano.function([x], conv1dlayer.output)

    a = np.random.uniform(0, 1, (2, 7, 2)).astype(np.float32)
    print a
    print f(a)

## 2D convolution on a matrix for distance matrix prediction
## no pooling is allowed here
class Conv2D4DistMatrix:
    def __init__(self, rng, input, n_in=1, n_hiddens=None, halfWinSize=0, mask=None, activation=T.nnet.relu):
	## input has shape (batchSize, nRows, nCols, n_in) where n_in is the number of input features at each element
	## here we assume that input[0], input[1], ...., input[batchSize -1 ] are aligned at the right bottom
        ## nrows and ncols are unknown
	## mask has shape (batchSize, #rows_to_be_masked, nCols), mask is a binary matrix with 0 indicates padding positions and 1 real data position

        ##hwsz is the half window size of the filter shape
	##mask is used to reduce noise introduced by padding. all the matrices in a batch shall be aligned at the right bottom corner 

	self.input = input
 	self.n_in = n_in
	self.n_hiddens = n_hiddens
	self.halfWinSize = halfWinSize
	self.mask = mask

        ##the window size of our filters is always odd
        wSize = 2*halfWinSize + 1
        
	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0


        ## reshape input so that it has shape (batchSize, n_in, nRows, nCols)
        output_in_last_layer = input.dimshuffle(0, 3, 1, 2)

        n_out_in_last_layer = n_in

        for i in xrange(len(n_hiddens)):
            W_shape = (n_hiddens[i], n_out_in_last_layer, wSize, wSize)
	
	    W_values = np.asarray(
                    rng.uniform( low = - np.sqrt(6. / (n_out_in_last_layer * wSize * wSize + n_hiddens[i])),
                                 high =  np.sqrt(6. / (n_out_in_last_layer * wSize * wSize + n_hiddens[i])),
		                 size = W_shape
                    ),
                    dtype = theano.config.floatX
                )
            if activation == T.nnet.relu:
	        W_values = np.asarray(
                    rng.normal( scale = np.sqrt(2. / (n_out_in_last_layer * wSize * wSize + n_hiddens[i])), size = W_shape),
                    dtype = theano.config.floatX
                )

	    b_shape = (n_hiddens[i], )
            b_values = np.asarray (rng.uniform(low = -.0, high =.0, size = b_shape), dtype=theano.config.floatX )

	    W = theano.shared (value = W_values, name = 'conv2d_W_' + str(i), borrow = True)
	    b = theano.shared (value = b_values, name = 'conv2d_b_' + str(i), borrow = True)

            shape_in_last_layer = (None, n_out_in_last_layer, None, None)

	    ## Please make sure that conv2d_out has shape (batchSize, n_hiddens[i], nRows, nCols)
            conv2d_out = T.nnet.conv2d(output_in_last_layer, W, input_shape = shape_in_last_layer, filter_shape=W_shape, border_mode='half')
	    if activation is not None:
	        out2 = activation(conv2d_out + b.dimshuffle('x', 0, 'x', 'x'))
	    else:
	        out2 = (conv2d_out + b.dimshuffle('x', 0, 'x', 'x'))


	    if mask is not None:
	        ## mask has shape (batchSize, #rows_to_be_masked, nCols)

                ## a subtensor of out2 along the horiz direction
                out2_sub_horiz = out2[:, :, :mask.shape[1], :]
                mask_horiz = mask.dimshuffle(0, 'x', 1, 2)
                out3 = T.set_subtensor(out2[:, :, :mask.shape[1], :], T.mul(out2_sub_horiz, mask_horiz) )

                ## a subtensor of out3 along the vertical direction
                out3_sub_vertical = out3[:, :, :, :mask.shape[1] ]
                mask_vertical = mask.dimshuffle(0, 'x', 2, 1)
		out4 = T.set_subtensor(out3[:, :, :, :mask.shape[1] ], T.mul(out3_sub_vertical, mask_vertical) )

                output_in_last_layer = out4

            else:
	        output_in_last_layer = out2

	    n_out_in_last_layer = n_hiddens[i]

	    self.params += [W, b]
            self.paramL1 += (abs(W).sum() +  abs(b).sum())
	    self.paramL2 += ((W**2).sum() + (b**2).sum())
	    

        ##change the shape of the final output to (batchSize, nRows, nCols, n_out)
        self.output = output_in_last_layer.dimshuffle(0, 2, 3, 1)
        self.n_out = n_out_in_last_layer

def TestConv2D4DistMatrix():
    x = T.tensor4('x')
    conv2dlayer = Conv2D4DistMatrix( rng=np.random.RandomState(), input=x, n_in=2, n_hiddens=[3], halfWinSize=1)
    f = theano.function([x], conv2dlayer.output)

    a = np.random.uniform(0, 1, (2, 7, 7, 2)).astype(np.float32)
    print a
    print f(a)

"""
def SeparateTrainByRange(rng, conv_out, n_in, numLabels, n_hiddens_logreg=None):
    
    batchSize, dataLen, dataLen, tmp_n = conv_out.shape

    n_in4logreg = n_in

    ##theano currently does not have triu for tensor3, so we have to do some work here
    M1s = T.ones((dataLen, dataLen), dtype=np.int8)
    M1s_3d = T.ones((batchSize, dataLen, dataLen), dtype=np.int8)
    Sep24Mat = T.triu(M1s, 24) + T.tril(M1s, -24)
    Sep12Mat = T.triu(M1s, 12) + T.tril(M1s, -12)
    Sep6Mat = T.triu(M1s, 6) + T.tril(M1s, -6)

    LRselection = Sep24Mat.dimshuffle('x',0,1)
    MRselection = (Sep12Mat - Sep24Mat).dimshuffle('x', 0, 1)
    SRselection = (Sep6Mat - Sep12Mat ).dimshuffle('x', 0, 1)
    LRselection = T.mul( LRselection, M1s_3d )
    MRselection = T.mul( MRselection, M1s_3d )
    SRselection = T.mul( SRselection, M1s_3d )

    selections = [ LRselection, MRselection, SRselection ]

    regressors = []
    output_2d = T.ones((batchSize, dataLen, dataLen), dtype=np.int8)*(-1)
    output_2d_prob = T.zeros((batchSize, dataLen, dataLen, numLabels), dtype=theano.config.floatX)

    for sel in selections:
        selected = conv_out[sel.nonzero()]
        regressor = MLLogReg(rng=rng, input=selected, n_in=n_in4logreg, n_out=numLabels, n_hiddens=n_hiddens_logreg)
        regressors.append(regressor)
        output_2d = T.set_subtensor(self.output_2d[ sel.nonzero() ], T.cast(regressor.y_pred, 'int8') )
        output_2d_prob = T.set_subtensor(self.output_2d_prob[ sel.nonzero() ], regressor.p_y_given_x )

    return selections, regressors, output_2d, output_2d_prob
"""

## this is the main class for the prediction of inter-residue distance, hydrogen bonding and beta-pairing
class ResNet4DistMatrix:

    def __init__(self, rng, seqInput, matrixInput, mask_seq=None, mask_matrix=None, embedInput=None, boundingbox=None, modelSpecs=None):
	"""
	seqInput has shape (batchSize, seqLen, n_in_seq)
	matrixInput has shape (batchSize, seqLen, seqLen, n_in_matrix)
	mask_seq has shape (batchSize, #cols_to_be_masked)
        mask_matrix has shape (batchSize, #rows_to_be_masked, seqLen)
	embedInput has shape (batchSize, seqLen, n_in2)
	boundingbox is a vector of 4 integer elements: top, left, bottom and right. boundingbox shall only be applied to the matrix converted from sequential features.
        """

	assert (modelSpecs is not None)

	self.modelSpecs = modelSpecs	
	self.responses = modelSpecs['responses']

	## set the number of hidden neurons and number of layers
	n_in_seq = modelSpecs['n_in_seq']
	n_in_matrix = modelSpecs['n_in_matrix']
	n_hiddens_seq = modelSpecs['conv1d_hiddens']
	n_hiddens_matrix = modelSpecs['conv2d_hiddens']
	n_hiddens_logreg = modelSpecs['logreg_hiddens']
	seq_repeats = modelSpecs['conv1d_repeats']
	matrix_repeats = modelSpecs['conv2d_repeats']

	## half win size for convolutional operation
	if modelSpecs['network'].startswith('DilatedResNet'):
		hwsz_matrix = modelSpecs['conv2d_hwszs']
		hwsz_seq = [modelSpecs['conv1d_hwsz'] ]* len(n_hiddens_seq)
		dilation_seq = [1] * len(n_hiddens_seq)
		dilation_matrix = modelSpecs['conv2d_dilations']
	else:
		hwsz_matrix=modelSpecs['halfWinSize_matrix']
		hwsz_seq=modelSpecs['halfWinSize_seq']


	## masks to reduce impact of padding zeros	
	self.mask_1d = mask_seq
	self.mask_2d = mask_matrix

	self.layers = []

	act = T.nnet.relu
	if modelSpecs['activation'] == 'TANH':
		act = T.tanh
        
	# sequence convolution 
	if modelSpecs['network'].startswith('DilatedResNet'):
        	#seqConv = DilatedResNet(rng, input=seqInput, n_in=n_in_seq, n_hiddens=n_hiddens_seq, n_repeats=seq_repeats, halfWinSize=hwsz_seq, dilation=dilation_seq, mask=mask_seq, activation=act, batchNorm=modelSpecs['batchNorm'], version=modelSpecs['network'])
        	seqConv = DilatedResNet(rng, input=seqInput, n_in=n_in_seq, n_hiddens=n_hiddens_seq, n_repeats=seq_repeats, halfWinSize=hwsz_seq, dilation=dilation_seq, mask=mask_seq, activation=act, modelSpecs=modelSpecs)
	else:
        	seqConv = ResNet(rng, input=seqInput, n_in=n_in_seq, n_hiddens=n_hiddens_seq, n_repeats=seq_repeats, halfWinSize=hwsz_seq, mask=mask_seq, activation=act, batchNorm=modelSpecs['batchNorm'], version=modelSpecs['network'])
	self.layers.append(seqConv)


	## transform 1d sequence to 2d matrix
        seq2matrixMode = modelSpecs['seq2matrixMode']
	seq2matrixLayers = []
	embedLayers = []

        ## determine if we shall use the sequential features or not. The sequential features include sequence profile (PSSM), predicted secondary structure and predicted solvent accessibility
	## useSequentialFeatures is True by default
	##useSequentialFeatures = ( modelSpecs.has_key('UseSequentialFeatures') and (modelSpecs['UseSequentialFeatures'] is True) )

	## use OuterConcatenation operation to convert sequence features into pairwise features
	if seq2matrixMode.has_key('OuterCat') and config.UseSequentialFeatures:

	    ##midpointfeature has shape (batchSize, seqLen, seqLen, n_midpoint_out)
	    midpointfeature, n_midpoint_out = MidpointFeature(seqConv.output, seqConv.n_out, box=boundingbox)

	    ##remove noise in midpointfeature
            ## mask_matrix is used to reduce noise introduced by padding positions
	    mid_subtensor = midpointfeature[:, :mask_matrix.shape[1], :, :]
	    midpointfeature = T.set_subtensor(mid_subtensor, T.mul(mask_matrix.dimshuffle(0,1,2,'x'), mid_subtensor) )
	    mid_subtensor2 = midpointfeature[:, :, :mask_matrix.shape[1], : ]
	    midpointfeature = T.set_subtensor(mid_subtensor2, T.mul(mask_matrix.dimshuffle(0,2,1,'x'), mid_subtensor2) )

	    ## here we use convolution with halfWinSize=0 to reduce model complexity
	    compressLayer = Conv2D4DistMatrix(rng, input=midpointfeature, n_in=n_midpoint_out, n_hiddens=seq2matrixMode['OuterCat'], halfWinSize=0, mask=mask_matrix )
	    #compressLayer = Conv2D4DistMatrix(rng, input=midpointfeature, n_in=n_midpoint_out, n_hiddens=seq2matrixMode['OuterCat'], halfWinSize=0, mask=None )
	    seq2matrixLayers.append(compressLayer)

	## embedding primary sequence and/or predicted secondary structure
	if embedInput is not None:
	    from EmbeddingLayer import EmbeddingLayer4AllRange

	    if seq2matrixMode.has_key('Seq+SS'):
	        n_out_embed = seq2matrixMode['Seq+SS']
	    elif seq2matrixMode.has_key('SeqOnly'):
	        n_out_embed = seq2matrixMode['SeqOnly']
	    else:
		print 'At least one of two embedding modes Seq+SS or SeqOnly shall be specified.'
		exit(1)

	    embedLayer = EmbeddingLayer4AllRange(embedInput, modelSpecs['n_in_embed'], n_out_embed, box=boundingbox)
	    seq2matrixLayers.append(embedLayer)
	    embedLayers.append(embedLayer)

	"""
	we do not use this profile embedding any more
	## embedding the sequence profile
	if seq2matrixMode.has_key('Profile') and useSequentialFeatures:
	    from EmbeddingLayer import ProfileEmbeddingLayer
	    pEmbedLayer = ProfileEmbeddingLayer(seqConv.output, seqConv.n_out, seq2matrixMode['Profile'])
	    seq2matrixLayers.append(pEmbedLayer)
	    embedLayers.append(pEmbedLayer)
	"""

        self.layers += seq2matrixLayers

	bUseCCMFnorm, bUseCCMsum, bUseCCMraw, bUseFullMI, bUseFullCov = config.ParseExtraCCMmode(modelSpecs)
	if (bUseCCMraw or bUseFullMI or bUseFullCov) and config.CompressMatrixInput(modelSpecs):
		## here we add a compress layer to reduce the #channels of the original matrix input. 
		n_hiddens4MatrixCompress = modelSpecs['hiddens4MatrixCompress']
		compressLayer4MatrixInput = Conv2D4DistMatrix(rng, input=matrixInput, n_in=n_in_matrix, n_hiddens=n_hiddens4MatrixCompress, halfWinSize=0, mask=mask_matrix )
		compressedMatrixInput = compressLayer4MatrixInput.output
		n_compressedMatrix = compressLayer4MatrixInput.n_out
		input_2d = T.concatenate( [ compressedMatrixInput ] + [ layer.output for layer in seq2matrixLayers], axis=3)
		n_input2d = n_compressedMatrix + sum([ layer.n_out for layer in seq2matrixLayers ])
	else:

		##old code for merging original matrix input and sequential input
		input_2d = T.concatenate( [ matrixInput ] + [ layer.output for layer in seq2matrixLayers], axis=3)
		n_input2d = n_in_matrix + sum([ layer.n_out for layer in seq2matrixLayers ])

	#print 'n_input2d=', n_input2d

	if modelSpecs['network'].startswith('ResNet'):
		matrixConv=ResNet(rng, input=input_2d, n_in=n_input2d, n_hiddens=n_hiddens_matrix, n_repeats=matrix_repeats, halfWinSize=hwsz_matrix, mask=mask_matrix, activation=act, batchNorm=modelSpecs['batchNorm'], version=modelSpecs['network'])

	elif modelSpecs['network'].startswith('DilatedResNet'):
		#matrixConv=DilatedResNet(rng, input=input_2d, n_in=n_input2d, n_hiddens=n_hiddens_matrix, n_repeats=matrix_repeats, halfWinSize=hwsz_matrix, dilation=dilation_matrix, mask=mask_matrix, activation=act, batchNorm=modelSpecs['batchNorm'], version=modelSpecs['network'])
		matrixConv=DilatedResNet(rng, input=input_2d, n_in=n_input2d, n_hiddens=n_hiddens_matrix, n_repeats=matrix_repeats, halfWinSize=hwsz_matrix, dilation=dilation_matrix, mask=mask_matrix, activation=act, modelSpecs=modelSpecs)
	else:
		print 'ERROR: Unimplemented deep network type: ', modelSpecs['network']
		exit(1)

	self.layers.append(matrixConv)

        conv_out = matrixConv.output

        selected = conv_out.dimshuffle(3, 0, 1, 2).flatten(2).dimshuffle(1, 0) 
	n_in4logreg = matrixConv.n_out 

	self.outputList = []
	self.output_probList = []
	self.predictors = []

	self.params4var = []
        self.paramL14var = 0
        self.paramL24var = 0

	for res in modelSpecs['responses']:

		labelType = Response2LabelType(res)
		predictor = None

		if labelType.startswith('Discrete'):
			assert GetResponseValueDims(res) == 1
                        predictor = NN4LogReg(rng=rng, input=selected, n_in=n_in4logreg, n_out=GetResponseProbDims(res), n_hiddens=n_hiddens_logreg)

		elif labelType.startswith('LogNormal') or labelType.startswith('Normal'):
			predictor = NN4Normal(rng=rng, input=selected, n_in=n_in4logreg, n_variables=GetResponseValueDims(res), n_out=GetResponseProbDims(res), n_hiddens=n_hiddens_logreg)

			## recording parameters specific for variance prediction
                        self.params4var += predictor.params4var
                        self.paramL14var += predictor.paramL14var
                        self.paramL24var += predictor.paramL24var

		else:
			print 'incorrect response name or label type: ', res
                        exit(1)

                self.layers.append(predictor)
                self.predictors.append(predictor)

            	## output in 2d matrix
	    	output_2d = predictor.y_pred.reshape( (conv_out.shape[0], conv_out.shape[1], conv_out.shape[2], GetResponseValueDims(res) ) )
	    	output_2d_prob = predictor.output.reshape( (conv_out.shape[0], conv_out.shape[1], conv_out.shape[2], GetResponseProbDims(res) ) )

	    	self.outputList.append( output_2d )
	    	self.output_probList.append( output_2d_prob )

	self.output = T.concatenate( self.outputList, axis=3 )
	self.output_prob = T.concatenate( self.output_probList, axis=3 ) 

	## collect all the model parameters and their norms
        self.params = []
	self.paramL2 = 0
	self.paramL1 = 0

	for layer in self.layers:
	    self.params  += layer.params
	    self.paramL2 += layer.paramL2
	    self.paramL1 += layer.paramL1

	"""
	## we do not use it temporarily
	## regularize parameters used in embedding to force the center of all the embedding vectors to be 0
	for embedLayer in embedLayers:
	    self.paramL2 += embedLayer.pcenters
	"""

    ## zList is a list of z and weightList is a list of w, each (z, w) pair corresponds to one response to be predicted
    ## w may not have the same shape as z, it contains weight for each element in z
    ## w has the same shape, i.e., (batchSize, seqLen, seqLen) , while z may have shape (batchSize, seqLen, seqLen, responseValueDims[response]) or (batchSize, seqLen, seqLen) if responseValueDims[response]==1
    ## this functions returns a vector of losses, in which each element is the loss of one response
    def loss(self, zList, weightList=None, useMeanOnly=False, trainByRefLoss=None):
	losses = []
        if weightList is not None and len(weightList)>0:
        	for predictor, z, w in zip(self.predictors, zList, weightList):
			if (z.ndim == 3 ):
                		zflat = z.flatten().dimshuffle(0, 'x')
			elif (z.ndim == 4):
                		zflat = z.dimshuffle(3, 0, 1, 2).flatten(2).dimshuffle(1, 0)
			else:
				print 'unsupported ndim for z in errors():', z.ndim
				exit(1)
			assert (w.ndim == 3)
                        wflat = w.flatten().dimshuffle(0, 'x')

			## zflat shall have shape (batchSize, valueDims) and wflat shall have shape (batchSize, 1)
                        losses.append(  predictor.loss(zflat, useMeanOnly=useMeanOnly, sampleWeight=wflat)   )
        else:
                for predictor, z in zip(self.predictors, zList):
			if (z.ndim == 3):
                		zflat = z.flatten().dimshuffle(0, 'x')
			elif (z.ndim == 4):
                		zflat = z.dimshuffle(3, 0, 1, 2).flatten(2).dimshuffle(1, 0)
			else:
				print 'unsupported ndim for z in errors():', z.ndim
				exit(1)

			## zflag shall have shape (batchSize, valueDims)
                        losses.append(  predictor.loss(zflat, useMeanOnly=useMeanOnly )  )

	if trainByRefLoss is not None:
		losses = [ trainByRefLoss * los for los in losses ]
	
	## losses is a list of scalar
        return T.stack(losses)

    ## calculate errors for one response with label Type Discrete12C, Discrete25C, Discrete52C, return a tensor with ndim=1
    def errors4one(self, z, out, weight=None, distLabelType='12C'):
	distBins = config.distCutoffs[distLabelType]
	label8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, distBins)
	label15 = DistanceUtils.LabelsOfOneDistance(config.InteractionLimit, distBins)

	z3C = T.cast( T.ge(z, label8), 'int32') + T.cast( T.ge(z, label15), 'int32')
	o3C = T.cast( T.ge(out, label8), 'int32') + T.cast( T.ge(out, label15), 'int32')

	if weight is not None:
            err = T.sum( T.mul(weight, T.neq(o3C, z3C) ) )*1./T.sum(weight)
	else:
            err = T.mean( T.neq(o3C , z3C) ) 

	## err is s scalar, convert it to a tensor with ndim=1
	return T.stack([err] )

    ## this function returns a vector of errors, the size of this vector is equal to the sum of ValueDims for all the responses
    def errors(self, zList, weightList=None):
	errs = []
	if weightList is not None and len(weightList)>0:
		for res, predictor, z, w, o in zip(self.responses, self.predictors, zList, weightList, self.outputList):
			labelType = Response2LabelType(res)
			numLabels = GetResponseProbDims(res)

		        ## if the label type is Discrete25C, Discrete52C, Discrete12C
			if res in config.allAtomPairNames and labelType.startswith('Discrete') and numLabels > 3:
				assert (z.ndim == 3 and GetResponseValueDims(res) == 1 )
				o2 = o.flatten(3)
				## here we convert 12C, 25C, and 52C to 3C for error calculation, which makes the result easier to interpret
				errs.append( self.errors4one(z, o2, weight=w, distLabelType=labelType[len('Discrete'): ] ) )
			else:
				## call the error function of each predictor
				if (z.ndim == 3 ):
                			zflat = z.flatten().dimshuffle(0, 'x')
				elif (z.ndim == 4 ):
                			zflat = z.dimshuffle(3, 0, 1, 2).flatten(2).dimshuffle(1, 0)
				else:
					print 'unsupported ndim for z in errors():', z.ndim
					exit(1)

				assert (w.ndim == 3)
                        	wflat = w.flatten().dimshuffle(0, 'x')
                                e = predictor.errors(zflat, sampleWeight=wflat)
				## e is a tensor with ndim=1
                                errs.append(e)

	else:
		for res, predictor, z, o in zip(self.responses, self.predictors, zList, self.outputList):
			labelType = Response2LabelType(res)
			numLabels = GetResponseProbDims(res)
			if res in config.allAtomPairNames and labelType.startswith('Discrete') and numLabels > 3 :
				assert (z.ndim == 3 and GetResponseValueDims(res) == 1 )
				o2 = o.flatten(3)
				errs.append( self.errors4one(z, o, distLabelType=labelType[len('Discrete'): ] ) )
			else:
				## call the error function of each predictor
				if (z.ndim == 3):
                			zflat = z.flatten().dimshuffle(0, 'x')
				elif (z.ndim == 4):
                			zflat = z.dimshuffle(3, 0, 1, 2).flatten(2).dimshuffle(1, 0)
				else:
					print 'unsupported ndim for z in errors():', z.ndim
					exit(1)
                                e = predictor.errors(zflat)
				## e is a tensor with ndim=1
                                errs.append(e)

	return T.concatenate(errs)

    """
    ## calculate the confusion matrix, not carefully checked, so please do not use it
    def confusionMatrix(self, zList):

	## here we consider only the first element of Zs to simplify the coding
	#z = Zs[0]

	##this function shall not be used any more
	def confusionMatrix6CSum(pred, truth):
	    truth_new = T.cast( T.gt(truth, 9) + T.gt(truth, 24), 'int32')
	    pred_new = T.cast( T.gt(pred, 9) + T.gt(pred, 24), 'int32')
	    return confusionMatrix3C(pred_new, truth_new)
	
	## both pred and truth have shape (batchSize, seqLen, seqLen), so is z
	def confusionMatrix12C(pred, truth):
	    ##convert pred, truth to 3C

	    distBins = config.distCutoffs[self.modelSpecs['distLabelType']]
	    label8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, distBins)
	    label15 = DistanceUtils.LabelsOfOneDistance(config.InteractionLimit, distBins)

            truth1 = T.cast( T.ge(truth, label8), 'int32')
            truth2 = T.cast( T.ge(truth, label15), 'int32')
            truth_new = truth1 + truth2
           
	    pred1 = T.cast( T.ge(pred, label8), 'int32')
	    pred2 = T.cast( T.ge(pred, label15), 'int32')
	    pred_new = pred1 + pred2

	    return confusionMatrix3C(pred_new, truth_new)

	## both pred and truth has shape (batchSize, seqLen, seqLen)
	def confusionMatrix3C(pred, truth):

	    selMatrix = T.ones_like(zList[0])
	    if self.mask_2d is not None:
	        mshape = self.mask_2d.shape
	        selMatrix = T.set_subtensor(selMatrix[:, :mshape[1], :], self.mask_2d)
	        selMatrix = T.set_subtensor(selMatrix[:, :, :mshape[1]], self.mask_2d.dimshuffle(0,2,1) )

            dataLen = truth.shape[1]
            triu_3d_matrices = []	
	    for sep in [24, 12, 6]:
	        triu_3d_matrices.append( T.triu(T.ones((dataLen, dataLen), dtype=np.int32), sep).dimshuffle('x', 0, 1) )

	    ##mask for long-range, medium-range and short-range
	    triu_LR = triu_3d_matrices[0]
	    triu_MR = triu_3d_matrices[1] - triu_3d_matrices[0]
	    triu_SR = triu_3d_matrices[2] - triu_3d_matrices[1]

	    ##theano.tensor does not have the triu operation for a 3D tensor, so we have to deal with this by first applying triu to a 2D tensor,
	    ##and then multiply the result with a 3D tensor (selMatrix) here to obtain the effect of 3D triu

	    pred_truth = truth * 3 + pred
	
	    confMs = []
	    for t in [ triu_LR, triu_MR, triu_SR ]:
	        sel = T.mul( t, selMatrix)
	        pred_truth_selected = pred_truth[ sel.nonzero() ]
	        predcount = T.bincount( pred_truth_selected, minlength=9 ).reshape((3,3))	
	        confMs.append(predcount)

	    return T.stacklists( confMs )

	distLabelType = self.modelSpecs['distLabelType']
	confusionMatrices = []

	if distLabelType == '3C':
		for out, z in zip(self.outputList, zList):
			confusionMatrices.append( confusionMatrix3C(out, z) )
	else:
		for out, z in zip(self.outputList, zList):
	    		confusionMatrices.append( confusionMatrix12C(out, z) )

	return T.stacklists(confusionMatrices)
    """
	    
    ## TopAccuracyByRange is only used for validation, but not for training. It is correct only when the input matrix is cut off along the diagonal line.
    ## this function returns the long-, medium-, long+medium- and short-range prediction accuracy when the top proteinLen*topRatio predicted contacts are evaluated
    ## Note that the average accuracy is only for one batch. To calculate the avg accuracy of multiple batches, u may have to take into consideration
    ## the fact that two batches may have different numbers of proteins
    ## this function returns a 4*2 matrix where rows 0, 1, 2 and 3 correspond to long-, medium-, long+medium- and short-range accuracy.
    ## column 0 is the accuracy and column 1 is the average number of predicted contacts among all the top proteinLen*topRatio predicted contacts
    ## the average number of predicted contacts may not be equal to avg(proteinLen)* topRatio when topRatio is big and proteinLen is small

    def TopAccuracyByRange(self, zList):

	currentResponse = None
	topRatio = 0.5

	## in this function, we assume that pred is a tensor3 of floatX and truth is a matrix
	## pred has shape (dataLen, dataLen, 2) and truth has shape (dataLen, dataLen)
	## we also assume that label 0 is positive and label 1 is negative
	## the result is not 100% accurate for non-symmetric response, e.g., hydrogen-bonding matrix
	def TopAccuracy2C(pred=None, truth=None, symmetric=False):

            	M1s = T.ones_like(truth, dtype=np.int8)
            	LRsel = T.triu(M1s, 24)
            	MLRsel = T.triu(M1s, 12)
	    	SMLRsel = T.triu(M1s, 6)
	    	MRsel = MLRsel - LRsel
	    	SRsel = SMLRsel - MLRsel
            
	    	dataLen = truth.shape[0]

		pred0 = pred[:,:,0]

		if symmetric:
	    		avg_pred = (pred0 + pred0.dimshuffle(1, 0) )/2.0
		else:
			avg_pred = pred0
		
            	#pred_truth = T.concatenate( (avg_pred, truth.dimshuffle(0, 1, 'x') ), axis=2)
            	pred_truth = T.stack( [avg_pred, T.cast(truth, 'int32')], axis=2 )

	    	accuracyList = []
	    	for Rsel in [ LRsel, MRsel, MLRsel, SRsel]:
                	selected_pred_truth = pred_truth[Rsel.nonzero()]

			## sort by the predicted value for label 0 from the largest to the smallest
	        	selected_pred_truth_sorted = selected_pred_truth[ (selected_pred_truth[:,0]).argsort()[::-1] ]

	        	#print 'topRatio =', topRatio
	        	numTops = T.minimum( T.iround(dataLen * topRatio), selected_pred_truth_sorted.shape[0])

			selected_sorted_truth = T.cast(selected_pred_truth_sorted[:, -1], 'int32')
                	numTruths = T.bincount(selected_sorted_truth, minlength=2 )
	        	numCorrects = T.bincount( selected_sorted_truth[0:numTops], minlength=2 )
			#numTops = T.minimum(numTops, numTruths[0])
	        	accuracyList.append( T.stack( [numCorrects[0] *1./(numTops + 0.001), numTops, numTruths[0] ], axis=0)  )

            	return T.stacklists( accuracyList )

	def TopAccuracyNormal(pred=None, truth=None, symmetric=True):
		truth_new = T.ge(truth, config.ContactDefinition )
		
		if pred.ndim == 2:
			pred_new = -pred.dimshuffle(0, 1, 'x')
		else:
			pred_new = -pred

		return TopAccuracy2C(pred = pred_new, truth = truth_new, symmetric=symmetric)

	def TopAccuracyLogNormal(pred=None, truth=None, symmetric=True):
		truth_new = T.ge(truth, T.log(config.ContactDefinition) )
		
		if pred.ndim == 2:
			pred_new = -pred.dimshuffle(0, 1, 'x')
		else:
			pred_new = -pred

		return TopAccuracy2C(pred = pred_new, truth = truth_new, symmetric=symmetric)

	## in this function, we assume that pred is tensor3 of float and truth is a matrix of int8 or int
	## pred has shape (dataLen, dataLen, numLabels), having the predicted probability of each label
	## truth has shape (dataLen, dataLen)
        def TopAccuracyMultiC(pred=None, truth=None, subType=None, symmetric=True):
            	## convert pred and truth to 2C

	    	distBins = config.distCutoffs[subType]
            	label8 = DistanceUtils.LabelsOfOneDistance(config.ContactDefinition, distBins)
            	label15 = DistanceUtils.LabelsOfOneDistance(config.InteractionLimit, distBins)

            	truth1 = T.cast( T.ge(truth, label8), 'int32')
            	truth_new = truth1
           
            	pred1 = T.sum(pred[:,:,:label8], axis=2, keepdims=True)
            	pred2 = T.sum(pred[:,:,label8:], axis=2, keepdims=True)
	    	pred_new = T.concatenate( (pred1, pred2), axis=2)

	    	return TopAccuracy2C(pred=pred_new, truth=truth_new, symmetric=symmetric)

	## this function calculates the top accuracy of predicted orientation angles. Here the accuracy is defined as the percentage of residue pairs that is correctly to be predicted to have a valid orientation angle
	## this is not a good definition, but this is what we can do here
	## pred has shape (dataLen, dataLen, numLabels), having the predicted probability of each label
	## truth has shape (dataLen, dataLen)
	def TopAccuracyOrientation(pred=None, truth=None, largestValidLabel=36, symmetric=True):

                pred1 = T.sum(pred[:,:,:largestValidLabel], axis=2, keepdims=True)
                pred2 = T.sum(pred[:,:,largestValidLabel:], axis=2, keepdims=True)
                pred_new = T.concatenate( (pred1, pred2), axis=2)

		truth_new = T.cast(T.ge(truth, largestValidLabel), 'int32')

                return TopAccuracy2C(pred=pred_new, truth=truth_new, symmetric=symmetric)


	##def EvaluateAccuracy(pred_prob, truth, pad_len, response):
	def EvaluateAccuracy(pred_prob, truth, pad_len):
	    	pred_in_correct_shape = T.cast( pred_prob[pad_len:, pad_len: ], dtype=theano.config.floatX)
	    	truth_in_correct_shape = truth[pad_len:, pad_len: ]
	
		labelName, labelType, subType = ParseResponse(currentResponse)
		symmetric = config.IsSymmetricLabel(labelName)

		if labelName in config.allOrientationNames:
			if not config.IsDiscreteLabel(labelType):
				print 'ERROR: unsupported label type for orientation matrix prediction: ', currentResponse
				exit(1)

			numLabels = GetResponseProbDims(currentResponse)
			if subType.endswith('Plus') or subType.endswith('Minus'):
				largestValidLabel = numLabels-2
			else:
				largestValidLabel = numLabels-1

			return TopAccuracyOrientation(pred=pred_in_correct_shape, truth=truth_in_correct_shape, largestValidLabel=largestValidLabel, symmetric=symmetric)

		if labelType.startswith('LogNormal'):
			return TopAccuracyLogNormal(pred=pred_in_correct_shape, truth=truth_in_correct_shape, symmetric=symmetric)

		elif labelType.startswith('Normal'):
			return TopAccuracyNormal(pred=pred_in_correct_shape, truth=truth_in_correct_shape, symmetric=symmetric)

		elif labelType.startswith('Discrete'):
			#subType = labelType[len('Discrete'): ]
			if subType.startswith('2C'):
				return TopAccuracy2C(pred=pred_in_correct_shape, truth=truth_in_correct_shape, symmetric=symmetric)
			else:
				return TopAccuracyMultiC(pred=pred_in_correct_shape, truth=truth_in_correct_shape, subType=subType, symmetric=symmetric)
		else:
			print 'ERROR: unsupported label type in EvaluateAccuracy: ', labelType
			exit(1)


	accuracyList = []
	for res, out_prob, z, ratio in zip(self.responses, self.output_probList, zList, self.modelSpecs['topRatios'] ):

		labelName, labelType, subType = config.ParseResponse(res)

		## currently TopAccuracy only works when the dimension of each z is 3
		assert z.ndim == 3

		if self.mask_1d is not None:
            		paddingLens = self.mask_1d.shape[1] - T.sum(self.mask_1d, axis=1)
		else:
	    		paddingLens = T.zeros_like(z[:,0,0], dtype=np.int32)

		currentResponse = res
		topRatio = ratio

		##here we use scan to calculate accuracy for each protein
        	result, updates = theano.scan(fn=EvaluateAccuracy, outputs_info=None, sequences=[out_prob, z, paddingLens] )
		accuracy = T.mean(result, axis=0)
		accuracyList.append(accuracy)

	return T.stacklists(accuracyList)

## When forTrain is True, the model is used for training and validation; otherwise, the model is used for prediction
def BuildModel(modelSpecs, forTrain=True):
        rng = np.random.RandomState()

        ## x is for sequential features and y for matrix (or pairwise) features
        x = T.tensor3('x')
        y = T.tensor4('y')

        ## mask for x and y, respectively
        xmask = T.bmatrix('xmask')
        ymask = T.btensor3('ymask')

        xem = None
        ##if any( k in modelSpecs['seq2matrixMode'] for k in ('SeqOnly', 'Seq+SS') ):
        if config.EmbeddingUsed(modelSpecs):
                xem = T.tensor3('xem')

	## bounding box for crop of a big protein distance matrix. This box allows crop at any position.
	box = None
	if forTrain:
		box = T.ivector('boundingbox')

	## trainByRefLoss can be either 1 or -1. When this variable exists, we train the model using both reference loss and the loss of real data
	trainByRefLoss = None
	if forTrain and config.TrainByRefLoss(modelSpecs):
		trainByRefLoss = T.iscalar('trainByRefLoss')

        distancePredictor = ResNet4DistMatrix( rng, seqInput=x, matrixInput=y, mask_seq=xmask, mask_matrix=ymask, embedInput=xem, boundingbox=box, modelSpecs=modelSpecs )

        ## labelList is a list of label tensors, each having shape (batchSize, seqLen, seqLen) or (batchSize, seqLen, seqLen, valueDims[response] )
        labelList = []
        if forTrain:
                ## when this model is used for training. We need to define the label variable
		for response in modelSpecs['responses']:
			labelType = Response2LabelType(response)
			rValDims = GetResponseValueDims(response)

			if labelType.startswith('Discrete'):
				if rValDims > 1:
					## if one response is a vector, then we use a 4-d tensor
					## wtensor is for 16bit integer
					labelList.append( T.wtensor4('Tlabel4' + response ) )
				else:
					labelList.append( T.wtensor3('Tlabel4' + response ) )
			else:
				if rValDims > 1:
					labelList.append( T.tensor4('Tlabel4' + response ) )
				else:
					labelList.append( T.tensor3('Tlabel4' + response ) )

        ## weightList is a list of label weight tensors, each having shape (batchSize, seqLen, seqLen)
        weightList = []
        if len(labelList)>0 and config.UseSampleWeight(modelSpecs):
                weightList = [ T.tensor3('Tweight4'+response) for response in modelSpecs['responses'] ]

	## for prediction, both labelList and weightList are empty
	if forTrain:
        	return distancePredictor, x, y, xmask, ymask, xem, labelList, weightList, box, trainByRefLoss
	else:
        	return distancePredictor, x, y, xmask, ymask, xem

def TestResNet4DistMatrix():

    x = T.tensor3('x')
    y = T.tensor4('y')
    xmask = T.bmatrix('xmask')
    ymask = T.btensor3('ymask')
    selection = T.wtensor3('selection')

    import cPickle
    fh = open('seqDataset4HF.pkl')
    data = cPickle.load(fh)
    fh.close()

    distancePredictor = ResNet4DistMatrix( rng=np.random.RandomState(), seqInput=x, matrixInput=y, n_in_seq=data[0][0].shape[2], n_in_matrix=data[1][0].shape[3], n_hiddens_seq=[3, 5], n_hiddens_matrix=[2], hwsz_seq=4, hwsz_matrix=4, mask_seq=xmask, mask_matrix=ymask)

    """
    f = theano.function([x, y, xmask, ymask], distancePredictor.output_1d)
    g = theano.function([x, y, xmask, ymask], distancePredictor.output_2d)
    """

    dataLen = 300
    batchSize = 60
    a = np.random.uniform(0, 1, (batchSize, dataLen, 20)).astype(np.float32)
    b = np.random.uniform(0, 1, (batchSize, dataLen, dataLen, 3)).astype(np.float32)
    amask = np.zeros((batchSize, 0)).astype(np.int8)
    bmask = np.zeros((batchSize, 0, dataLen)).astype(np.int8)
    sel = np.ones((batchSize, dataLen, dataLen)).astype(np.int8)
    #print a
    #print b
    c = np.random.uniform(0, 3, (batchSize, dataLen, dataLen)).round().astype(np.int8)
    np.putmask(c, c>=2, 2)
    """
    c[0, 1, 13]=1
    c[0, 2, 15]=1
    c[0, 4, 16]=1

    c[0, 1, 27]=1
    c[0, 2, 28]=1
    c[0, 4, 29]=1

    c[1, 0, 13]=2
    c[1, 1, 15]=2
    c[1, 3, 16]=2

    c[2, 0, 23]=2
    c[2, 1, 25]=2
    c[2, 3, 26]=2
    """
    #sel = c

    #out1d = f(a, b, amask, bmask)
    #out2d = g(a, b, amask, bmask)

    #print out1d
    #print out2d

    z = T.btensor3('z')
    loss = distancePredictor.loss(z, selection)
    errs = distancePredictor.ErrorsByRange(z)
    accs = distancePredictor.TopAccuracyByRange(z)
    confM = distancePredictor.confusionMatrix(z)

    h = theano.function([x, y, xmask, ymask, selection, z], confM, on_unused_input='ignore')
    #l, e, accu = h(a, b, amask, bmask, sel, c)

    cms = []
    for i in np.arange(5):
    	cm = h(data[0][i], data[1][i], data[2][i], data[3][i], data[4][i], data[5][i])
        print cm
	cms.append(cm)

    sumofcms = np.sum(cms, axis=0)*1.

    for i in xrange(sumofcms.shape[0]):
        sumofcms[i] /= np.sum(sumofcms[i])

    confusions = sumofcms
    print confusions
    print np.sum(confusions[0])
    print np.sum(confusions[1])
    print np.sum(confusions[2])

    """
    params = contactPredictor.params
    cost = loss + 0.1 * contactPredictor.paramL2
    grads = [ T.grad(cost, p) for p in params ]
    fg = theano.function([x, y, xmask, ymask, selection, z], grads)
    g0 = fg(a, b, amask, bmask, sel, c)

    param_shapes = [ p.get_value().shape for p in params ]
    print g0
    """

if __name__ == "__main__":

    #TestOuterConcatenate() 

    #TestConv1D2Matrix()

    #TestConv2D4DistMatrix()

    TestResNet4DistMatrix()
    
    #TestMidpointFeature()
