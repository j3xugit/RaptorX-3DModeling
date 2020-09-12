"""
This script conducts 1D convolutional operation on a protein sequence, which is represented as a matrix of dimesion L * n_in
where L is the sequence length and n_in is the number of input features at each residue 

the input is a tensor with shape (numProteins, seqLen, n_in). The output has a shape (numProteins, seqLen, numOfOutFeatures)
two parameters: W for filter and b for bias vector
W has a shape (numOfOutFeatures, numOfInFeatures, 1, WindowSize) and b has a shape (numOfOutFeatures,)
"""
import pickle as cPickle
import gzip
import os
import sys
import time

import numpy
from numpy import random as rng

import theano
import theano.tensor as T

##this class conducts 1D convolution of one input
class Conv1DLayer(object):

    def __init__(self, rng, input, numOfInFeatures=0, numOfOutFeatures=0, halfWinSize=0, activation=T.nnet.relu, mask=None):

	## input is a tensor3, has shape (batchSize, seqLen, numOfInFeatures)
	self.input = input
	self.n_in = numOfInFeatures
	self.n_out = numOfOutFeatures
        windowSize = 2*halfWinSize + 1
	self.filter_size = windowSize


        # reshape input for conv2d
        # reshape input to shape (batchSize, numOfInFeatures, numRows, numColumns) where numColumns is equal to seqLen and numRows=1
        in4conv2D = input.dimshuffle(0,2,'x',1)
        inshape = (None, numOfInFeatures, 1, None)
        
        # initialize the filter
        # define the filter shape, which shall be (numOutFeatures, numInFeatures, 1, winSize)
        w_shp = (numOfOutFeatures, numOfInFeatures, 1, windowSize)
        W_values = numpy.asarray(
                rng.uniform(
                    low = -numpy.sqrt(6. / (numOfInFeatures*windowSize + numOfOutFeatures)),
                    high = numpy.sqrt(6. / (numOfInFeatures*windowSize + numOfOutFeatures)),
                    size = w_shp
                ),
                dtype=theano.config.floatX
            )
        if activation == T.nnet.sigmoid:
            W_values *= 4

        self.W = theano.shared(value=W_values, name='conv1d_W', borrow=True)
        
	b_shp = (numOfOutFeatures,)
	self.b = theano.shared(numpy.asarray( rng.uniform(low=-.0, high=.0, size=b_shp), dtype=input.dtype), name ='conv1d_b', borrow=True)


	# computes the convolution of input with filters in w
        # conv_out and conv_out_bias have shape (batch_size, numOfOutFeatures, 1, numColumns)
	conv_out = T.nnet.conv2d(in4conv2D, self.W, input_shape=inshape, filter_shape=w_shp, border_mode='half')
	if activation is not None:
	    conv_out_bias = activation(conv_out + self.b.dimshuffle('x',0,'x','x'))
	else:
	    conv_out_bias = (conv_out + self.b.dimshuffle('x',0,'x','x'))

        ## out2 now has shape (batchSize, numColumns, numOfOutFeatures)
        out2 = conv_out_bias.dimshuffle(0, 3, 1, 2)[:, :, :, 0]

        if mask is not None:
            ## since we did zero padding at some left columns of the input tensor,
            ## we need to reset these positions to 0 again after convolution to avoid introducing noise
	    ## mask has shape (batchSize, #positions_to_be_masked)

	    ##take the subtensor of out2 that needs modification
            out2_sub = out2[:, :mask.shape[1], :]
            mask_new = mask.dimshuffle(0, 1,'x')
            self.output = T.set_subtensor(out2[:, :mask.shape[1], :], T.mul(out2_sub, mask_new) )
        else:
            self.output = out2

        # parameters of the model
	self.params=[self.W, self.b]

        self.paramL1 = abs(self.W).sum() + abs(self.b).sum()
        self.paramL2 = (self.W**2).sum() + (self.b**2).sum()

def testConv1DLayer():

    rng = numpy.random.RandomState()

    input = T.tensor3('input')

    #windowSize = 3
    n_in = 4
    n_hiddens = [10,10,5]
    #convR = Conv1DR(rng, input, n_in, n_hiddens, windowSize/2)
    convLayer = Conv1DLayer(rng, input, n_in, 5, halfWinSize=1)
    
    #f = theano.function([input],convR.output)    
    #f = theano.function([input],[convLayer.output, convLayer.out2, convLayer.convout, convLayer.out3])    
    f = theano.function([input], convLayer.output)    

    numOfProfiles=6
    seqLen = 10
    profile = numpy.random.uniform(0,1, (numOfProfiles, seqLen,n_in))
    
    out = f(profile)
    print out.shape
    print out

if __name__ == '__main__':
    testConv1DLayer()
