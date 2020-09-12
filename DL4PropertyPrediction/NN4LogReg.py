"""
multi-layer neural network for classification using Theano.
"""

import os
import sys
import time

import numpy

import theano
import theano.tensor as T

#from Optimizers import AdaGrad, AdaDelta, SGDMomentum, GD
from Adams import Adam

class HiddenLayer(object):
    def __init__(self, rng, input, n_in, n_out, W=None, b=None, activation=T.tanh):
        """
        rng: a random number generator used to initialize weights

        input: a symbolic tensor of shape (batchSize, n_in)
        n_in: dimensionality of input
        n_out: number of hidden units

        activation: Non linearity to be applied in the hidden layer
        """
        self.input = input
        self.n_in = n_in
        self.n_out = n_out

        if W is None:
            W_values = numpy.asarray( rng.uniform( low = -numpy.sqrt(6./(n_in + n_out)), high = numpy.sqrt(6./(n_in + n_out)), size=(n_in, n_out)), dtype=theano.config.floatX )
            if activation == T.nnet.sigmoid:
                W_values *= 4

            W = theano.shared(value=W_values, name='HL_W', borrow=True)

        if b is None:
            b_values = numpy.zeros((n_out,), dtype=theano.config.floatX)
            b = theano.shared(value=b_values, name='HL_b', borrow=True)

        self.W = W
        self.b = b

        lin_output = T.dot(input, self.W) + self.b
        self.output = ( lin_output if activation is None else activation(lin_output) )

        # parameters of the model
        self.params = [self.W, self.b]
        self.paramL1 = abs(self.W).sum() + abs(self.b).sum()
        self.paramL2 = (self.W**2).sum() + (self.b**2).sum()

## a simple logistic regression for classification
class LogRegLayer(object):

    def __init__(self, rng, input, n_in, n_out):
        """ 
        input: symbolic variable that describes the input of the architecture (one minibatch). It has shape (batchSize, n_in)
        n_in: number of input units, the dimension of the space in which the datapoints lie
        n_out: number of output units, the dimension of the space in which the labels lie

        """
        self.n_in = n_in
        self.n_out = n_out

        # initialize with 0 the weights W as a matrix of shape (n_in, n_out)
        value_bound = numpy.sqrt(6. /(n_in + n_out))
        W_values = numpy.asarray(rng.uniform( low = - value_bound, high = value_bound, size=(n_in, n_out) ), dtype=theano.config.floatX)
        self.W = theano.shared (value = W_values, name ='LogReg_W', borrow=True)

        # initialize the baises b as a vector of n_out 0s
        self.b = theano.shared( value=numpy.zeros( (n_out,), dtype=theano.config.floatX ), name='LogReg_b', borrow=True)

        self.pre_act = T.dot(input, self.W) + self.b
        self.p_y_given_x = T.nnet.softmax(self.pre_act)
        self.y_pred = T.argmax(self.p_y_given_x, axis=1)
        self.output = self.p_y_given_x

       # parameters of the model
        self.params = [self.W, self.b]
        self.paramL1 = abs(self.W).sum() + abs(self.b).sum()
        self.paramL2 = (self.W**2).sum() + (self.b**2).sum()

    def NLL(self, y, sampleWeight=None):
        ###Return the mean of the negative log-likelihood of the prediction of this model under a given target distribution.

        if sampleWeight is not None:
            return -T.sum(T.mul(sampleWeight, T.log(self.p_y_given_x)[T.arange(y.shape[0]), y] ) )/T.sum(sampleWeight)
        else:
            return -T.mean(T.log(self.p_y_given_x)[T.arange(y.shape[0]), y])

    def errors(self, y, sampleWeight=None):
        ###Return the 0-1 error rate in a minibatch y: a vector of true labels

        # check if y has same dimension of y_pred
        if y.ndim != self.y_pred.ndim:
            raise TypeError(
                'y should have the same shape as self.y_pred',
                ('y', y.type, 'y_pred', self.y_pred.type)
            )
        # check if y is of the correct datatype
        if y.dtype.startswith('int'):
            # the T.neq operator returns a vector of 0s and 1s, where 1 represents a mistake in prediction
            if sampleWeight is not None:
                return T.sum( T.mul(sampleWeight, T.neq(self.y_pred, y) ) ) * 1./T.sum(sampleWeight)
            else:
                return T.mean(T.neq(self.y_pred, y))
        else:
            raise NotImplementedError()


### A neural network Logistic Regression for Classification
class NN4LogReg(object):

    def __init__(self, rng, input, n_in, n_out, n_hiddens=[], activation=T.nnet.relu):
        """Initialize the parameters for the multilayer perceptron

        rng: a random number generator used to initialize weights
	
	input has shape (batchSize, n_in)
	n_in is the number of input features
	n_out is the number of classes (or labels)

        n_hidden: a tuple defining the number of hidden units at each hidden layer
	activation: the nonlinear function for the hidden layers

        """
        self.input = input
        self.n_in = n_in
        self.n_hiddens = n_hiddens

        self.hlayers = []
	self.layers = []

        output_in_last_layer = input
        n_out_in_last_layer = n_in

        for i in xrange(len(n_hiddens)):

            hiddenLayer = HiddenLayer( rng=rng, input=output_in_last_layer, n_in=n_out_in_last_layer, n_out=n_hiddens[i], activation=activation) 
            self.hlayers.append(hiddenLayer)

            output_in_last_layer = hiddenLayer.output
            n_out_in_last_layer = n_hiddens[i]


	## add the final logistic regression layer
        linLayer = LogRegLayer(rng, output_in_last_layer, n_out_in_last_layer, n_out)
	self.linLayer = linLayer
	self.layers = self.hlayers + [ self.linLayer ]

	self.pre_act = linLayer.pre_act
	self.p_y_given_x = linLayer.p_y_given_x

	## here we make self.y_pred have shape (batchSize, 1) instead of (batchSize, )
	self.y_pred = linLayer.y_pred.dimshuffle(0, 'x')

	## self.output has shape (batchSize, n_out)
        self.output = self.p_y_given_x
	self.n_out = n_out

        self.paramL1 =0
        self.paramL2 =0
        self.params = []

	for layer in self.layers:
		self.paramL1 += layer.paramL1
        	self.paramL2 += layer.paramL2
        	self.params += layer.params

    ## Both y and sampleWeight shall have shape (batchSize, 1) instead of (batchSize,)
    ## this function returns a scalar
    ## useMeanOnly here shall always be set to False, it is used only for placeholder
    def NLL(self, y, useMeanOnly=False, sampleWeight=None):
	assert (y.ndim == 2)
	##convert to 1d vector
	y0 = y[:,0]

	if sampleWeight is None:
		return self.linLayer.NLL(y0)

	assert (sampleWeight.ndim == 2)
	w = sampleWeight[:, 0]
		
        return self.linLayer.NLL(y0, w)

    ## Both y and sampleWeight shall have shape (batchSize, 1) instead of (batchSize,)
    ## this function returns a vector
    def errors(self, y, sampleWeight=None):
	assert (y.ndim == 2)
	err = T.neq(self.y_pred, y)
	if sampleWeight is None:
		return T.mean(err, axis=0)

	assert (sampleWeight.ndim == 2)
	return T.sum( T.mul(err, sampleWeight), axis=0)/T.sum(sampleWeight)

    def loss(self, y, useMeanOnly=False, sampleWeight=None):
        return self.NLL(y, useMeanOnly, sampleWeight)



def testNN4LogReg(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=2000, n_hiddens=[200, 200], trainData=None, testData=None):

    	## generate some random train and test data
	batchSize = 200000
	nFeatures = 30
    	trainX = numpy.random.uniform(0, 1, (batchSize, nFeatures)).astype(numpy.float32)
    	trainXsum = numpy.sum(trainX**2, axis=1, keepdims=True)
    	trainY = numpy.zeros((batchSize, 1), dtype=numpy.int32 )
    	numpy.putmask(trainY, trainXsum>5, 1)
    	numpy.putmask(trainY, trainXsum>10, 2)
    	numpy.putmask(trainY, trainXsum>15, 3)

	testBatchSize = 50
    	testX = numpy.random.uniform(0, 1, (testBatchSize, nFeatures)).astype(numpy.float32)
    	testXsum = numpy.sum(testX**2, axis=1, keepdims=True)
    	testY = numpy.zeros((testBatchSize, 1), dtype=numpy.int32 )
    	numpy.putmask(testY, testXsum>5, 1)
    	numpy.putmask(testY, testXsum>10, 2)
    	numpy.putmask(testY, testXsum>15, 3)

    	######################
    	# BUILD ACTUAL MODEL #
    	######################
    	print('... building the model')


    	x = T.matrix('x')  # the data is presented as rasterized images
    	y = T.imatrix('y')  # the labels 

    	rng = numpy.random.RandomState()

    	regressor = NN4LogReg(rng, input=x, n_in=trainX.shape[1], n_hiddens=n_hiddens, n_out=4, activation=T.nnet.relu)
	loss = regressor.loss(y)
	error = regressor.errors(y)
	cost = loss + L1_reg * regressor.paramL1 + L2_reg * regressor.paramL2

    	gparams = [T.grad(cost, param) for param in regressor.params]
    	param_shapes = [ param.shape.eval() for param in regressor.params ]
    	updates, others = Adam(regressor.params, gparams) 

    	train = theano.function( inputs=[x,y], outputs=[loss, error, regressor.paramL1, regressor.paramL2], updates=updates)
    	test = theano.function( inputs=[x,y], outputs=error)
	calculate = theano.function( inputs=[x], outputs=regressor.output )

        step = 200
        numEpochs = 30
        for j in range(0, numEpochs):
                results = []
                for i in range(0, trainX.shape[0], step):
                        los, err, l1, l2 = train(trainX[i:i+step, :], trainY[i:i+step, :])
                        results.append( los )
                        if i%5000 == 0:
                                print 'i=', i, ' loss=', los, ' error=', err, ' L1norm=', l1, ' L2norm=', l2

                print 'j=', j, ' avgLos, avgErr=', numpy.mean(results, axis=0)


        out = calculate(testX)
	print numpy.around( numpy.concatenate( (out, testY.reshape(testBatchSize, 1) ), axis=1), 2)

    	print 'err=', test(testX, testY)

if __name__ == '__main__':
    testNN4LogReg()
