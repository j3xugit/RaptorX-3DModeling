"""
multi-layer neural network for single-varibale or two-variable normal distribution.
"""

import os
import sys
import time
import math

import numpy

import theano
import theano.tensor as T

#from Optimizers import AdaGrad, AdaDelta, SGDMomentum, GD
from Adams import Adam
#from LogReg import LogisticRegression as LogReg

# start-snippet-1
class HiddenLayer(object):
    def __init__(self, rng, input, n_in, n_out, W=None, b=None, activation=T.tanh):
        """
        Typical hidden layer of a MLP: units are fully-connected and have
        sigmoidal activation function. Weight matrix W is of shape (n_in,n_out)
        and the bias vector b is of shape (n_out,).

        NOTE : The nonlinearity used here is tanh

        Hidden unit activation is given by: tanh(dot(input,W) + b)

        :type rng: numpy.random.RandomState
        :param rng: a random number generator used to initialize weights

        :type input: theano.tensor.dmatrix
        :param input: a symbolic tensor of shape (n_examples, n_in)

        :type n_in: int
        :param n_in: dimensionality of input
	self.n_in = n_in

        :type n_out: int
        :param n_out: number of hidden units
	self.n_out = n_out

        :type activation: theano.Op or function
        :param activation: Non linearity to be applied in the hidden
                           layer
        """
        self.input = input
        self.n_in = n_in
        self.n_out = n_out


        if W is None:
            	W_values = numpy.asarray( rng.uniform( low = -numpy.sqrt(6. / (n_in + n_out)), high = numpy.sqrt(6. / (n_in + n_out)), size=(n_in, n_out) ), dtype=theano.config.floatX )
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

"""
## x is a matrix with shape (batchSize, 2)
## this function returns T.prod(x, axis=1, keepdims=True)
## we reimplement this because theano has bugs with T.prod()
def MyProd(x):
        y = T.mul(x[:, 0], x[:, 1])
        return y.dimshuffle(0, 'x')
"""

class NN4PhiPsi(object):
    	"""A deep neural network for phi, psi prediction.
	 We model them using a distribution described in https://salilab.org/modeller/9.18/manual/node492.html .

    	A multi-layer feedforward artificial neural network that has one layer or more of hidden units and nonlinear activations.
	This neural network predicts the 5 parameters of the phi, psi distribution by maximum likelihood
    	"""

    	def __init__(self, rng, input=None, n_in=1, n_variables=2, n_out=5, n_hiddens=[], mymean=None, sigma_sqr_min=numpy.float32(0.0001), rho_abs_max=numpy.float32(0.99)):
        	"""
        	rng: a random number generator used to initialize weights
	
		input has shape (batchSize, n_in)
		n_in is the number of input features

		n_variables (=2) indicates the number of variables, i.e., phi and psi
		output has shape (batchSize, n_out)
		n_out (=5) is the number of parameters defining the distribution. Meanwhile, two for mean, two for variance and one for correlation

        	n_hidden: a list containing the number of hidden units at each hidden layer

		if you already have the mean, but just want to estimate variance and rho, then input a mean by mymean

		sigma_sqr_min is the minimum value of sigma square

        	"""

		assert (n_variables == 2)
                assert (n_out == 2 or n_out == 4 or n_out == 5)

		self.n_variables = n_variables
        	self.input = input
        	self.n_in = n_in
		self.n_out = n_out
        	self.n_hiddens = n_hiddens

        	self.params = []
        	self.paramL1 =0
        	self.paramL2 =0

        	self.hlayers = []
		self.layers = []

        	output_in_last_layer = input
        	n_out_in_last_layer = n_in

            	## add hidden layers
        	for i in xrange(len(n_hiddens)):
            		hiddenLayer = HiddenLayer( rng = rng, input = output_in_last_layer, n_in = n_out_in_last_layer, n_out = n_hiddens[i], activation = T.nnet.relu ) 
            		self.hlayers.append(hiddenLayer)
            		output_in_last_layer = hiddenLayer.output
            		n_out_in_last_layer = n_hiddens[i]

		self.layers = self.hlayers

		self.mean = None
		self.sigma_sqr = None
		self.corr = None

		self.params4var = []
		self.paramL14var = 0
		self.paramL24var = 0

		if mymean is not None:
			self.mean = mymean
		else:
			## calculate the mean
			uLayer = HiddenLayer( rng = rng, input = output_in_last_layer, n_in = n_out_in_last_layer, n_out = 2, activation = T.tanh )
			self.mean = uLayer.output * numpy.pi
			self.layers.append(uLayer)

		if n_out >= (2 * n_variables):
			##calculate sigma_sqr, which is bounded, so we use tanh here
			sigmaLayer = HiddenLayer( rng = rng, input = output_in_last_layer, n_in = n_out_in_last_layer, n_out = 2, activation = T.tanh )
			self.sigma_sqr = (1. + sigmaLayer.output)/2. * numpy.square(numpy.pi) + sigma_sqr_min
			self.layers.append(sigmaLayer)

			self.params4var += sigmaLayer.params
			self.paramL14var += sigmaLayer.paramL1
			self.paramL24var += sigmaLayer.paramL2

		if n_out == 5:
			##calculate correlation, need to make sure that correlation falls into [-1, 1]
			corrLayer = HiddenLayer( rng = rng, input = output_in_last_layer, n_in = n_out_in_last_layer, n_out = 1, activation = T.tanh )
			self.corr = corrLayer.output * rho_abs_max
			self.layers.append(corrLayer)

			self.params4var += corrLayer.params
			self.paramL14var += corrLayer.paramL1
			self.paramL24var += corrLayer.paramL2

		for layer in self.layers:
        		self.params += layer.params
			self.paramL1 += layer.paramL1
        		self.paramL2 += layer.paramL2

		self.y_pred = self.mean 

		outputList = [ self.mean ]
                if self.sigma_sqr is not None:
                        outputList.append(self.sigma_sqr)

                if self.corr is not None:
                        outputList.append(self.corr)

                self.output = T.concatenate(outputList, axis=1)


	## y has shape (batchSize, n_variables), sampleWeight has shape (batchSize, 1) instead of (batchSize,)
	## this function returns a scalar
	## useMeanOnly: only use mean to calculate NLL, ignore the variance and correlation
    	def NLL(self, y, useMeanOnly=False, sampleWeight=None):

		assert (y.ndim == 2)
		pi = numpy.pi

		err = y - self.mean
		cos_err = 1. -  T.cos(err)

		if useMeanOnly==True or (self.sigma_sqr is None):
			sig_sqr = T.ones_like(cos_err)
		else:
			sig_sqr = self.sigma_sqr

		e = T.sum(cos_err/sig_sqr, axis=1, keepdims=True )

		sig = T.sqrt(sig_sqr)
		sin_err = T.sin(err)
		f = T.prod(sin_err/sig, axis=1, keepdims=True)

		if useMeanOnly==True or (self.corr is None):
			rho = T.zeros_like(f)
		else:
			rho = self.corr

		g = e - T.mul(rho, f)
		rho_sqr = T.sqr(rho)
		h = g/(1 - rho_sqr )

		nll = h + numpy.log(2*pi) + T.sum(T.log(sig_sqr), axis=1, keepdims=True)/2. + T.log(1 - rho_sqr)/2.

		if sampleWeight is None:
			return T.mean(nll)
		return T.sum(T.mul(nll, sampleWeight) )/T.sum(sampleWeight)


	## y has shape (batchSize, n_variables), sampleWeight shall have shape (batchSize, 1) instead of (batchSize,)
	## this function returns a vector
    	def errors(self, y, sampleWeight=None):
		assert (y.ndim == 2)
		err1 = T.abs_( y - self.y_pred )
		err2 = numpy.float32(2 * numpy.pi) - err1
		err = T.minimum(err1, err2)

		if sampleWeight is None:
			return T.mean(err, axis=0 )
		else:
			assert (sampleWeight.ndim==2)
			weight = T.concatenate( [sampleWeight, sampleWeight], axis=1 )
			return T.sum(T.mul(err, weight), axis=0)/ T.sum(sampleWeight) 

	## y has shape (batchSize, n_variables), sampleWeight shall have shape (batchSize, 1) instead of (batchSize,)
	## useMeanOnly: only use mean to calculate the loss, ignore the variance and correlation
    	def loss(self, y, useMeanOnly=False, sampleWeight=None):
		if useMeanOnly == True:
			return self.NLL(y, useMeanOnly=True, sampleWeight=sampleWeight)
        	return self.NLL(y, sampleWeight=sampleWeight )



def testNN4PhiPsi(learning_rate=0.01, L1_reg=0.00, L2_reg=0.01, n_epochs=2000, n_hiddens=[100, 200, 200], trainData=None, testData=None):

    	## generate some random train and test data
	batchSize = 200000
	nFeatures = 50

    	trainX = numpy.random.uniform(0, 12, (batchSize, nFeatures)).astype(numpy.float32)
    	u1 = numpy.sum( trainX[:,:30], axis=1, keepdims=True)
    	u2 = numpy.sum( trainX[:,21:], axis=1, keepdims=True)

	numpy.putmask(u1, u1>180, u1-360)
	numpy.putmask(u2, u2>180, u2-360)
	u1 = u1/180.0 * numpy.pi
	u2 = u2/180.0 * numpy.pi

	trainY = (numpy.random.normal(0, 0.4, (batchSize, 2)) +  numpy.concatenate( (u1, u2), axis=1) ).astype(numpy.float32)
	numpy.putmask(trainY, trainY > numpy.pi, trainY - 2*numpy.pi)
	numpy.putmask(trainY, trainY < -numpy.pi, trainY + 2*numpy.pi)

	testBatchSize = 500
    	testX = numpy.random.uniform(0, 12, (testBatchSize, nFeatures)).astype(numpy.float32)
    	testu1 = numpy.sum(testX[:,:30], axis=1, keepdims=True)
    	testu2 = numpy.sum(testX[:,21:], axis=1, keepdims=True)
	testCorr = numpy.sum(testX[:,21:30], axis=1, keepdims=True)/numpy.sum(testX, axis=1, keepdims=True)

	numpy.putmask(testu1, testu1>180, testu1-360)
	numpy.putmask(testu2, testu2>180, testu2-360)
	testu1 = testu1/180.0 * numpy.pi
	testu2 = testu2/180.0 * numpy.pi
	testY = (numpy.random.normal(0, 0.4, (testBatchSize, 2)) + numpy.concatenate( (testu1, testu2), axis=1) ).astype(numpy.float32)
	numpy.putmask(testY, testY > numpy.pi, testY - 2*numpy.pi)
	numpy.putmask(testY, testY < -numpy.pi, testY + 2*numpy.pi)

    	######################
    	# BUILD ACTUAL MODEL #
    	######################
    	print('... building the model')

    	# allocate symbolic variables for the data

    	x = T.matrix('x')  # the input feature
    	y = T.matrix('y')  # the response

    	rng = numpy.random.RandomState()

    	regressor = NN4PhiPsi(rng, input=x, n_in=trainX.shape[1], n_hiddens=n_hiddens)

	loss = regressor.loss(y)
    	cost = loss + L1_reg * regressor.paramL1 + L2_reg * regressor.paramL2
	error = regressor.errors(y)

    	gparams = [T.grad(cost, param) for param in regressor.params]
    	param_shapes = [ param.shape.eval() for param in regressor.params ]
    	updates, others = Adam(regressor.params, gparams) 

    	train = theano.function( inputs=[x,y], outputs=[loss, error, regressor.paramL1, regressor.paramL2], updates=updates)
    	calculate = theano.function( inputs=[x], outputs=regressor.output )

    	#test = theano.function( inputs=[x,y], outputs=regressor.NLL(y))
    	test = theano.function( inputs=[x,y], outputs=regressor.errors(y))

    	step = 500
	numEpochs = 12
	for j in range(0, numEpochs):
    		for i in range(0,trainX.shape[0], step):
        		los, err, l1, l2 = train(trainX[i:i+step, :], trainY[i:i+step, :])
			if i%10000 == 0:
				print 'j=', j, ' i=', i, ' loss=', los, ' error=', err, ' L1norm=', l1, ' L2norm=', l2


	out = calculate(testX)
	print numpy.concatenate( (out, testCorr, testY), axis=1)
	print 'err=', test(testX, testY)
	corr = numpy.concatenate( (out[:,4:5], testCorr), axis=1)
	print numpy.corrcoef( numpy.transpose(corr) )

	import scipy
	print scipy.stats.mstats.spearmanr(corr[:,0], corr[:,1])

if __name__ == '__main__':
    	testNN4PhiPsi()
