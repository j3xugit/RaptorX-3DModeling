import os
import sys
import numpy as np
import theano
import theano.tensor as T

# this function returns the mean of each channel of input x. The returned value is a 1d vector.
def AvgPool(x, mask=None):
        ## This function calculates the mean of each channel by excluding the impact of masked positions in the input
        ## x is the input and shall have shape (batchSize, n_in, nCols) or (batchSize, n_in, nRows, nCols)
        ## mask shall have shape (batchSize, #cols_to_be_masked) or (batchSize, #rows_to_be_masked, nCols)
	## n_in: the number of channels in x
        ## Assume that all the entries in x to be masked have already been set to 0

        if mask is not None:
                if x.ndim == 4:
                        x_sum = x.sum(axis=[0,2,3])

                        ## x_num is the number of effective elements in x excluding the positions to be masked
                        x_num = mask.sum(dtype=theano.config.floatX) * 2 + T.cast( x.shape[0] * (x.shape[2] - mask.shape[1]) * (x.shape[3] - mask.shape[1]), theano.config.floatX)

                        ## remove repeated count. the 1s in mask[:, :mask.shape[1], :mask.shape[1]] are counted twice, so we need to do adjustment
                        x_num = x_num - mask[:, :mask.shape[1], :mask.shape[1]].sum(dtype=theano.config.floatX)
                        x_mean = x_sum / x_num

                elif x.ndim == 3:
                        x_sum = x.sum(axis=[0,2])
                        x_num = mask.sum(dtype=theano.config.floatX) + T.cast( x.shape[0] * (x.shape[2] - mask.shape[1]), theano.config.floatX)
                        x_mean = x_sum / x_num
                else:
                        print 'ERROR: the ndim of input for AvgPool can only be 3 or 4!'
                        exit(1)

                return x_mean

        ##the below code will be excuted when mask is None
        if x.ndim == 4:
                x_mean = x.mean(axis=[0,2,3])
        elif x.ndim == 3:
                x_mean = x.mean(axis=[0,2])
        else:
                print 'ERROR: the ndim of input for AvgPool shall be 3 or 4!'
                exit(1)

        return x_mean

# this function returns the max of each channel of input x. The returned value is a 1d vector.
# x is the input and shall have shape (batchSize, n_in, nCols) or (batchSize, n_in, nRows, nCols)
# Assume that all the entries in x to be masked have already been set to 0
def MaxPool(x):
	## we calculate the max of each channel for each instance and then average them
        if x.ndim == 4:
                x_max = x.max(axis=[2,3]).mean(axis=[0])
        elif x.ndim == 3:
                x_max = x.max(axis=[2]).mean(axis=[0])
        else:
                print 'ERROR: the ndim of input for MaxPool shall be 3 or 4!'
                exit(1)

        return x_max

## input has shape (n_in,) and output has the same shape
class FullConnectionLayer(object):
	def __init__(self, input, n_in, activation=T.nnet.sigmoid, W=None):
		self.input = input

		if W is None:
            		W_values = np.random.uniform(low=-np.sqrt(6./(n_in+n_in)), high=np.sqrt(6./(n_in+n_in)), size=(n_in, n_in)).astype(theano.config.floatX)
            		if activation == T.nnet.sigmoid:
                		W_values *= 4
            		W = theano.shared(value=W_values, name='FC_W', borrow=True)
        	self.W = W

		input2 = input.dimshuffle('x', 0)
        	output = T.dot(input2, self.W)
		if activation is not None:
        		output2 = activation(output)
		self.output = output2.dimshuffle(1, 0)[:,0]

        	# parameters of the model
        	self.params = [self.W]
        	self.paramL1 = abs(self.W).sum() 
        	self.paramL2 = (self.W**2).sum()

## avg and maxhave shape (n_in,) and output has the same shape
class FullConnectionLayer2(object):
	def __init__(self, avg, max, n_in, activation=T.nnet.sigmoid, W=None):
		self.input = [avg, max]

		if W is None:
            		W_values = np.random.uniform(low=-np.sqrt(6./(n_in+n_in)), high=np.sqrt(6./(n_in+n_in)), size=(n_in, n_in)).astype(theano.config.floatX)
            		if activation == T.nnet.sigmoid:
                		W_values *= 4
            		W = theano.shared(value=W_values, name='FC_W', borrow=True)
        	self.W = W

		input_avg = avg.dimshuffle('x', 0)
		input_max = max.dimshuffle('x', 0)
        	output_avg = T.dot(input_avg, self.W)
        	output_max = T.dot(input_max, self.W)

		if activation is not None:
        		output_avg = activation(output_avg)
        		output_max = activation(output_max)

		self.output = (output_avg + output_max).dimshuffle(1, 0)[:,0] /np.float32(2)

        	# parameters of the model
        	self.params = [self.W]
        	self.paramL1 = abs(self.W).sum() 
        	self.paramL2 = (self.W**2).sum()

#The input is an 1d vector with shape (n_in,) and The output has the same shape
class SimpleConvLayer(object):
	def __init__(self, input, activation=T.nnet.sigmoid, halfWinSize=2, W=None):

                self.input = input
                # reshape input to shape (1, 1, nRows=n_in, nCols=1)
                in4conv2D = input.dimshuffle('x', 'x', 0, 'x')

		if W is None:
                	windowSize = 2*halfWinSize + 1
                	w_shp = (1, 1, windowSize, 1)
                	W_values = np.asarray(np.random.uniform(low=-np.sqrt(6./windowSize), high=np.sqrt(6./windowSize), size=w_shp), dtype=theano.config.floatX)
                	if activation == T.nnet.sigmoid:
                        	W_values *= 4
		else:
			w_shp = W.shape

                self.W = theano.shared(value=W_values, name='SimpleConv1d_W', borrow=True)

                # conv_out has shape (1, 1, n_in, 1)
                conv_out = T.nnet.conv2d(in4conv2D, self.W, filter_shape=w_shp, border_mode='half')
                if activation is not None:
                        conv_out = activation(conv_out)

                ## output has shape (n_in,)
                self.output = conv_out.dimshuffle(2, 3, 0, 1)[:, 0, 0, 0]

        	# parameters of the model
        	self.params = [self.W]
        	self.paramL1 = abs(self.W).sum() 
        	self.paramL2 = (self.W**2).sum()

class SimpleConvLayer2(object):
	def __init__(self, avg, max, activation=T.nnet.sigmoid, halfWinSize=2, W=None):

                self.input = [avg, max]
                # reshape input to shape (1, 1, nRows=n_in, nCols=1)
                input_avg = avg.dimshuffle('x', 'x', 0, 'x')
                input_max = max.dimshuffle('x', 'x', 0, 'x')

		if W is None:
                	windowSize = 2*halfWinSize + 1
                	w_shp = (1, 1, windowSize, 1)
                	W_values = np.asarray(np.random.uniform(low=-np.sqrt(6./windowSize), high=np.sqrt(6./windowSize), size=w_shp), dtype=theano.config.floatX)
                	if activation == T.nnet.sigmoid:
                        	W_values *= 4
		else:
			w_shp = W.shape

                self.W = theano.shared(value=W_values, name='SimpleConv1d_W', borrow=True)

                # conv_out has shape (1, 1, n_in, 1)
                out_avg = T.nnet.conv2d(input_avg, self.W, filter_shape=w_shp, border_mode='half')
                out_max = T.nnet.conv2d(input_max, self.W, filter_shape=w_shp, border_mode='half')

                if activation is not None:
                        out_avg = activation(out_avg)
                        out_max = activation(out_max)

                ## output has shape (n_in,)
                self.output = (out_avg + out_max).dimshuffle(2, 3, 0, 1)[:, 0, 0, 0] /np.float32(2)

        	# parameters of the model
        	self.params = [self.W]
        	self.paramL1 = abs(self.W).sum() 
        	self.paramL2 = (self.W**2).sum()

## input has shape (batchSize, n_in, nRows, nCols)
## output shall have the same shape
class AttentionLayer(object):
	def __init__(self, input, n_in, mask=None, UseAvg=True, UseMax=False, UseFC=True):
		if (not UseAvg) and (not UseMax):
			print 'ERROR: at least AvgPool or MaxPool shall be used'
			exit(1)

		self.input = input

		if UseAvg:	
			input_avg = AvgPool(input, mask)
		if UseMax:
			input_max = MaxPool(input)

		if UseFC:
			if UseAvg and UseMax: 
				fcLayer = FullConnectionLayer2(input_avg, input_max, n_in)
			elif UseAvg:
				fcLayer = FullConnectionLayer(input_avg, n_in)
			else:
				fcLayer = FullConnectionLayer(input_max, n_in)
		else:
			if UseAvg and UseMax: 
				fcLayer = SimpleConvLayer2(input_avg, input_max)
			elif UseAvg:
				fcLayer = SimpleConvLayer(input_avg)
			else:
				fcLayer = SimpleConvLayer(input_max)

		attn = fcLayer.output
		self.attn = attn

		if input.ndim==4:
			attn = attn.dimshuffle(['x', 0, 'x', 'x'])
		else:
			attn = attn.dimshuffle(['x', 0, 'x'])

		self.output = T.mul(input, attn)

		self.params = fcLayer.params
		self.paramL1 = fcLayer.paramL1
		self.paramL2 = fcLayer.paramL2


## for test only
if __name__ == "__main__":
	x = T.tensor4('x')
	y = AvgPool(x)
	z = MaxPool(x)
	avgPooling = theano.function([x], y)
	maxPooling = theano.function([x], z)

	batchSize = 2
	nFeatures = 4
	nRows = 5
	nCols = 5
	a = np.random.uniform(0, 1, (batchSize, nFeatures, nRows, nCols)).astype(theano.config.floatX)
	print a.shape
	print a

	b = avgPooling(a)
	print 'avgPool=', b

	c = maxPooling(a)
	print 'maxPool=', c

	attnLayer = AttentionLayer(x, nFeatures, UseFC=False)
	attnFunc = theano.function([x], [attnLayer.output, attnLayer.attn])
	d, w = attnFunc(a)
	print 'w_attn_conv= ', w
	print 'result='
	print d

	attnLayer = AttentionLayer(x, nFeatures)
	attnFunc = theano.function([x], [attnLayer.output, attnLayer.attn])
	e, w = attnFunc(a)
	print 'w_attn_fc=', w
	print 'result='
	print e
