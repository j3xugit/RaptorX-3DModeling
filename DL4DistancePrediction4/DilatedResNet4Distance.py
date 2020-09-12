import numpy as np
import theano
import theano.tensor as T

import config
from AttentionLayer import AttentionLayer

##a new implementation of 1D convolution 
class ResConv1DLayer(object):

    	def __init__(self, rng, input, n_in=0, n_out=0, halfWinSize=0, dilation=1, activation=T.nnet.relu, mask=None):
       		""" Initialize the parameters
        	The input has shape (batchSize, n_in, seqLen)
		The output will have shape (batchSize, n_out, seqLen)
        	"""

        	self.input = input
        	self.n_in = n_in
        	self.n_out = n_out
		self.halfWinSize = halfWinSize

        	windowSize = 2*halfWinSize + 1
        	self.filter_size = windowSize
		self.dilation = (1, dilation)

        	# reshape input to shape (batchSize, n_in, nRows=1, nCols=seqLen) 
        	in4conv2D = input.dimshuffle(0, 1, 'x', 2)

        	# initialize the filter, i.e., convoltion kernel
        	w_shp = (n_out, n_in, 1, windowSize)

		if activation == T.nnet.relu:
           		W_values = np.asarray(rng.normal(scale=np.sqrt(2. / (n_in*windowSize + n_out)), size = w_shp), dtype=theano.config.floatX )
		else:
            		W_values = np.asarray(
                		rng.uniform(low = - np.sqrt(6. / (n_in*windowSize + n_out)), high = np.sqrt(6. / (n_in*windowSize + n_out)), size = w_shp),
                		dtype=theano.config.floatX
            		)
            	if activation == T.nnet.sigmoid:
            		W_values *= 4

        	self.W = theano.shared(value=W_values, name='ResConv1d_W', borrow=True)

        	b_shp = (n_out,)
        	self.b = theano.shared(np.asarray( rng.uniform(low=-.0, high=.0, size=b_shp), dtype=input.dtype), name ='ResConv1d_b', borrow=True)

        	# conv_out and conv_out_bias have shape (batch_size, n_out, 1, nCols)
		if dilation > 1:
			conv_out = T.nnet.conv2d(in4conv2D, self.W, filter_shape=w_shp, border_mode='half', filter_dilation=(1, dilation) )
		else:
        		conv_out = T.nnet.conv2d(in4conv2D, self.W, filter_shape=w_shp, border_mode='half')
        	if activation is not None:
            		conv_out_bias = activation(conv_out + self.b.dimshuffle('x', 0, 'x', 'x'))
        	else:
            		conv_out_bias = (conv_out + self.b.dimshuffle('x', 0, 'x', 'x'))

		## out2 has shape (batchSize, n_out, nCols)
        	out2 = conv_out_bias.dimshuffle(0, 1, 3, 2)[:, :, :, 0]

        	if mask is not None:
            		## Zeros are padded at left columns of the input tensor, we need to reset these positions to 0 again after convolution to avoid noise
            		## mask has shape (batchSize, #positions_to_be_masked)

            		##take the subtensor of out2 that needs modification
            		out2_sub = out2[:, :, :mask.shape[1] ]
            		mask_new = mask.dimshuffle(0, 'x', 1)
            		self.output = T.set_subtensor(out2_sub, T.mul(out2_sub, mask_new) )
        	else:
            		self.output = out2

		##self.output has shape (batchSize, n_out, nCols)

        	# parameters of the model
        	self.params=[self.W, self.b]

        	self.paramL1 = abs(self.W).sum() + abs(self.b).sum()
        	self.paramL2 = (self.W**2).sum() + (self.b**2).sum()


## 2D convolution on a matrix, no pooling is done here
class ResConv2DLayer:
    	def __init__(self, rng, input, n_in=0, n_out=0, halfWinSize=0, dilation=1, activation=T.nnet.relu, mask=None):
        	## input has shape (batchSize, n_in, nRows, nCols) where n_in is the number of input features at each element
        	## here we assume that input[0], input[1], ...., input[batchSize -1 ] are aligned at the right bottom
        	## nrows and ncols are unknown
        	## mask has shape (batchSize, #rows_to_be_masked, nCols), mask is a binary matrix with 0 indicates padding positions and 1 real data position

        	##halfWinSize is the half window size of the filter shape
        	##Mask is used to reduce noise introduced by padding. all the matrices in a batch shall be aligned at the right bottom corner

		##output shall have shape (batchSize, n_out, nRows, nCols)

       		self.input = input
        	self.n_in = n_in
        	self.halfWinSize = halfWinSize
        	self.mask = mask

        	##the window size of our filters is always odd
        	wSize = 2*halfWinSize + 1
		self.filter_size = wSize
		self.dilation = (dilation, dilation)

        	assert n_out >= 0
        	assert n_in >= 0

        	W_shape = (n_out, n_in, wSize, wSize)
	
		if activation == T.nnet.relu:
            		W_values = np.asarray( rng.normal( scale = np.sqrt(2. / (n_in * wSize * wSize + n_out )), size = W_shape), dtype = theano.config.floatX)
		else:
            		W_values = np.asarray(
            			rng.uniform( low = - np.sqrt(6. / (n_in * wSize * wSize + n_out )), high =  np.sqrt(6. / (n_in * wSize * wSize + n_out )), size = W_shape),
            			dtype = theano.config.floatX
           		)
            	if activation == T.nnet.sigmoid:
            		W_values *= 4

        	W = theano.shared (value = W_values, name = 'ResConv2dlayer_W', borrow = True)

        	b_shape = (n_out, )
        	b_values = np.asarray (rng.uniform(low = -.0, high =.0, size = b_shape), dtype=theano.config.floatX )
        	b = theano.shared (value = b_values, name = 'ResConv2dlayer_b', borrow = True)

        	## conv2d_out and out2 have shape (batchSize, n_out, nRows, nCols)
		if dilation > 1:
        		conv2d_out = T.nnet.conv2d(input, W, filter_shape=W_shape, border_mode='half', filter_dilation=(dilation, dilation) )
		else:
        		conv2d_out = T.nnet.conv2d(input, W, filter_shape=W_shape, border_mode='half')
        	if activation is not None:
           		out2 = activation(conv2d_out + b.dimshuffle('x', 0, 'x', 'x'))
        	else:
            		out2 = conv2d_out + b.dimshuffle('x', 0, 'x', 'x')
       
        	if mask is not None:
            		## mask has shape (batchSize, #rows_to_be_masked, nCols)

            		## a subtensor of out2 along the horiz direction
            		out2_sub_horiz = out2[:, :, :mask.shape[1], :]
            		mask_horiz = mask.dimshuffle(0, 'x', 1, 2)
            		out3 = T.set_subtensor(out2_sub_horiz, T.mul(out2_sub_horiz, mask_horiz) )

            		## a subtensor of out3 along the vertical direction
            		out3_sub_vertical = out3[:, :, :, :mask.shape[1] ]
            		mask_vertical = mask.dimshuffle(0, 'x', 2, 1)
            		out4 = T.set_subtensor(out3_sub_vertical, T.mul(out3_sub_vertical, mask_vertical) )

            		self.output = out4

        	else:
            		self.output = out2

        	self.params = [W, b]
        	self.paramL1 = (abs(W).sum() +  abs(b).sum())
        	self.paramL2 = ((W**2).sum() + (b**2).sum())

        	self.n_out = n_out

## 2D convolution with attention on a matrix, no pooling is done here
class AttnConv2DLayer:
    	#def __init__(self, rng, input, n_in=0, n_out=0, halfWinSize=0, dilation=1, activation=T.nnet.relu, mask=None, attnOptions=None):
    	def __init__(self, rng, input, n_in=0, n_out=0, halfWinSize=0, dilation=1, activation=T.nnet.relu, mask=None):
        	## input has shape (batchSize, n_in, nRows, nCols) where n_in is the number of input features at each element
        	## here we assume that input[0], input[1], ...., input[batchSize -1 ] are aligned at the right bottom
        	## nrows and ncols are unknown
        	## mask has shape (batchSize, #rows_to_be_masked, nCols), mask is a binary matrix with 0 indicates padding positions and 1 real data position
		## attnOptions is a dict() for attention mechanism

        	##halfWinSize is the half window size of the filter shape
        	##Mask is used to reduce noise introduced by padding. all the matrices in a batch shall be aligned at the right bottom corner

		##output shall have shape (batchSize, n_out, nRows, nCols)

       		self.input = input
        	self.n_in = n_in
        	self.halfWinSize = halfWinSize
        	self.mask = mask
		#self.attnOptions = attnOptions

        	##the window size of our filters is always odd
        	wSize = 2*halfWinSize + 1
		self.filter_size = wSize
		self.dilation = (dilation, dilation)

        	assert n_out >= 0
        	assert n_in >= 0

        	W_shape = (n_out, n_in, wSize, wSize)
	
		if activation == T.nnet.relu:
            		W_values = np.asarray(rng.normal(scale=np.sqrt(2./(n_in * wSize * wSize + n_out )), size=W_shape), dtype=theano.config.floatX)
		else:
            		W_values = np.asarray(rng.uniform(low=-np.sqrt(6./(n_in * wSize * wSize + n_out )), high=np.sqrt(6./(n_in * wSize * wSize + n_out )), size=W_shape), dtype = theano.config.floatX)

            		if activation == T.nnet.sigmoid:
            			W_values *= 4

        	W = theano.shared (value = W_values, name = 'AttnConv2dlayer_W', borrow=True)

        	b_shape = (n_out, )
        	b_values = np.asarray (rng.uniform(low=-.0, high=.0, size=b_shape), dtype=theano.config.floatX )
        	b = theano.shared (value=b_values, name='AttnConv2dlayer_b', borrow = True)

        	## conv2d_out and out2 have shape (batchSize, n_out, nRows, nCols)
		if dilation > 1:
        		conv2d_out = T.nnet.conv2d(input, W, filter_shape=W_shape, border_mode='half', filter_dilation=(dilation, dilation) )
		else:
        		conv2d_out = T.nnet.conv2d(input, W, filter_shape=W_shape, border_mode='half')

            	conv_out = conv2d_out + b.dimshuffle('x', 0, 'x', 'x')

		#attnLayer = AttentionLayer(conv_out, n_out, mask=mask, UseAvg=attnOptions['UseAvg'], UseMax=attnOptions['UseMax'], UseFC=attnOptions['UseFC'])
		attnLayer = AttentionLayer(conv_out, n_out, mask=mask, UseAvg=True, UseMax=True, UseFC=True)
		attnOut = attnLayer.output

        	if activation is not None:
           		out2 = activation(attnOut)
		else:
			out2 = attnOut
       
        	if mask is not None:
            		## mask has shape (batchSize, #rows_to_be_masked, nCols)

            		## a subtensor of out2 along the horiz direction
            		out2_sub_horiz = out2[:, :, :mask.shape[1], :]
            		mask_horiz = mask.dimshuffle(0, 'x', 1, 2)
            		out3 = T.set_subtensor(out2_sub_horiz, T.mul(out2_sub_horiz, mask_horiz) )

            		## a subtensor of out3 along the vertical direction
            		out3_sub_vertical = out3[:, :, :, :mask.shape[1] ]
            		mask_vertical = mask.dimshuffle(0, 'x', 2, 1)
            		out4 = T.set_subtensor(out3_sub_vertical, T.mul(out3_sub_vertical, mask_vertical) )

            		self.output = out4
        	else:
            		self.output = out2

        	self.params = [W, b] + attnLayer.params
        	self.paramL1 = abs(W).sum() +  abs(b).sum()  + attnLayer.paramL1
        	self.paramL2 = (W**2).sum() + (b**2).sum() + attnLayer.paramL2

        	self.n_out = n_out

def batch_norm(x, n_in, mask=None, eps=1e-6):
	## This batch_norm calculates the mean and standard deviation by excluding the impact of masked positions in the input
    	## x is the input and shall have shape (batchSize, n_in, nCols) or (batchSize, n_in, nRows, nCols)
    	## mask shall have shape (batchSize, #cols_to_be_masked) or (batchSize, #rows_to_be_masked, nCols)
    	## Assume that all the entries in x to be masked have already been set to 0

    	gamma = theano.shared(np.asarray(np.ones((n_in,)), dtype=theano.config.floatX), borrow=True)
    	bias = theano.shared(np.asarray(np.zeros((n_in,)), dtype=theano.config.floatX), borrow=True)
    	if mask is not None:
       		if x.ndim == 4:
	    		x_sum = x.sum(axis=[0,2,3], keepdims=True)

	    		## x_num is the number of effective elements in x by excluding the positions to be masked
	    		x_num = mask.sum(dtype=theano.config.floatX) * 2 + T.cast( x.shape[0] * (x.shape[2] - mask.shape[1]) * (x.shape[3] - mask.shape[1]), theano.config.floatX)

	    		## remove repeated count, the 1s in mask[:, :mask.shape[1], :mask.shape[1]] are counted twice, so we need to do adjustment
	    		x_num = x_num - mask[:, :mask.shape[1], :mask.shape[1]].sum(dtype=theano.config.floatX)
	    		x_mean = x_sum / x_num

			##calculate the average of x^2
	    		x_sqr = T.sqr(x)
	    		x_sq_sum = x_sqr.sum(axis=[0,2,3], keepdims=True)
	    		x_sq_mean = x_sq_sum / x_num
			
			##calculate the standard devition, eps makes sure that sqrt() will not produce NaN
	    		x_std = (x_sq_mean - T.sqr(x_mean) + eps ).sqrt()

            		y_temp = T.nnet.bn.batch_normalization(x, gamma[None,:,None,None], bias[None,:,None,None], x_mean, x_std, mode = "low_mem")	

	    		## reset the zero-padded entries to 0 along the horizonal direction
            		y_temp_sub_horiz = y_temp[:, :, :mask.shape[1], :]
            		mask_horiz = mask.dimshuffle(0, 'x', 1, 2)
            		y_temp2 = T.set_subtensor(y_temp_sub_horiz, T.mul(y_temp_sub_horiz, mask_horiz) )

            		## reset the zero-padded entries to 0 along the vertical direction
            		y_temp_sub_vertical = y_temp2[:, :, :, :mask.shape[1] ]
            		mask_vertical = mask.dimshuffle(0, 'x', 2, 1)
            		y = T.set_subtensor(y_temp_sub_vertical, T.mul(y_temp_sub_vertical, mask_vertical) )

		elif x.ndim == 3:
	    		x_sum = x.sum(axis=[0,2], keepdims=True)
	    		x_num = mask.sum(dtype=theano.config.floatX) + T.cast( x.shape[0] * (x.shape[2] - mask.shape[1]), theano.config.floatX)
	    		x_mean = x_sum / x_num

            		x_sq_sum = T.sqr(x).sum(axis=[0,2], keepdims=True)
	    		x_sq_mean = x_sq_sum / x_num
	    		x_std = (x_sq_mean - T.sqr(x_mean) + eps ).sqrt()

	    		y_intermediate = T.nnet.bn.batch_normalization(x, gamma[None,:,None], bias[None,:,None], x_mean, x_std, mode = "low_mem")

            		## reset the zero-padded entries to 0 again
	    		y_intermediate_sub = y_intermediate[:, :, :mask.shape[1] ]
            		mask_new = mask.dimshuffle(0, 'x', 1)
            		y = T.set_subtensor(y_intermediate_sub, T.mul(y_intermediate_sub, mask_new) )

		else:
	    		print 'the ndim of input for batch_norm can only be 3 or 4!'
	    		sys.exit(-1)

        	return y, [gamma, bias]

	##the below code will be excuted when mask is None
    	if x.ndim == 4:
       		x_mean = x.mean(axis=[0,2,3], keepdims=True)
        	x_std = (x.var(axis=[0,2,3], keepdims=True) + eps ).sqrt()
        	y = T.nnet.bn.batch_normalization(x, gamma[None,:,None,None], bias[None,:,None,None], x_mean, x_std, mode = "low_mem")
    	elif x.ndim == 3:
       		x_mean = x.mean(axis=[0,2], keepdims=True)
        	x_std = (x.var(axis=[0,2], keepdims=True) + eps ).sqrt()
        	y = T.nnet.bn.batch_normalization(x, gamma[None,:,None], bias[None,:,None], x_mean, x_std, mode = "low_mem")
    	else:
		print 'the ndim of input for batch_norm can only be 3 or 4!'
		sys.exit(-1)

    	return y, [gamma, bias]

class BatchNormLayer:
	def __init__(self, input, n_in, mask=None):
       		self.input = input
		self.n_in = n_in
		bnout, bnparams = batch_norm(input, n_in, mask=mask)
		self.output = bnout
		self.params = bnparams
		self.n_out = n_in

		self.paramL1 = abs(bnparams[0]).sum() + abs(bnparams[1]).sum()
		self.paramL2 = (bnparams[0]**2).sum() + (bnparams[1]**2).sum()

def TestBatchNorm3D():
    n_in = 20
    maxSeqLen = 200
    minSeqLen = 190
    batchSize = 3

    matrices = []

    container = np.zeros( (batchSize, n_in, maxSeqLen), dtype=np.float32 )
    masks = np.zeros( (batchSize, maxSeqLen - minSeqLen), dtype=np.float32 )

    for i in xrange(batchSize):
        nRows = np.random.randint(minSeqLen, maxSeqLen)
        nCols = nRows
	a = np.random.uniform(0, 2, (n_in, nCols))
	matrices.append(a)
	container[i, :, maxSeqLen - nCols: ] = a

	masks[i, maxSeqLen-nCols: ].fill(1)

    ## calculate the mean of all the matrices

    b = np.zeros((n_in, 1))
    b2 = np.zeros((n_in, 1))
    numElements = 0

    for a in matrices:
	b = b + np.sum(a, axis=1, keepdims=True)
	b2 = b2 + np.sum( np.square(a), axis=1, keepdims=True)
	numElements = numElements + a.shape[1]

    m = b/numElements
    m2 = b2/numElements
    std = np.sqrt(m2 - m*m)

    newmatrices = []
    for a in matrices:
	b = (a - m ) /std
	newmatrices.append(b)

    x = T.tensor3('x')
    mask = T.matrix('mask')
    y, [gamma, bias] = batch_norm(x, n_in, mask)
    f = theano.function([x, mask], y)
    c = f(container, masks)

    for i in xrange(batchSize):
	a = newmatrices[i]
	c_sub = c[i][:, c[i].shape[1]-a.shape[1]: ]
	c[i][:, c[i].shape[1]-a.shape[1]: ] = c_sub - a
	print np.absolute(c[i]).max()

def TestBatchNorm4D():
    n_in = 20
    maxSeqLen = 200
    minSeqLen = 190
    batchSize = 5

    matrices = []

    container = np.zeros( (batchSize, n_in, maxSeqLen, maxSeqLen), dtype=np.float32 )
    masks = np.zeros( (batchSize, maxSeqLen  - minSeqLen, maxSeqLen), dtype=np.float32 )

    for i in xrange(batchSize):
        nRows = np.random.randint(minSeqLen, maxSeqLen)
        nCols = nRows
	a = np.random.uniform(0, 6, (n_in, nRows, nCols))
	matrices.append(a)
	container[i, :, maxSeqLen - nRows:, maxSeqLen - nCols: ] = a
	masks[i, maxSeqLen-nRows:, maxSeqLen - nRows: ].fill(1)

    ## calculate the mean of all the matrices

    b = np.zeros((n_in, 1, 1))
    b2 = np.zeros((n_in, 1, 1))
    numElements = 0

    for a in matrices:
	b = b + np.sum(a, axis=(1, 2), keepdims=True)
	b2 = b2 + np.sum( np.square(a), axis=(1, 2), keepdims=True)
	numElements = numElements + a.shape[1]*a.shape[2]

    m = b/numElements
    m2 = b2/numElements
    std = np.sqrt(m2 - m*m + 1e-6)

    newmatrices = []
    for a in matrices:
	b = (a - m ) /std
	newmatrices.append(b)

    x = T.tensor4('x')
    mask = T.tensor3('mask')
    y, [gamma, bias] = batch_norm(x, n_in, mask)
    f = theano.function([x, mask], y)
    c = f(container, masks)

    for i in xrange(batchSize):
	a = newmatrices[i]
	#print a
	c_sub = c[i][:, c[i].shape[1]-a.shape[1]:, c[i].shape[2]-a.shape[2]: ]
	#print c[i]
	c[i][:, c[i].shape[1]-a.shape[1]:, c[i].shape[2]-a.shape[2]: ] = c_sub - a
	print np.absolute(c[i]).max()
	#print c[i]

class BottleneckBlock:

    def __init__(self, rng, input, n_in, halfWinSize=1, mask=None, n_out=None, n_bottleneck=None, activation=T.tanh, dim_inc_method='partial_projection', batchNorm=False):
	## please make sure that the input has shape (batchSize, n_in, nRows, nCols) or (batchSize, n_in, nCols)

	##if n_out is not None, then it shall be no smaller than n_in
        if n_out is not None:
            assert n_out >= n_in
            self.n_out = n_out
        else:
            self.n_out = n_in

	if n_bottleneck is not None:
	    assert n_bottleneck < n_in
	    self.n_bottleneck = n_bottleneck
	else:
	    self.n_bottleneck = n_in / 2

        self.n_in = n_in
        self.input = input
        self.mask = mask
        self.halfWinSize = halfWinSize

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	elif input.ndim == 4:
	    ConvLayer = ResConv2DLayer
	else:
	    print 'the ndim of input can only be 3 or 4'
	    sys.exit(-1)

        if batchNorm:
	    l1_pre = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_bottleneck, halfWinSize=0, mask=mask, activation=activation)
	    l1 = BatchNormLayer(l1_pre.output, l1_pre.n_out)
	    l2_pre = ConvLayer(rng, input=l1.output, n_in=self.n_bottleneck, n_out=self.n_bottleneck, halfWinSize=halfWinSize, mask=mask, activation=activation)
	    l2 = BatchNormLayer(l2_pre.output, l2_pre.n_out)
	    l3_pre = ConvLayer(rng, input=l2.output, n_in=self.n_bottleneck, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    l3 = BatchNormLayer(l3_pre.output, l3_pre.n_out)
	    self.layers = [l1_pre, l1, l2_pre, l2, l3_pre, l3]
	else:
            ## if needed, we do dimension increase at the first layer
	    l1 = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_bottleneck, halfWinSize=0, mask=mask, activation=activation)
	    l2 = ConvLayer(rng, input=l1.output, n_in=self.n_bottleneck, n_out=self.n_bottleneck, halfWinSize=halfWinSize, mask=mask, activation=activation)
	    l3 = ConvLayer(rng, input=l2.output, n_in=self.n_bottleneck, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    self.layers = [l1, l2, l3]

	## intermediate has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	intermediate = l3.output

	if dim_inc_method == 'full_projection':
	    ## we do 1*1 convolution here without any nonlinear transformation
	    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    intermediate = intermediate + linlayer.output
	    self.layers.append(linlayer)
	    #print 'projection is True'

	elif dim_inc_method == 'identity':
	    if input.ndim == 3:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :], input)
	    elif input.ndim == 4:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :, :], input)
	    else:
	    	print 'the ndim of input can only be 3 or 4'
	    	sys.exit(-1)
	    print 'only projection is supported'
	    sys.exit(-1)

	elif dim_inc_method == 'partial_projection':
	    if self.n_out == n_in:
		intermediate = intermediate + input
	    else:
	        if batchNorm:
		    linlayer_pre = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		    linlayer = BatchNormLayer(linlayer_pre.output, linlayer_pre.n_out)
	            self.layers += [linlayer_pre, linlayer]	
	        else:
		    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		    self.layers.append(linlayer)

		intermediate = intermediate + T.concatenate([input, linlayer.output], axis=1)
	else:
	    print 'unsupported dimension increase method: ', dim_inc_method
	    exit(1)

	if activation is not None:
	    self.output = activation(intermediate)
	else:
	    self.output = intermediate

	## self.output has shape (batchSize, n_out, nRows, nCols)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0

	for layer in self.layers:
	    self.params += layer.params
	    self.paramL1 += layer.paramL1
	    self.paramL2 += layer.paramL2

class ResBlockV2:

    def __init__(self, rng, input, n_in, halfWinSize=0, mask=None, n_out=None, activation=T.nnet.relu, dim_inc_method='partial_projection', batchNorm=False, dropout=False):
	## please make sure that the input has shape (batchSize, n_in, nRows, nCols) or (batchSize, n_in, nCols)

	##if n_out is not None, then it shall be no smaller than n_in
        if n_out is not None:
            assert n_out >= n_in
            self.n_out = n_out
        else:
            self.n_out = n_in

        self.n_in = n_in
        self.input = input
        self.mask = mask
        self.halfWinSize = halfWinSize

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	elif input.ndim == 4:
	    ConvLayer = ResConv2DLayer
	else:
	    print 'the ndim of input can only be 3 or 4'
	    exit(1)

        if batchNorm:
	    bnlayer1 = BatchNormLayer(input, n_in, mask=mask)
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    bnlayer2 = BatchNormLayer(l1.output, l1.n_out, mask=mask)
	    input2 = activation(bnlayer2.output)

	    ## add dropout code here
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    self.layers = [bnlayer1, l1, bnlayer2, l2]
	    #self.layers = [ l1, bnlayer2, l2]
	else:
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    input2 = activation(l1.output)

	    ## add dropout code here
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    self.layers = [l1, l2]

	## intermediate has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	intermediate = l2.output

	if dim_inc_method == 'full_projection':
	    ## we do 1*1 convolution here without any nonlinear transformation
	    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    intermediate = intermediate + linlayer.output
	    self.layers.append(linlayer)
	    #print 'projection is True'

	elif dim_inc_method == 'identity':
	    if input.ndim == 3:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :], input)
	    elif input.ndim == 4:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :, :], input)
	    else:
	    	print 'the ndim of input can only be 3 or 4'
	    	exit(1)
	    print 'only projection is supported'
	    sys.exit(-1)

	elif dim_inc_method == 'partial_projection':
	    if self.n_out == n_in:
		intermediate = intermediate + input
	    else:
		linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		self.layers.append(linlayer)

		intermediate = intermediate + T.concatenate([input, linlayer.output], axis=1)
	else:
	    print 'unsupported dimension increase method: ', dim_inc_method
	    exit(1)

	self.output = intermediate

	## self.output has shape (batchSize, n_out, nRows, nCols)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0

	for layer in self.layers:
	    self.params += layer.params
	    self.paramL1 += layer.paramL1
	    self.paramL2 += layer.paramL2

## ResBlockV23 is almost same as ResBlockV2 except that the former has removed unused batchNormLayers and accordingly model parameters
class ResBlockV23:

    def __init__(self, rng, input, n_in, halfWinSize=0, mask=None, n_out=None, activation=T.nnet.relu, dim_inc_method='partial_projection', batchNorm=False, dropout=False):
	## please make sure that the input has shape (batchSize, n_in, nRows, nCols) or (batchSize, n_in, nCols)

	##if n_out is not None, then it shall be no smaller than n_in
        if n_out is not None:
            assert n_out >= n_in
            self.n_out = n_out
        else:
            self.n_out = n_in

        self.n_in = n_in
        self.input = input
        self.mask = mask
        self.halfWinSize = halfWinSize

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	elif input.ndim == 4:
	    ConvLayer = ResConv2DLayer
	else:
	    print 'the ndim of input can only be 3 or 4'
	    exit(1)

        if batchNorm:
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    bnlayer2 = BatchNormLayer(l1.output, l1.n_out, mask=mask)
	    input2 = activation(bnlayer2.output)

	    ## add dropout code here
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    ##self.layers = [bnlayer1, l1, bnlayer2, l2]
	    self.layers = [ l1, bnlayer2, l2]
	else:
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    input2 = activation(l1.output)

	    ## add dropout code here
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    self.layers = [l1, l2]

	## intermediate has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	intermediate = l2.output

	if dim_inc_method == 'full_projection':
	    ## we do 1*1 convolution here without any nonlinear transformation
	    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    intermediate = intermediate + linlayer.output
	    self.layers.append(linlayer)
	    #print 'projection is True'

	elif dim_inc_method == 'identity':
	    if input.ndim == 3:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :], input)
	    elif input.ndim == 4:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :, :], input)
	    else:
	    	print 'the ndim of input can only be 3 or 4'
	    	exit(1)
	    print 'only projection is supported'
	    exit(1)

	elif dim_inc_method == 'partial_projection':
	    if self.n_out == n_in:
		intermediate = intermediate + input
	    else:
		linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		self.layers.append(linlayer)

		intermediate = intermediate + T.concatenate([input, linlayer.output], axis=1)
	else:
	    print 'unsupported dimension increase method: ', dim_inc_method
	    exit(1)

	self.output = intermediate

	## self.output has shape (batchSize, n_out, nRows, nCols)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0

	for layer in self.layers:
	    self.params += layer.params
	    self.paramL1 += layer.paramL1
	    self.paramL2 += layer.paramL2

## ResBlockV22 has two batch normalization layers in each residual block while ResBlockV2 and ResBlockV23 have only one
class ResBlockV22:

    def __init__(self, rng, input, n_in, halfWinSize=0, mask=None, n_out=None, activation=T.nnet.relu, dim_inc_method='partial_projection', batchNorm=False, dropout=False):
	## please make sure that the input has shape (batchSize, n_in, nRows, nCols) or (batchSize, n_in, nCols)

	##if n_out is not None, then it shall be no smaller than n_in
        if n_out is not None:
            assert n_out >= n_in
            self.n_out = n_out
        else:
            self.n_out = n_in

        self.n_in = n_in
        self.input = input
        self.mask = mask
        self.halfWinSize = halfWinSize

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	elif input.ndim == 4:
	    ConvLayer = ResConv2DLayer
	else:
	    print 'the ndim of input can only be 3 or 4'
	    exit(1)

        if batchNorm:
	    bnlayer1 = BatchNormLayer(input, n_in, mask=mask)
	    input1 = activation(bnlayer1.output)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    bnlayer2 = BatchNormLayer(l1.output, l1.n_out, mask=mask)
	    input2 = activation(bnlayer2.output)

	    ## add dropout code here
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    self.layers = [bnlayer1, l1, bnlayer2, l2]
	else:
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    input2 = activation(l1.output)

	    ## add dropout code here
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)

	    self.layers = [l1, l2]

	## intermediate has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	intermediate = l2.output

	if dim_inc_method == 'full_projection':
	    ## we do 1*1 convolution here without any nonlinear transformation
	    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    intermediate = intermediate + linlayer.output
	    self.layers.append(linlayer)
	    #print 'projection is True'

	elif dim_inc_method == 'identity':
	    if input.ndim == 3:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :], input)
	    elif input.ndim == 4:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :, :], input)
	    else:
	    	print 'the ndim of input can only be 3 or 4'
	    	sys.exit(-1)
	    print 'only projection is supported'
	    exit(1)

	elif dim_inc_method == 'partial_projection':
	    if self.n_out == n_in:
		intermediate = intermediate + input
	    else:
		linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		self.layers.append(linlayer)

		intermediate = intermediate + T.concatenate([input, linlayer.output], axis=1)
	else:
	    print 'unsupported dimension increase method: ', dim_inc_method
	    exit(1)

	self.output = intermediate

	## self.output has shape (batchSize, n_out, nRows, nCols)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0

	for layer in self.layers:
	    self.params += layer.params
	    self.paramL1 += layer.paramL1
	    self.paramL2 += layer.paramL2
    
class ResBlockV1:

    def __init__(self, rng, input, n_in, halfWinSize=0, mask=None, n_out=None, activation=T.tanh, dim_inc_method='partial_projection', batchNorm=False):
	## please make sure that the input has shape (batchSize, n_in, nRows, nCols) or (batchSize, n_in, nCols)

	##if n_out is not None, then it shall be no smaller than n_in
        if n_out is not None:
            assert n_out >= n_in
            self.n_out = n_out
        else:
            self.n_out = n_in

        self.n_in = n_in
        self.input = input
        self.mask = mask
        self.halfWinSize = halfWinSize

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	elif input.ndim == 4:
	    ConvLayer = ResConv2DLayer
	else:
	    print 'the ndim of input can only be 3 or 4'
	    exit(1)

        if batchNorm:
	    l1_pre = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=activation)
	    l1 = BatchNormLayer(l1_pre.output, l1_pre.n_out)
	    l2_pre = ConvLayer(rng, input=l1.output, n_in=self.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)
	    l2 = BatchNormLayer(l2_pre.output, l2_pre.n_out)
	    self.layers = [l1_pre, l1, l2_pre, l2]
	else:
            ## if needed, we do dimension increase at the first layer
	    l1 = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=activation)
	    l2 = ConvLayer(rng, input=l1.output, n_in=self.n_out, n_out=self.n_out, halfWinSize=halfWinSize, mask=mask, activation=None)
	    self.layers = [l1, l2]

	## intermediate has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	intermediate = l2.output

	if dim_inc_method == 'full_projection':
	    ## we do 1*1 convolution here without any nonlinear transformation
	    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    intermediate = intermediate + linlayer.output
	    self.layers.append(linlayer)
	    #print 'projection is True'

	elif dim_inc_method == 'identity':
	    if input.ndim == 3:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :], input)
	    elif input.ndim == 4:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :, :], input)
	    else:
	    	print 'the ndim of input can only be 3 or 4'
	    	sys.exit(-1)
	    print 'only projection is supported'
	    exit(1)

	elif dim_inc_method == 'partial_projection':
	    if self.n_out == n_in:
		intermediate = intermediate + input
	    else:
	        if batchNorm:
		    linlayer_pre = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		    linlayer = BatchNormLayer(linlayer_pre.output, linlayer_pre.n_out)
	            self.layers += [linlayer_pre, linlayer]	
	        else:
		    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		    self.layers.append(linlayer)

		intermediate = intermediate + T.concatenate([input, linlayer.output], axis=1)
	else:
	    print 'unsupported dimension increase method: ', dim_inc_method
	    exit(1)

	if activation is not None:
	    self.output = activation(intermediate)
	else:
	    self.output = intermediate

	## self.output has shape (batchSize, n_out, nRows, nCols)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0

	for layer in self.layers:
	    self.params += layer.params
	    self.paramL1 += layer.paramL1
	    self.paramL2 += layer.paramL2

class DilatedResBlock:
    #def __init__(self, rng, input, n_in, halfWinSize=0, dilation=1, mask=None, n_out=None, activation=T.nnet.relu, dim_inc_method='partial_projection', batchNorm=False, dropout=False, modelSpecs=None):
    def __init__(self, rng, input, n_in, halfWinSize=0, dilation=1, mask=None, n_out=None, activation=T.nnet.relu, dim_inc_method='partial_projection', modelSpecs=None):
	## please make sure that the input has shape (batchSize, n_in, nRows, nCols) or (batchSize, n_in, nCols)

	##if n_out is not None, then it shall be no smaller than n_in
        if n_out is not None:
            assert n_out >= n_in
            self.n_out = n_out
        else:
            self.n_out = n_in

        self.n_in = n_in
        self.input = input
        self.mask = mask
        self.halfWinSize = halfWinSize
	self.dilation = dilation

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	elif input.ndim == 4:
		attnFlag = config.ParseAttentionMode(modelSpecs)
		if attnFlag is not None:
			"""
			bUseAvg, bUseMax, bUseFC = attnFlag
			attnOptions = dict()
			attnOptions['UseAvg'] = bUseAvg
			attnOptions['UseMax'] = bUseMax
			attnOptions['UseFC'] = bUseFC
			"""
			ConvLayer = AttnConv2DLayer
		else:
	    		ConvLayer = ResConv2DLayer
	else:
	    print 'ERROR: the ndim of input shall be 3 or 4'
	    exit(1)

	batchNorm = modelSpecs['batchNorm']
        if batchNorm:
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, dilation=dilation, mask=mask, activation=None)

	    bnlayer2 = BatchNormLayer(l1.output, l1.n_out, mask=mask)
	    input2 = activation(bnlayer2.output)

	    ## add dropout code here if needed
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, dilation=dilation, mask=mask, activation=None)

	    self.layers = [ l1, bnlayer2, l2]
	else:
	    input1 = activation(input)
	    l1 = ConvLayer(rng, input=input1, n_in=n_in, n_out=self.n_out, halfWinSize=halfWinSize, dilation=dilation, mask=mask, activation=None)
	    input2 = activation(l1.output)
	    l2 = ConvLayer(rng, input=input2, n_in=l1.n_out, n_out=self.n_out, halfWinSize=halfWinSize, dilation=dilation, mask=mask, activation=None)

	    self.layers = [l1, l2]

	## intermediate has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	intermediate = l2.output

	"""
	## add attention here
	attnFlag = config.ParseAttentionMode(modelSpecs)
	if attnFlag is not None:
		bUseAvg, bUseMax, bUseFC = attnFlag
		attnLayer = AttentionLayer(intermediate, self.n_out, mask, UseAvg=bUseAvg, UseMax=bUseMax, UseFC=bUseFC)
		intermediate = attnLayer.output
		self.layers.append(attnLayer)
	"""

	if dim_inc_method == 'full_projection':
	    ## we do 1*1 convolution here without any nonlinear transformation
	    linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out, halfWinSize=0, mask=mask, activation=None)
	    intermediate = intermediate + linlayer.output
	    self.layers.append(linlayer)

	elif dim_inc_method == 'identity':
	    if input.ndim == 3:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :], input)
	    elif input.ndim == 4:
	    	intermediate = T.inc_subtensor(intermediate[:, :n_in, :, :], input)
	    else:
	    	print 'ERROR: the ndim of input shall be 3 or 4'
	    	exit(1)
	    print 'ERROR: only projection is supported'
	    exit(1)

	elif dim_inc_method == 'partial_projection':
	    if self.n_out == n_in:
		intermediate = intermediate + input
	    else:
	    	## we do 1*1 convolution here without any nonlinear transformation
		linlayer = ConvLayer(rng, input=input, n_in=n_in, n_out=self.n_out - n_in, halfWinSize=0, mask=mask, activation=None)
		self.layers.append(linlayer)

		intermediate = intermediate + T.concatenate([input, linlayer.output], axis=1)
	else:
	    print 'ERROR: unsupported dimension increase method: ', dim_inc_method
	    exit(1)

	self.output = intermediate

	## self.output has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0

	for layer in self.layers:
	    self.params += layer.params
	    self.paramL1 += layer.paramL1
	    self.paramL2 += layer.paramL2

class DilatedResNet:
    #def __init__(self, rng, input, n_in, halfWinSize=[0], dilation=[1], mask=None, n_hiddens=None, n_repeats=None, activation=T.nnet.relu, dim_inc_method='partial_projection', batchNorm=True, version='DilatedResNet', modelSpecs=None):
    def __init__(self, rng, input, n_in, halfWinSize=[0], dilation=[1], mask=None, n_hiddens=None, n_repeats=None, activation=T.nnet.relu, dim_inc_method='partial_projection', modelSpecs=None):
	"""
	This DilatedResNet consists of the following components:
	1) a start layer with input as input, output has n_hiddens[0] features
	2) in total there are len(n_repeats) stacks
	3) each stack has 1 + n_repeats[i] blocks 
	4) the first block of each stack has n_hiddens[i-1] or n_in input features and n_hiddens[i] output features
	5) the other blocks of each stack has #inputfeatures = #outputfeatures = n_hiddens[i]

	input has shape (batchSize, nRows, nCols, n_in) or (batchSize, nCols, n_in)
	output has shape (batchSize, nRows, nCols, n_hiddens[-1]) or (batchSize, nCols, n_hiddens[-1])
	n_hiddens needs to be an increasing sequence
	n_repeats shall be non-negative
	"""

	assert n_hiddens is not None
	assert n_repeats is not None
	assert len(n_hiddens) > 0
	assert len(n_hiddens) == len(n_repeats)
	assert len(halfWinSize) == len(n_hiddens)
	assert len(dilation) == len(n_hiddens)

	version = modelSpecs['network']
	assert version.startswith('DilatedResNet')

	if input.ndim == 3:
	    ConvLayer = ResConv1DLayer
	    input2 = input.dimshuffle(0, 2, 1)
	elif input.ndim == 4:
	    ConvLayer = ResConv2DLayer
	    input2 = input.dimshuffle(0, 3, 1, 2)
	else:
	    print 'the ndim of input can only be 3 or 4!'
	    exit(1)

	ResBlock = DilatedResBlock

	blocks = []
	startLayer = ConvLayer(rng, input=input2, n_in=n_in, n_out=n_hiddens[0], halfWinSize=halfWinSize[0], dilation=dilation[0], mask=mask, activation=activation)
	blocks.append(startLayer)

	## repeat a block n_repeats[i] times
	curr_block = startLayer
	for j in range(n_repeats[0]):
	    assert (curr_block.n_out == n_hiddens[0])
            #bnflag =(batchNorm and j==0)
	    #new_block = ResBlock(rng, input=curr_block.output, n_in=n_hiddens[0], mask=mask, halfWinSize=halfWinSize[0], dilation=dilation[0], activation=activation, dim_inc_method=dim_inc_method, batchNorm=batchNorm)
	    new_block = ResBlock(rng, input=curr_block.output, n_in=n_hiddens[0], mask=mask, halfWinSize=halfWinSize[0], dilation=dilation[0], activation=activation, dim_inc_method=dim_inc_method, modelSpecs=modelSpecs)
	    blocks.append(new_block)
	    curr_block = new_block

	for i in range(1, len(n_hiddens)):
	    ## the start block is in charge of dimension increase
	    assert (curr_block.n_out == n_hiddens[i-1])

	    #new_block = ResBlock(rng, input=curr_block.output, n_in=n_hiddens[i-1], n_out=n_hiddens[i], mask=mask, halfWinSize=halfWinSize[i], dilation=dilation[i], activation=activation, dim_inc_method=dim_inc_method, batchNorm=batchNorm)
	    new_block = ResBlock(rng, input=curr_block.output, n_in=n_hiddens[i-1], n_out=n_hiddens[i], mask=mask, halfWinSize=halfWinSize[i], dilation=dilation[i], activation=activation, dim_inc_method=dim_inc_method, modelSpecs=modelSpecs)
	    blocks.append(new_block)
	    curr_block = new_block

	    ## repeat a block n_repeats[i] times
	    for j in range(n_repeats[i]):
		assert (curr_block.n_out == n_hiddens[i])
		#new_block = ResBlock(rng, input=curr_block.output, n_in=n_hiddens[i], mask=mask, halfWinSize=halfWinSize[i], dilation=dilation[i], activation=activation, dim_inc_method=dim_inc_method, batchNorm=batchNorm)
		new_block = ResBlock(rng, input=curr_block.output, n_in=n_hiddens[i], mask=mask, halfWinSize=halfWinSize[i], dilation=dilation[i], activation=activation, dim_inc_method=dim_inc_method, modelSpecs=modelSpecs)
		blocks.append(new_block)
	        curr_block = new_block

	out2 = curr_block.output
	self.n_out = curr_block.n_out

	## out2 has shape (batchSize, n_out, nRows, nCols) or (batchSize, n_out, nCols)
	## change the output shape back to (batchSize, nRows, nCols, n_out) or (batchSize, nCols, n_out)
	if input.ndim == 3:
	    self.output = out2.dimshuffle(0, 2, 1)
	elif input.ndim == 4:
	    self.output = out2.dimshuffle(0, 2, 3, 1)
	else:
	    print 'the ndim of input can only be 3 or 4!'
	    exit(1)

	self.params = []
	self.paramL1 = 0
	self.paramL2 = 0
	for block in blocks:
	    self.params += block.params
	    self.paramL1 += block.paramL1
	    self.paramL2 += block.paramL2

	self.layers = blocks

def TestConvLayers():

    from Conv1d import Conv1DLayer
    from Conv2d import Conv2DLayer

    rng = np.random.RandomState()
    n_in = 3
    n_out = 4
    seqLen = 15
    nRows = 20
    nCols = 20

    m = T.tensor3('m')
    x = T.tensor4('x')
    x2 = x.dimshuffle(0, 3, 1, 2)
    #l1 = ResConv2DLayer(rng, input=x2, n_in=n_in, n_out=n_out, halfWinSize=1, activation=T.nnet.relu, mask=m )
    l1 = ResConv2DLayer(rng, input=x2, n_in=n_in, n_out=n_out, halfWinSize=1, activation=T.nnet.relu )
    y = l1.output.dimshuffle(0, 2, 3, 1)
    #f = theano.function([x], y)

    cost = T.mean( (y)**2 )
    params = l1.params
    gparams = T.grad(cost, params)
    #h = theano.function([x, m], gparams)
    h = theano.function([x], gparams)

    bSize = 2
    a = np.random.uniform(0, 1, (bSize, nRows, nCols, n_in))
    m_value = np.ones( (bSize, 1, nCols) )
    #gs = h(a, m_value)
    gs = h(a)

    for g in gs:
	print g.shape

    """
    l2 = Conv2DLayer(rng, input=x, n_in=n_in, n_out=n_out, halfWinSize=1, activation=T.nnet.relu)
    z = l2.output
    g = theano.function([x], z)

    for p, q in zip(l1.params, l2.params):
	q.set_value(p.get_value() )

    bSize = 2
    a = np.random.uniform(0, 1, (bSize, nRows, nCols, n_in))
    b = f(a)
    c = g(a)
    print b - c
    """

def TestResNet():
    rng = np.random.RandomState()
    n_in = 3
    n_out = 4
    nRows = 20
    nCols = 20
    m = T.tensor3('m')
    x = T.tensor4('x')
    net = ResNet(rng, input = x, n_in=n_in, halfWinSize=1, n_hiddens=[30, 35, 40, 45], n_repeats=[5, 1, 0, 2], activation=T.nnet.relu, mask=m)
    #net2 = ResNet(rng, input = net1.output, n_in=3, halfWinSize=1, n_hiddens=[4], n_repeats=[1], activation=T.tanh, mask=m)
    y = net.output
    f = theano.function([x, m], y, on_unused_input='warn')

    bSize = 2
    a = np.random.uniform(0, 1, (bSize, nRows, nCols, n_in))
    m_value = np.ones( (bSize, 2, nCols) )

    """
    b = f(a, m_value)
    print b
    print b.shape
    """

    paramL2 = net.paramL2

    loss = T.mean( (y)**2 )
    cost = loss + 0.001 * paramL2
    #params = net1.params + net2.params
    params = net.params 
    gparams = T.grad(cost, params)
    h = theano.function([x, m], gparams, on_unused_input='warn')
    gs = h(a, m_value)

    for g in gs:
	print g.shape
	#print g

    updates = [ (p, p - 0.03 * g) for p, g in zip(params, gparams) ]
    train = theano.function([x, m], [cost, loss, paramL2], updates=updates)

    for i in xrange(1500):
	c, los, l2 = train(a, m_value)
	print c, los, l2

if __name__ == "__main__":

    #TestConvLayers()
    #TestResNet()
    TestBatchNorm4D()
