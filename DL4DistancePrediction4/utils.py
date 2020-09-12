import numpy as np
import random
import theano
import theano.tensor as T
import math
import string

## this script contains some helper functions that depend on some theano modules

## d is a list of np.darray. This function converts each element in d from any float to floatX
def ToFloatX(d):
	newD = []
	for e in d:
		if e.dtype == theano.config.floatX:
			newD.append(e)
		elif e.dtype==np.float16 or e.dtype==np.float32 or e.dtype==np.float64:
			newD.append(e.astype(theano.config.floatX) )
		else:
			newD.append(e)
	return newD

##input has a shape (batchSize, seqLen, n_in)
##this function returns a output with shape (batchSize, seqLen, seqLen, 3*n_in) such as
## output[:,i ,j, :] = input[:, (i+j)/2, :], input[:, i,:], input[:, j,:]
## box is a bounding box (1-d ivector). When specified, only return a subtenor with shape (batchSize, top:bottom, left:right, 3*n_in) where box = [top, left, bottom, right]
def MidpointFeature(input, n_in, box=None):

    	if box is not None:
		top = box[0]
		left = box[1]
		bottom = box[2]
		right = box[3]
	else:
    		seqLen = input.shape[1]
		top = 0
		left = 0
		bottom = seqLen
		right = seqLen

	if box is None:
		x = T.mgrid[0:seqLen, 0:seqLen]
	else:
    		x = T.mgrid[top:bottom, left:right]
    	y1 = x[0] 
    	y2 = (x[0] + x[1])/2
    	y3 = x[1]

    	input2 = input.dimshuffle(1, 0, 2)

    	out1 = input2[y1]
    	out2 = input2[y2]
    	out3 = input2[y3]

    	out = T.concatenate([out1, out2, out3], axis=3)
    	final_out = out.dimshuffle(2, 0, 1, 3)
    	n_out = 3 * n_in

    	return final_out, n_out


def TestMidpointFeature():
    x = T.tensor3('x')
    y = MidpointFeature(x)
    f= theano.function([x], y)
    a = np.random.uniform(0, 1, (3, 10, 2)).astype(theano.config.floatX)
    b,c  = f(a)
    print c
    #return
    print '**********0*********'
    print a[0]
    print b[0][0]
    print '********4*******'
    print a[0]
    print b[0][4]
    print '**********9******'
    print a[0]
    print b[0][9]


def OuterConcatenate(input):
    ##an operation similar to the outer product of two vectors, but here we do concatenation instead of product of one vector with itself
    ##input has a shape (batchSize, seqLen, n_in), output has shape (batchSize, seqLen, seqLen, 2*n_in)

    seqLen = input.shape[1]
    input2 = input.dimshuffle(1, 0, 2)
    x = T.mgrid[0:seqLen, 0:seqLen]
    out = input2[x]
    output = T.concatenate((out[0], out[1]), axis=3)

    return output.dimshuffle(2, 0, 1, 3)

def TestOuterConcatenate():

    x = T.tensor3('x')
    f = theano.function([x], OuterConcatenate(x))

    a = np.random.uniform(0, 1, (2, 10, 2)).astype(np.float32)
    print a[0]
    print f(a)[0]

## Replace each vector (the last dim) of x by the linear combination of this vector and patterns 
## x is tensor4 with shape (bSize, nRows, nCols, n_in) where n_in = pshape[0]
## x can be interpreted as a set of reduced contact maps, where each contact map has size (nRows, nCols)
## the resultant matrix shall have shape (batchSize, nRows * pshape[1], nCols* pshape[2], 2), indicating the predicted prob of contacts and non-contacts
def ExpandBy3dPattern(x, patterns):

    ## patterns is a binary tensor3 with shape (numPatterns, nRows, nCols), where the 2nd and 3rd dims indicate the size of a pattern
    ## for each pattern element, 1 indicates a non-contact while  0 indicates a contact
    pshape = patterns.shape

    ## y1 has shape (batchSize, nRows * pshape[1], nCols * pshape[2], pshape[0])
    y1 = MyRepeat(x, (pshape[1], pshape[2]), axes=[1, 2])

    ##expand each pattern to a big matrix
    ## expandedPatterns has shape (1, nRows * pshape[1], nCols * pshape[2], pshape[0])
    #expandedPatterns = MyTile(patterns.dimshuffle(1, 2, 0), (x.shape[1], x.shape[2]) ).dimshuffle('x', 0, 1, 2)
    expandedPatterns = T.tile(patterns, (1, x.shape[1], x.shape[2]) ).dimshuffle('x', 1, 2, 0)


    ## calculate linear combination and the prob of non-contacts
    y2 = T.mul( y1, expandedPatterns)
    y3 = T.sum( y2, axis=3, keepdims=True)

    ##calculate the prob of contacts
    y4 = 1 - y3

    ##both y3 and y4 have shape (bSize, nRows*pshape[1], nCols*pshape[2])
    ## y3 contains the probability of non-contacts
    ## y4 contains the probability of contacts
    return T.concatenate([y4, y3], axis=3)

#this function returns a tensor with shape (bSize, nRows*nPatternRows, nCols*nPatternCols, numLabels)
def ExpandBy4dPattern(x, patterns):
    ##patterns has shape (numPatterns, nPatternRows, nPatternCols, numLabels)
    ##each element is between 0 and 1 and the sum of the vector patterns[i, j, k, :] is equal to 1
    pshape = patterns.shape

    ## y1 has shape (batchSize, nRows * pshape[1], nCols * pshape[2], pshape[0])
    y1 = MyRepeat(x, (pshape[1], pshape[2]), axes=[1, 2])
    expandedPatterns = T.tile(patterns, (1, x.shape[1], x.shape[2], 1) ).dimshuffle('x', 1, 2, 0, 3)

    ylist = []
    for i in xrange(pshape[3]):
	y2 = T.mul( y1, expandedPatterns[:, :, :, :, i] )
	y3 = T.sum( y2, axis=3, keepdims=True)
	ylist.append(y3)
    return T.concatenate( ylist, axis=3)

def ExpandByPattern(x, patterns):
    if patterns.ndim == 3:
	return ExpandBy3dPattern(x, patterns)
    elif patterns.ndim == 4:
	return ExpandBy4dPattern(x, patterns)
    else:
	print 'unsupported ndim of patterns: ', patterns.ndim
	exit(1)

## x has shape (bSize, nRows, nCols, numPatterns)
## patterns has shape (numPatterns, patternshape, numLabels)
def ConvByPattern(x, patterns, mask=None):
    W = np.transpose(patterns, (3, 0, 1, 2))
    out2 = T.nnet.conv2d(x.dimshuffle(0, 3, 1, 2), W, filter_shape=W.shape, border_mode='half')
    if mask is not None:
        ## mask has shape (batchSize, #rows_to_be_masked, nCols)

        ## a subtensor of out2 along the horiz direction
        out2_sub_horiz = out2[:, :, :mask.shape[1], :]
        mask_horiz = mask.dimshuffle(0, 'x', 1, 2)
        out3 = T.set_subtensor(out2_sub_horiz, T.mul(out2_sub_horiz, mask_horiz) )

        ## a subtensor of out3 along the vertical direction
        out3_sub_vertical = out3[:, :, :, :mask.shape[1] ]
        mask_vertical = mask.dimshuffle(0, 'x', 2, 1)
        y = T.set_subtensor(out3_sub_vertical, T.mul(out3_sub_vertical, mask_vertical) )
    else:
	y = out2

    y = y.dimshuffle(0, 2, 3, 1)

    return y/np.prod(patterns.shape[1:3])

  
def TestConvByPattern():
    x = T.tensor4('x')
    
    numPatterns = 10
    psize1 = 3
    psize2 = 3
    numLabels = 3
    bSize = 2
    nRows = 20
    nCols = 20

    pshape = (numPatterns, psize1, psize2, numLabels)
    patterns = np.random.uniform(0, 1, pshape).astype(theano.config.floatX)
    psum = np.sum(patterns, axis=3, keepdims=True)
    patterns = patterns / psum

    f = theano.function([x], ConvByPattern(x, patterns) )

    xshape = (bSize, nRows, nCols, numPatterns)
    a = np.random.uniform(0, 1, xshape).astype(theano.config.floatX)
    asum = np.sum(a, axis=3, keepdims=True)
    a = a / asum

    b = f(a)
    print b.shape
    print b[0][1][5]
    print b[1][3][3]

    print np.sum(b, axis=3)
    
def TestExpandByPattern():
    x = T.tensor4('x')
    pa = np.array([ [ 0, 0, 0], [ 0, 0, 1], [ 0, 1, 0] ]).astype(np.float32)
    pb = np.array([ [ 0, 1, 0], [ 0, 1, 1], [ 0, 1, 1] ]).astype(np.float32)
    pc = np.array([ [ 1, 1, 0], [ 1, 1, 1], [ 1, 0, 1] ]).astype(np.float32)
    pd = np.array([ [ 1, 1, 0], [ 1, 1, 1], [ 1, 0, 1] ]).astype(np.float32)
    pe = np.array([ [ 1, 1, 0], [ 1, 1, 1], [ 1, 0, 1] ]).astype(np.float32)

    patterns = np.array([ pa, pb, pc, pd, pe])

    pshape = patterns.shape

    ashape = (2, 4, 4, pshape[0])

    f = theano.function([x], ExpandByPattern(x, patterns) )

    a = np.random.uniform(0, 1, ashape).astype(theano.config.floatX)

    ##normalize a so that the sum of its last dim is 1
    b = np.sum(a, axis=3, keepdims=True)
    a = a/b

    b = f(a)
    print b
    print b.shape

    patterns = np.random.uniform(0, 1, (5, 3, 3, 3)).astype(np.float32)
    p2 = np.sum(patterns, axis=3, keepdims=True)
    p3 = patterns / p2
    f = theano.function([x], ExpandByPattern(x, p3) )
    b = f(a)
    print b
    print b.shape

def TestStack():
    x = T.matrix('x')
    y = T.matrix('y')
    z = T.matrix('z')

    f = theano.function([x, y, z], T.stacklists([ x, y, z]) )

    a = np.ones((5,4), dtype=np.float32)
    b = np.ones((5,4), dtype=np.float32)
    c = np.ones((5,4), dtype=np.float32)

    d = f(a, b, c)
    print d.shape
    print d


if __name__ == "__main__":
    #TestMyRepeat()
    #TestStack()
    #TestMyTile()
    #TestExpandByPattern()
    TestConvByPattern()

