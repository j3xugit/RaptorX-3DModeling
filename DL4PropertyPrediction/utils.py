import numpy as np
import random
import theano
import theano.tensor as T
import math

def LoadFASTAFile(seqFile):

        ##load in the seq file
        seqfh = open(seqFile, 'r')
        content = [ line.strip() for line in list(seqfh) ]
        seqfh.close()
        content2 = [ c for c in content if not c.startswith('>') and not c.startswith('#') ]
        sequence = ''.join(content2)

        return sequence


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

##pretty print a 2D array or array_like objects such as list and tuple
def PrettyPrint(count):
	for row in count:
		line = ' '.join( '%.6f' % v for v in row)
		print line


##calculate the row-wise outer product of two matrix
## A has shape n*m, B has shape n*l, the result has shape n*(m*l)
def RowWiseOuterProduct(A, B):
    a = A[:, :, np.newaxis ]
    b = B[:, np.newaxis, : ]
    c = a * b
    return c.reshape( (A.shape[0], A.shape[1] * B.shape[1]) )

##check if two lists have the same shape
## list1 contains a set of theano shared variables
## list2 contains a set of np.ndarray
def Compatible(list1, list2):
    if len(list1) != len(list2):
	return False

    for l1, l2 in zip(list1, list2):
	if type(l1.get_value()) != type(l2):
		return False

	if np.isscalar(l1.get_value()):
	    continue

	if l1.get_value().shape != l2.shape:
	    return False

    return True

##generate the tile of a small tensor x, the first 2 dims will be expanded
## x is a small matrix or tensor3 to be tiled, y is a tuple of 2 elements
## This function generates a tile of x by copying it y*y times
## The resultant matrix shall have dimension ( x.shape[0]*y, x.shape[1]*y), consisting of y*y copies of x

def MyTile(x, y):
    p = T.mgrid[ 0: x.shape[0] * y[0], 0: x.shape[1] * y[1] ]
    q0 = p[0] % x.shape[0]
    q1 = p[1] % x.shape[1]
    z = x[(q0, q1)]
    return z

def TestMyTile():
    x = T.matrix('x')
    y = T.iscalar('y')
    z = MyTile(x, (y, y))

    f=theano.function([x, y], z)
    a = np.array([ [11, 12, 13, 14], [21, 22, 23, 24], [31, 32, 33, 34] ] ).astype(theano.config.floatX)
    b = 3
    c = f(a, b)

    print a
    print c
    print c.shape

    x3d = T.tensor3('x3d')
    g = theano.function([x3d, y], MyTile(x3d, (y, y)) )
    a = np.random.uniform(0, 1, (3, 3, 2)).astype(theano.config.floatX)
    c = g(a, b)
    print c.shape
    print c[:,:,0]
    print c[:,:,1]


##repeat each element in x by reps at the axis given by axes
def MyRepeat(x, reps, axes):
    assert len(reps) == len(axes)
    y = x
    for r, a in zip(reps, axes):
	y = T.repeat(y, [r], axis=a)
    return y


def ExpandByPattern(x, patterns):
    if patterns.ndim == 3:
	return ExpandBy3dPattern(x, patterns)
    elif patterns.ndim == 4:
	return ExpandBy4dPattern(x, patterns)
    else:
	print 'unsupported ndim of patterns: ', patterns.ndim
	sys.exit(-1)

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

def TestMyRepeat():

    x=T.tensor4('x')
    y=MyRepeat(x, (2, 2), axes=[1, 2])
    f=theano.function([x], y)

    a=np.random.uniform(0,1,(2, 3, 3, 2)).astype(theano.config.floatX)
    b=f(a)

    print a
    print b

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

