import numpy as np
import random
import math
import string

## this utils script contains some helper functions without using any theano module

def GenRandomString(stringLength=10):
	"""Generate a random string of fixed length """
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(stringLength))

## split a list into sublists of equal length
def SplitList(alist, wanted_parts=1):
    	length = len(alist)
    	return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]

## pretty display of a scalar, a list, a tuple or a 1D array
def str_display(ls, precision=4):
	fmtStr = '{0:.' + str(precision) + 'f}'
        if not isinstance(ls, (list, tuple, np.ndarray)):
                str_ls = fmtStr.format(ls)
                return str_ls

        str_ls = [fmtStr.format(v) for v in ls ]
        str_ls2 = '[' + ' '.join(str_ls) + ']'
        return str_ls2

##pretty print a 2D array or array_like objects such as list and tuple
def PrettyPrint(count):
	for row in count:
		line = ' '.join( '%.6f' % v for v in row)
		print line

## test if a string s is a number or not
def IsNumber(s):
        try:
                a=float(s)
                return True
        except ValueError:
                return False

## calculate the outer product of A and B row-by-row
def RowWiseOuterProduct(A, B):
	a = A[:, :, np.newaxis ]
	b = B[:, np.newaxis, : ]
	c = a * b
	return c.reshape( (A.shape[0], A.shape[1] * B.shape[1]) )

##sample a bounding box from a matrix. Currently only square box is sampled.
##the bounding box is saved as a 1d vector [top, left, bottom, right]
def SampleBoundingBox(original_shape, sizelimit):

	assert original_shape[0] == original_shape[1]
	#sample a bounding box such that its size is <= sizelimit
	if np.prod(original_shape) <= sizelimit:
		return np.array([0, 0, original_shape[0], original_shape[1] ]).astype(np.int32)

	if (sizelimit < 20*np.array(original_shape) ).any():
		print 'sizelimit is too small or the orignal shape is too big!'
		exit(1)

	## due to the difficulty in dealing with 2d mask, currently the sampled bounding box has to be a square
    	nRows = int( math.floor( math.sqrt(sizelimit) ) )
    	assert nRows <= original_shape[0]
    	nCols = nRows
    	assert nCols <= original_shape[1]

    	top = random.randrange(0, original_shape[0] - nRows + 1)
    	bottom = top + nRows
    	#left = top
    	left = random.randrange(0, original_shape[1] - nCols + 1)
    	right = left + nCols
    	return np.array([top, left, bottom, right]).astype(np.int32)
	
    
##check if the content of two lists have the same shape
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

#if __name__ == "__main__":
