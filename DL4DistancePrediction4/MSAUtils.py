import os
import sys
import numpy as np
import cPickle
from ctypes import *
import datetime

lib4MSA = None
MSALibName = 'CLib4MSA.so'

simThreshold = 0.62
eps = np.finfo(np.float32).eps

# read a2m from *.a2m file
def LoadA2MFile(infile):
	with open(infile, 'r') as fh:
		seqs = [ s.strip() for s in list(fh) ]

	if len(seqs) < 1:
		print 'there is no sequence in a2m file: ', infile
		exit(1)

	## check to make sure all sequences have the same length here
	check = [ len(seqs[0]) == len(s) for s in seqs ]
	if not all(check):
		print 'inconsistent sequence length in a2m file: ', infile
		exit(1)

	return seqs

## read a2m and weight from *.a2m.pkl file
def LoadA2MPKLFile(infile):
	with open(infile, 'rb') as fh:
		data = cPickle.load(fh)

	seqs = data['seqs']
	weight = data['weight']

	return seqs, weight

# this function shall not be used now
## a2m is a list of strings
#3 this function converts it to a binary matrix of numSeqs * (seqLen*21)
def MSA2Matrix(sequences):
	m = []
	for seq in sequences:
		m.append( SeqOneHotEncodingWithGaps(seq).flatten() )

	return np.array(m, dtype=np.int16)

def CalcRealLengths(sequences):
	realLens = []
	for seq in sequences:
		realLen = len( seq.replace('-', '') )
		realLens.append(realLen)
	return realLens

##python code for weight calculation, which is slow
def CalcSeqWeight(sequences):
	numSeqs = len(sequences)
	realLens = CalcRealLengths(sequences)
	#print realLens
	requiredIDs =  [ simThreshold * rLen for rLen in realLens ]

	simMatrix = np.zeros( (numSeqs, numSeqs), np.int8)
	np.fill_diagonal(simMatrix, 1)

	for s1, i, l1 in zip(sequences, range(numSeqs), requiredIDs):
		for s2, j, l2 in zip(sequences, range(numSeqs), requiredIDs):
			if i>=j: 
				continue

			numIDs = [ a==b and a!='-' for a, b in zip(s1, s2) ].count(True)
			if numIDs > l1:
				simMatrix[i, j] = 1

			if numIDs > l2:
				simMatrix[j, i] = 1

	simScore = np.sum(simMatrix, axis=1).astype(np.int32)
	#print 'simScore: ', simScore

	weights = (1./simScore).astype(np.float32)
	#print weights

	return weights

## external C code for weight calculation
def CCalcSeqWeight(sequences):
	global lib4MSA
	if lib4MSA is None:
		libfile = os.path.join( os.path.dirname(os.path.abspath(__file__)), MSALibName)
		lib4MSA = cdll.LoadLibrary(libfile)

	CalcSeqWeight = lib4MSA.CalcSeqWeight
	CalcSeqWeight.argtypes = [ POINTER(c_char_p), c_int, c_int, POINTER(c_float) ]

	numSeqs = len(sequences)
	seqLen = len(sequences[0])

	weights = (c_float * numSeqs)()
	
	seqs = [ c_char_p(s) for s in sequences ]
	seqs = (c_char_p * numSeqs)(*seqs)
	
	CalcSeqWeight(seqs, numSeqs, seqLen, weights)

	return np.ctypeslib.as_array(weights, shape=(numSeqs,))

## generate a a2m.pkl file from an a2m file. the PKL file is a dict with two keys: seqs and weight
def GenA2MPKLFile(a2mfile, savefolder='./'):
	seqs = LoadA2MFile(a2mfile)
	weight = CCalcSeqWeight(seqs).astype(np.float16)
	result = dict()
	result['name'] = os.path.basename(a2mfile).split('.')[0]
	result['seqs'] = seqs
	result['weight'] = weight
	savefile = os.path.join(savefolder, os.path.basename(a2mfile) + '.pkl')
	with open(savefile, 'wb') as fh:
		cPickle.dump(result, fh, protocol=cPickle.HIGHEST_PROTOCOL)

## external C code for the calcuation of covariance/mutual information matrix
def CCalcPairwiseMatrix(sequences, weight=None, box=None, matrixType='covariance'):

	global lib4MSA
	if lib4MSA is None:
		libfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), MSALibName)
		lib4MSA = cdll.LoadLibrary(libfile)

	CalcPairwiseMatrix = lib4MSA.CalcPairMatrixFromMSA

	## it shall has the following arguments: seqs, weight, numSeqs, seqLen, bounds, flag, pairwiseMatrix, simpleMatrix
	CalcPairwiseMatrix.argtypes = [ POINTER(c_char_p), POINTER(c_float), c_int, c_int, POINTER(c_int), c_int, POINTER(c_float), POINTER(c_float)]

	numSeqs = len(sequences)
	seqLen = len(sequences[0])

	seqs = [ c_char_p(s) for s in sequences ]
	seqs = (c_char_p * numSeqs)(*seqs)

	if weight is None:
		weights = (c_float*numSeqs)()
	else:
		weights = (c_float*numSeqs)(*weight)

	## bounds = [ top, left, bottom, right]
	if box is None:
		bounds = (c_int*4)(0, 0, seqLen, seqLen)
	else:
		bounds = (c_int*4)(*box)

	flag = 1
	if matrixType != 'covariance':
		flag = 0

	blockSize = seqLen* seqLen * 441
	pairMatrix = (c_float * blockSize)()

	simpleMatrix = (c_float * (seqLen**2) )()

	CalcPairwiseMatrix(seqs, weights, numSeqs, seqLen, bounds, flag, pairMatrix, simpleMatrix)

	## add code here to convert the result back to ndarray
	resultMatrix = np.ctypeslib.as_array(pairMatrix, shape=(blockSize,)).reshape((bounds[2]-bounds[0], bounds[3]-bounds[1], 441))

	return resultMatrix

## python code for covariance calculation. This code shall not be used
def CalcCovMatrix(sequences, seqWeights):
	numSeqs = len(sequences)
	print '#seqs in MSA: ', numSeqs
	assert numSeqs > 0

	print 'Calc numeric represenation of MSA at ', datetime.datetime.now()

	seqLen = len( sequences[0] )
	seqMatrix = MSA2Matrix(sequences)
	assert seqMatrix.shape == (numSeqs, seqLen*21)

	print 'Calc covariance matrix at', datetime.datetime.now()
	cov = np.cov(seqMatrix.astype(np.float32), rowvar=False, aweights=seqWeights)
	assert cov.shape == (seqLen*21 , seqLen*21)

	cov = cov.reshape( (seqLen, 21, seqLen, 21) )
	cov = np.transpose(cov, (0, 2, 1, 3) )
	cov = cov.reshape( (seqLen, seqLen, 21*21) )

	print 'Finish calc cov matrix at ', datetime.datetime.now()

	return cov

def CalcPairMatrixFromFile(infile, box=None, matrixType='covariance'):

	basename = os.path.basename(infile)
	weight = None
	if basename.endswith('.a2m.pkl'):
		seqs, weight = LoadA2MPKLFile(infile)
	else:
		seqs = LoadA2MFile(infile)
		#weight = [1.] * len(seqs)

	pmatrix = CCalcPairwiseMatrix(seqs, weight=weight, matrixType=matrixType)

	return pmatrix

def Usage():
	print 'command: a2mFile [ResultFile] '

if __name__ == '__main__':
	if len(sys.argv) < 2:
		Usage()
		exit(1)

	starttime = datetime.datetime.now()
	a2mFile = sys.argv[1]
	cov = CalcPairMatrixFromFile(a2mFile)
	#mi = CalcPairMatrixFromFile(a2mFile, matrixType='mi').astype(np.float16)

	print 'time spent for calc pairwise matrix from file: ', datetime.datetime.now() - starttime

	if len(sys.argv) > 2:
		resultFile = sys.argv[2]
	else:
		resultFile = os.path.basename(a2mFile).split('.')[0] + '.covmi.pkl'

	with open(resultFile, 'wb') as fh:
		#cPickle.dump((cov,mi), fh, protocol=cPickle.HIGHEST_PROTOCOL)
		cPickle.dump(cov, fh, protocol=cPickle.HIGHEST_PROTOCOL)
