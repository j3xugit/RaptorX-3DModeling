import os
import sys
import numpy as np
import cPickle
import datetime

from Common.SequenceUtils import SeqOneHotEncodingWithGaps

simThreshold = 0.62
eps = np.finfo(np.float32).eps

def LoadA2MFile(infile):
	with open(infile, 'r') as fh:
		seqs = [ s.strip() for s in list(fh) ]
	return seqs

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

def CalcSeqWeights(sequences):
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

			#numIDs = np.sum( [ a==b and a!='-' for a, b in zip(s1, s2) ] )
			numIDs = [ a==b and a!='-' for a, b in zip(s1, s2) ].count(True)
			#numIDs = 50
			if numIDs > l1:
				simMatrix[i, j] = 1

			if numIDs > l2:
				simMatrix[j, i] = 1

	simScore = np.sum(simMatrix, axis=1).astype(np.int32)
	#print 'simScore: ', simScore

	weights = (1./simScore).astype(np.float32)
	#print weights

	## normalize weights
	weights = weights / np.mean(weights)
	#print weights

	return weights

def CalcCovMatrix(a2mFile):
	sequences = LoadA2MFile(a2mFile)
	numSeqs = len(sequences)
	print '#seqs in MSA: ', numSeqs
	assert numSeqs > 0

	print 'Calc numeric represenation of MSA at ', datetime.datetime.now()

	seqLen = len( sequences[0] )
	seqMatrix = MSA2Matrix(sequences)
	assert seqMatrix.shape == (numSeqs, seqLen*21)

	print 'Calc sequence weights at ', datetime.datetime.now()

	seqWeights = CalcSeqWeights(sequences)

	print 'Calc covariance matrix at', datetime.datetime.now()
	cov = np.cov(seqMatrix.astype(np.float32), rowvar=False, aweights=seqWeights)
	assert cov.shape == (seqLen*21 , seqLen*21)

	## normalize cov to corrcoef
	selfcov = np.diagonal(cov)
	cov = cov / np.sqrt(np.outer(selfcov, selfcov) + eps) 
	
	cov = cov.reshape( (seqLen, 21, seqLen, 21) )
	cov = np.transpose(cov, (0, 2, 1, 3) )
	cov = cov.reshape( (seqLen, seqLen, 21*21) )

	print 'Finish calc cov matrix at ', datetime.datetime.now()

	return cov

def Usage():
	print 'command: a2mFile [ResultFile] '

if __name__ == '__main__':
	if len(sys.argv) < 2:
		Usage()
		exit(1)

	a2mFile = sys.argv[1]
	cov = CalcCovMatrix(a2mFile).astype(np.float16)

	if len(sys.argv) > 2:
		resultFile = sys.argv[2]
	else:
		resultFile = os.path.basename(a2mFile).split('.')[0] + '.cov.pkl'

	with open(resultFile, 'wb') as fh:
		cPickle.dump(cov, fh, protocol=cPickle.HIGHEST_PROTOCOL)
