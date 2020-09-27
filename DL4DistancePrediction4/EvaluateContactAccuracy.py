import sys
import os
import numpy as np
import cPickle

from ContactUtils import TopAccuracy
from ContactUtils import LoadContactMatrix
from utilsNoT import str_display

if __name__ == "__main__":

	if len(sys.argv) < 3:
    		print 'python EvaluateContactAccuracy.py predMatrixFile groudTruthFile [targetName]'
		print '\tThis script evaluates the accuracy of predicted contact matrix by comparing it with ground truth matrix'
		print '\tpredMatrixFile: the predicted matrix file in text or cPickle format, determined by its suffix'
		print '\t\twhen the suffix is .txt or .ccmpred, it is a text matrix file with L rows and L columns where L is protein length and each entry shall be predicted confidence/probability of being a contact'
		print '\t\twhen the suffix is .predictedDistMatrix.pkl, it is a file in cPickle format containing a tuple of at least 6 items: name, sequence, predictedDistProbMatrix, predictedContactProbMatrix, labelWeight and labelDistribution'
		print '\tgroundTruthFile: the ground truth file in text or cPickle format, usually ending with .native.pkl'
		print'\t\twhen it is a text format, its content shall be distance matrix instead of contact matrix, i.e., a large value indicates a non-contact and -1 indicates an invalid entry'
    		exit(1)

	predFile = sys.argv[1]
	nativeFile = sys.argv[2]

	if len(sys.argv)>3:
		target = sys.argv[3]
	else:
		target = os.path.basename(nativeFile).split('.')[0]

	if predFile.endswith('.txt') or predFile.endswith('.ccmpred'):
		predCbCbContactMatrix = LoadContactMatrix(predFile)
	elif predFile.endswith('.pkl'):
		with open(predFile, 'rb') as fh:
			pred = cPickle.load(fh)
		predCbCbContactMatrix = pred[3]['CbCb']
	else:
		print 'ERROR: predFile shall end with .txt or .pkl'
		exit(1)

	if nativeFile.endswith('.txt'):
		nativeCbCbDistMatrix = LoadContactMatrix(nativeFile)
	elif nativeFile.endswith('.pkl'):
		with open(nativeFile, 'rb') as fh:
			truth = cPickle.load(fh)
		nativeCbCbDistMatrix = truth['atomDistMatrix']['CbCb']
	else:
		print 'ERROR: ground truth file shall end with .txt or .pkl'
		exit(1)

	accs = TopAccuracy(predCbCbContactMatrix, nativeCbCbDistMatrix)
	## accs[0]: precision of extra-long-range, long-range, medium-range, long+medium-range, short-range
	## accs[1]: recall of extra-long-range, long-range, medium-range, long+medium-range, short-range
	## accs[2]: F1 of extra-long-range, long-range, medium-range, long+medium-range, short-range
	## for each range, top L, L/2, L/5 and L/10 predicted contacts are evaluted

	seqLen = predCbCbContactMatrix.shape[0]
	for metric, values in zip(['precision', 'recall', 'F1'], accs):
		#resultStr = ' '.join( [target, str(seqLen), metric] + [ str(a) for a in values ] )
		resultStr = ' '.join( [target, str(seqLen), metric, str_display(values) ] )
		print resultStr

