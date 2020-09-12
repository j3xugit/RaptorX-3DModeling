import sys
import os
import numpy as np
import cPickle

from ContactUtils import TopAccuracy
from utilsNoT import str_display

if __name__ == "__main__":

	if len(sys.argv) < 3:
    		print 'python EvaluateContactAccuracyPKL.py predDistMatrixFile_pkl groudTruthFile_pkl [targetName]'
		print '	This script evaluates the accuracy of predicted contact matrix by comparing it with native contact matrix'
		print '	predDistMatrixFile: the predicted distance/orientation matrix file in cPickle format, usually ending with .predictedDistMatrix.pkl '
		print '	     this file shall be a tuple of at least 6 items: name, sequence, predictedDistProbMatrix, predictedContactProbMatrix, labelWeight and labelDistribution'
		print '	groundTruthFile: the ground truth distance/orientation matrix file in cPickle format, usually ending with .native.pkl'
    		exit(1)

	predFile = sys.argv[1]
	nativeFile = sys.argv[2]

	if len(sys.argv)>3:
		target = sys.argv[3]
	else:
		target = os.path.basename(nativeFile).split('.')[0]

	with open(predFile, 'rb') as fh:
		pred = cPickle.load(fh)
	predCbCbContactMatrix = pred[3]['CbCb']
	
	with open(nativeFile, 'rb') as fh:
		truth = cPickle.load(fh)
	nativeCbCbDistMatrix = truth['atomDistMatrix']['CbCb']

	accs = TopAccuracy(predCbCbContactMatrix, nativeCbCbDistMatrix)
	## accs[0]: precision of extra-long-range, long-range, medium-range, long+medium-range, short-range
	## accs[1]: recall of extra-long-range, long-range, medium-range, long+medium-range, short-range
	## accs[2]: F1 of extra-long-range, long-range, medium-range, long+medium-range, short-range
	## for each range, the performance of top L, L/2, L/5 and L/10 contacts is evaluted

	seqLen = predCbCbContactMatrix.shape[0]
	for metric, values in zip(['precision', 'recall', 'F1'], accs):
		#resultStr = ' '.join( [target, str(seqLen), metric] + [ str(a) for a in values ] )
		resultStr = ' '.join( [target, str(seqLen), metric, str_display(values) ] )
		print resultStr

