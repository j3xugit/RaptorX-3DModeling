import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np
import getopt

import config
import PropertyUtils

def Usage():
    	print 'python EvaluatePropertyAccuracy.py prediction_PKL ground_truth_PKL'
	print '  This script evaluate property prediction accuracy for a single protein'
    	print '  prediction_PKL: a predicted property file with name like XXX.predictedProperty.pkl'
    	print '     This file contains a tuple of 4 items: name, primary sequence, pred4prob, pred'
	print '     Both pred4prob and pred are dict(). pred4prob contains probability for property prediction and pred has the final prediction value'
	print '  ground_truth_PKL: a native property file with name like XXX.nativeProperty.pkl'
	print '  This script will output absolute error, Q8 and Q3 '

def main(argv):
	if len(argv)<2:
		Usage()
		exit(1)

	predFile = argv[0]
	nativeFile = argv[1]

	if not os.path.isfile(predFile):
		print 'ERROR: invalid predicted property file', predFile
		exit(1)

	if not os.path.isfile(nativeFile):
		print 'ERROR: invalid native property file', nativeFile
		exit(1)

	with open(predFile, 'rb') as fh:
		pred = cPickle.load(fh)

        avgerrors = PropertyUtils.EvaluateSinglePropertyPrediction(pred[0], nativeFile)

	target = os.path.basename(nativeFile)
	for apt, value in avgerrors.iteritems():
		print target, apt, value

if __name__ == "__main__":
    	main(sys.argv[1:])
