import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
#import ContactUtils
#from MergePredictedContactMatrix import MergeAndSaveOneProtein

import getopt

def Usage():
    	print 'python EvaluateDistanceAccuracy.py predictedDistMatrix_PKL ground_truth_PKL [-s minSeqSep ] [-c distCutoff for prediction] [ -d distCutoff for native]'
	print '  This script evaluate distance bound accuracy for a protein with its predicted distance/orientation matrix file '
    	print '  predictedDistMatrixPKL: a predicted distance matrix file with name like XXX.predictedDistMatrix.pkl'
    	print '     This file contains a tuple of at least 4 items: name, primary sequence, predictedDistMatrix, predictedContactMatrix'
	print '  ground_truth_PKL: a native ground truth file with name like XXX.native.pkl (default) or .atomDistMatrix.pkl'
	print '  -s: specify the minimum sequence separation (default 12)'
	print '  -c: specify the max dist cutoff for prediction (default 15.0)'
	print '  -d: specify the max dist cutoff for native (default 15.0)'
	print '  This script outputs absolute error of distance prediction, relative error, precision, recall, F1, Ri and GDT'
	print '		Ri is the ratio of predicted distance (<=15) with absolute prediction error no more than i Angstrom where i=0.5, 1, 2, 4, 8'

def main(argv):

	minSeqSep = 12
	distCutoff4Pred = 15.0
	distCutoff4Native = 15.0

	if len(argv)<2:
		Usage()
		exit(1)

	predFile = argv[0]
	nativeFile = argv[1]

	if not os.path.isfile(predFile):
		print 'ERROR: the predicted dist/orientation matrix file does not exist: ', predFile
		exit(1)

	if not os.path.isfile(nativeFile):
		print 'ERROR: the ground truth file does not exist: ', nativeFile
		exit(1)

	try:
                opts, args = getopt.getopt(argv[2:], "s:c:d:", ["minSeqSep=", "distCutoff4Pred=", "distCutoff4Native="])
                #print opts, args
        except getopt.GetoptError as err:
                print err
                Usage()
                exit(1)


	for opt, arg in opts:
                if opt in ("-s", "--minSeqSep"):
			minSeqSep = np.int32(arg)
			if minSeqSep < 1: 
				print 'ERROR: minSeqSep shall be a positive integer'
				exit(1)

		elif opt in ("-c", "--distCutoff4Pred"):
			distCutoff4Pred = np.float32(arg)
			if distCutoff4Pred < 1.0 :
				print 'ERROR: dist cutoff for predicted is too small '
				exit(1)

		elif opt in ("-d", "--distCutoff4Native"):
			distCutoff4Native = np.float32(arg)
			if distCutoff4Native < 1.0 :
				print 'ERROR: dist cutoff for native is too small '
				exit(1)
		else:
			Usage()
			exit(1)

	if distCutoff4Pred > distCutoff4Native:
		print 'WARNING: dist cutoff for pred is larger than dist cutoff for native'
		
	with open(predFile, 'rb') as fh:
		pred = cPickle.load(fh)
	assert len(pred)>=4
	bounds = EstimateDistanceBounds(pred[2])

	with open(nativeFile, 'rb') as fh:
		native = cPickle.load(fh)
	if nativeFile.endswith('native.pkl'):
		native = native['atomDistMatrix']

        acc = DistanceUtils.EvaluateDistanceBoundAccuracy(bounds, native, minSeqSep=minSeqSep, distCutoff4Pred=distCutoff4Pred, distCutoff4Native=distCutoff4Native)

	target = os.path.basename(nativeFile)
	for apt, value in acc.iteritems():
		print target, apt, value

if __name__ == "__main__":
    	main(sys.argv[1:])
