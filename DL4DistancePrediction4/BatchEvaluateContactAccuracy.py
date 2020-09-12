import cPickle
import sys
import os
#import scipy.stats.mstats
import numpy as np

import ContactUtils
from utilsNoT import str_display

#import config
#import DistanceUtils
#from MergePredictedContactMatrix import MergeAndSaveOneProtein
#import getopt

def Usage():
    	print 'python BatchEvaluateContactAccuracy.py poteinListFile PKL_folder ground_truth_folder [fileSuffix]'
	print '  This script evaluates contact prediction accuracy for a list of proteins in their predicted dist or contact matrix files '
    	print '  PKL_folder: a folder containing predicted distance/orientation matrix files named after XXX.predictedDistMatrix.pkl'
    	print '     A predicted distance matrix file contains a tuple of at least 6 items: name, sequence, predicted distance prob matrix, predicted contact prob matrix, labelWeights, labelDistribution'
	print '  ground_truth_folder: folder for native distance/orientation matrix files'
    	print '  fileSuffix: suffix for native dist/contact matrix file: .atomDistMatrix.pkl or .native.pkl (default)'
	print '		two different suffices represent two different formats'

def main(argv):

	if len(argv)<3:
		Usage()
		exit(1)

	proteinListFile = argv[0]
	predFolder = argv[1]
	nativefolder = argv[2]

	fileSuffix = ".native.pkl"
	if len(argv)>=4:
		fileSuffix = argv[3]

	if not os.path.isfile(proteinListFile):
		print 'ERROR: the protein list file does not exist: ', proteinListFile
		exit(1)

	if not os.path.isdir(predFolder):
		print 'ERROR: the folder for predicted matrix files does not exist: ', predFolder
		exit(1)

	if not os.path.isdir(nativefolder):
		print 'ERROR: the folder for native distance/orientation matrix files does not exist: ', nativefolder
		exit(1)
        print 'nativeFolder=', nativefolder

	with open(proteinListFile, 'r') as fh:
		proteins = [ line.strip() for line in list(fh) ]

	## load predictions
	predictions = dict()
	#allnatives = dict()
	for protein in proteins:
		predFile = os.path.join(predFolder, protein + '.predictedDistMatrix.pkl' ) 
		if not os.path.isfile(predFile):
			print 'ERROR: the prediction file does not exist: ', predFile
			exit(1)

		with open(predFile, 'rb') as fh:
			pred = cPickle.load(fh)
		assert len(pred) >= 6
		predContactMatrix = pred[3]
		predictions[ protein ] = predContactMatrix

        avgacc, allacc = ContactUtils.EvaluateContactPredictions(predictions, nativefolder, fileSuffix=fileSuffix)
	print '******************average and detailed contact prediction accuracy*********************'
	ContactUtils.PrintAllContactAccuracy(avgacc, allacc)

if __name__ == "__main__":
    	main(sys.argv[1:])
