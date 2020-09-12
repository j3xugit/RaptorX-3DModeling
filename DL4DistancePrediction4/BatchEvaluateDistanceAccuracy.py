import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
from EstimateAtomDistBounds import EstimateDistanceBounds

from utilsNoT import str_display

import getopt

def Usage():
    	print 'python BatchEvaluateDistanceAccuracy.py poteinList predDistMatrix_PKL_folder ground_truth_folder [ -s minSeqSep]  [-c distCutoff for prediction] [ -d distCutoff for native]'
	print '  This script evaluate distance bound accuracy for a list of proteins in their predicted distance matrix files '
    	print '  predDistMatrix_PKL_folder: a folder containing predicted distance/orientation matrix with name like XXX.predictedDistMatrix.pkl'
    	print '     This file contains a tuple of at least 3 items: name, primary sequence, predDistMatrix'
	print '  -s: optional. The minimum sequence separation between two residues for which its distance is evaluated. default 12'
	print '  -c: specify the max dist cutoff for prediction (default 15.0)'
        print '  -d: specify the max dist cutoff for native (default 15.0)'
	print '  This script will output absolute error of distance prediction, relative error, precision, recall, F1, Ri and GDT'
	print '		Ri is the ratio of predicted distance (<=15) with absolute prediction error no more than i Angstrom, i=0.5, 1, 2, 4, 8'

def main(argv):
	if len(argv)<3:
		Usage()
		exit(1)

	proteinListFile = argv[0]
	predFolder = argv[1]
	nativefolder = argv[2]
	fileSuffix = '.predictedDistMatrix.pkl'

	minSeqSep = 12
	distCutoff4Pred = 15.0
        distCutoff4Native = 15.0

	try:
                opts, args = getopt.getopt(argv[3:], "s:c:d:", ["minSeqSep=", "distCutoff4Pred=", "distCutoff4Native="])
                #print opts, args
        except getopt.GetoptError as err:
                print err
                Usage()
                exit(1)

        for opt, arg in opts:
                if opt in ("-s", "--minSeqSep"):
                        minSeqSep = np.int32(arg)
                        if minSeqSep < 2:
                                print 'ERROR: minSeqSep too small'
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

	if not os.path.isfile(proteinListFile):
		print 'ERROR: the protein list file does not exist: ', proteinListFile
		exit(1)

	if not os.path.isdir(predFolder):
		print 'ERROR: the folder for predicted distance matrix files does not exist: ', predFolder
		exit(1)

	if not os.path.isdir(nativefolder):
		print 'ERROR: the folder for native distance matrix files does not exist: ', nativefolder
		exit(1)

	with open(proteinListFile, 'r') as fh:
		proteins = [ line.strip() for line in list(fh) ]

	AccPerProtein = dict()
	accs = dict()
	for protein in proteins:
		predFile = os.path.join( predFolder, protein + fileSuffix ) 
		with open(predFile, 'rb') as fh:
			pred = cPickle.load(fh)
		assert len(pred)>=3
		predDistMatrix = pred[2]

		nativeFile = os.path.join(nativefolder, protein + '.native.pkl' )
		if not os.path.isfile(nativeFile):
			nativeFile = os.path.join(nativefolder, protein + '.atomDistMatrix.pkl' )

		with open(nativeFile, 'rb') as fh:
			native = cPickle.load(fh)
		if nativeFile.endswith('native.pkl'):
			native = native['atomDistMatrix']

		bounds = EstimateDistanceBounds(predDistMatrix)
                acc = DistanceUtils.EvaluateDistanceBoundAccuracy(bounds, native, minSeqSep=minSeqSep, distCutoff4Pred=distCutoff4Pred, distCutoff4Native=distCutoff4Native)
		AccPerProtein[protein] = acc

		for k, v in acc.iteritems():
			if not accs.has_key(k):
				accs[k] = [ v ]
			else:
				accs[k].append(v)

	for k, v in accs.iteritems():
		accs[k] = np.average(v, axis=0)
		print 'average', k, str_display(accs[k])

	for k, v in AccPerProtein.iteritems():
		for apt, value in v.iteritems():
			print k, apt, str_display(value)

if __name__ == "__main__":
    	main(sys.argv[1:])
