import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np

import config
import DistanceUtils
import ContactUtils
from MergePredictedContactMatrix import MergeAndSaveOneProtein

import getopt

def Usage():

    	print 'python BatchMergePredictedContactMatrix.py [-m method | -s fileSuffix | -g ground_truth_folder | -c ] poteinList PKL_folders'
	print '  This script builds contact prediction matrices for a list of proteins from their predicted contact or distance matrix files '
	print '  -m: algorithm for merge, amean for arithmetic mean (default) and gmean for geometric mean'
    	print '  -s: suffix for the predicted dist matrix file, e.g., .predictedDistMatrix.pkl (default) or .fixedDistMatrix.pkl'
	print '  -g: specify the folder containing the ground truth. If provided, contact prediction accuracy will be calculated'
    	print '  -c: output the contact prob matrix in txt file. The Cb-Cb contact file has suffix .gcnn. The Ca-Ca contact file has suffix .CaCa.gcnn .'
    	print '  PKL_folders: a list of folders containing predicted distance matrix files with name like XXX.predictedDistMatrix.pkl or XXX.fixedDistMatrix.pkl'
    	print '     An input distance matrix file contains a tuple of at least 6 items: name, primary sequence, predicted distance prob matrix, predicted contact prob matrix, labelWeights, reference probabilities'

	print '  This script will output one or a few files. The resultant .predictedContactMatrix.pkl file is a python dictionary containing protein name, sequence and predicted contact matrix '
        print '  If the -c option is used, the resultant .gcnn file is a contact matrix in text format.'

def main(argv):

	methodpool = set(['amean', 'gmean'])

    	inputFiles = None
    	targetName = None
	method = 'amean'

	nativefolder = None
    	printContactMatrix = False

    	try:
        	opts, args = getopt.getopt(argv,"cm:g:s:",["contact=", "method=", "nativefolder=", "fileSuffix="])
        	print opts, args

    	except getopt.GetoptError:
        	Usage()
        	exit(1)

	if len(args)<3:
		Usage()
		exit(1)

	proteinListFile = args[0]
	#fileSuffix = args[1]
	inputFolders = args[1:]

	fileSuffix = ".predictedDistMatrix.pkl"

    	for opt, arg in opts:
		if opt in ("-c", "--contact"):
			printContactMatrix = True

		elif opt in ("-s", "--fileSuffix"):
			fileSuffix = arg

		elif opt in ("-g", "--nativefolder"):
			nativefolder = arg
			if not os.path.isdir(nativefolder):
                        	print 'The specified folder does not exist or is not accessible:', nativefolder
                                exit(1)

		elif opt in ("-m", "--method"):
			method = arg.strip().lower()

			if method not in methodpool:
				print 'ERROR: please specify a correct method for merge. It can only be amean or gmean.'
				exit(1)

		else:
	    		print Usage()
	    		exit(1)

	fh = open(proteinListFile, 'r')
	proteins = [ line.strip() for line in list(fh) ]
	fh.close()

	results = dict()
	for protein in proteins:
		inputFiles = [ os.path.join( folder, protein + fileSuffix ) for folder in inputFolders ]
		##print inputFiles
		result = MergeAndSaveOneProtein(inputFiles, method, printContactMatrix)
		results[ protein ] = result['predContactMatrix']

	if nativefolder is not None:
                print 'nativeFolder=', nativefolder
                avgacc, allacc = ContactUtils.EvaluateContactPredictions(results, nativefolder)
		ContactUtils.PrintAllContactAccuracy(avgacc, allacc)


if __name__ == "__main__":
    	main(sys.argv[1:])
