import cPickle
import sys
import os
import scipy.stats.mstats
import numpy as np
import getopt

import config
import DistanceUtils
import ContactUtils

from MergePredictedDistMatrix import MergeOneProtein

def Usage():

    	print 'python BatchMergePredictedDistMatrix.py [-m method | -s fileSuffix | -g ground_truth_folder | -d savefolder | -c ] proteinListFile predicteDistMatrix_folder1 predictedDistMatrix_folder2 ... '
	print '  -m: specify merge method: amean for arithmetic mean (default) while gmean for geometric mean'
	print '  -s: suffix for a predicted dist matrix file, e.g. .predictedDistMatrix.pkl (default) '
	print '  -g: the folder for ground truth. If provided, contact prediction accuracy will be calculated'
    	print '  -c: output predicted Cb-Cb contact matrix in text file, which shall have suffix .CM.txt and .CASP.rr'
	print '  proteinListFile: a file contains a list of protein names, each in one line'
	print '  predictedDistMatrix_folderXX: folders for predicted dist matrix'
    	print '  	A predicted dist matrix file has a name like XXX.predictedDistMatrix.pkl, e.g., 1f20A.predictedDistMatrix.pkl'
    	print '     	This file shall contain a tuple of at least 6 items: name, primary sequence, predicted distance prob, predicted contact prob, labelWeights and labelDistributions'

def main(argv):

	savefolder = os.getcwd()
    	printContactMatrix = False

    	inputFiles = None
    	targetName = None
	suffix = '.predictedDistMatrix.pkl'
	nativefolder = None

	method = 'amean'
	methodpool = set(['amean', 'gmean'])

    	try:
        	opts, args = getopt.getopt(argv,"cm:g:s:d:",["contact=", "method=", "nativefolder=", "fileSuffix=", "savefolder="])
        	print opts, args
    	except getopt.GetoptError:
        	Usage()
        	exit(1)

	if len(args) < 3:
		Usage()
		exit(1)

	proteinListFile = args[0]

	if not os.path.isfile(proteinListFile):
		print 'ERROR: invalid protein list file', proteinListFile
		exit(1)

	inputFolders = args[1:]
	for folder in inputFolders:
		if not os.path.isdir(folder):
			print 'ERROR: invalid folder for predicted dist matrices:', folder
			exit(1)

    	for opt, arg in opts:
		if opt in ("-c", "--contact"):
			printContactMatrix = True

		elif opt in("-s", "--fileSuffix"):
			suffix = arg

		elif opt in("-d", "--savefolder"):	
			savefolder = arg
			if not os.path.isdir(savefolder):
				os.mkdir(savefolder)

		elif opt in ("-g", "--nativefolder"):
			nativefolder = arg
			if not os.path.isdir(nativefolder):
                        	print 'ERROR: the specified native folder does not exist or is not accessible: ', nativefolder
                                exit(1)

		elif opt in ("-m", "--method"):
			method = arg.strip().lower()

			if method not in methodpool:
				print 'ERROR: please specify a correct method for merge. It can only be amean or gmean.'
				exit(1)

		else:
	    		Usage()
	    		exit(1)

	with open(proteinListFile, 'r') as fh:
		proteins = [ name.strip() for name in list(fh) ]

	contactPredictions = dict()
	for protein in proteins:
		inputFiles = []
		for folder in inputFolders:
			file = os.path.join(folder, protein + suffix )
			inputFiles.append(file)

		print 'Merging the predicted distance matrices for ', protein, '...'

		predictedContact, content4save = MergeOneProtein(inputFiles, method)
		contactPredictions[protein] = predictedContact 

		savefile = protein + suffix
		savefile = os.path.join(savefolder, savefile)
        	with open(savefile, 'wb') as fh:
        		cPickle.dump( content4save, fh, protocol = cPickle.HIGHEST_PROTOCOL)

    		if not printContactMatrix:
			continue

		for apt, m in predictedContact.iteritems():
			if apt == 'CbCb':
    				contactFileName = os.path.join(savefolder, targetName + '.CM.txt')
				contactCASPFileName = os.path.join(savefolder, targetName + '.CASP.rr')
			else:
				continue
    				contactFileName = targetName + '.' + apt + '.CM.txt'
				contactCASPFileName = targetName + '.' + apt + '.CASP.rr'

    			np.savetxt(contactFileName, m, fmt='%1.6f', delimiter=' ')
			SaveContactMatrixInCASPFormat(m, contactCASPFileName)

	if nativefolder is not None:
                print 'nativeFolder=', nativefolder
                avgAcc, allAccs = ContactUtils.EvaluateContactPredictions(contactPredictions, nativefolder)
		ContactUtils.PrintAllContactAccuracy(avgAcc, allAccs)

if __name__ == "__main__":
    	main(sys.argv[1:])
